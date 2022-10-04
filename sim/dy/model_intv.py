import numpy as np
import pandas as pd
from sim.intv import Demography, BgACF, VulACF
from sim.dy import Model
import sim.dy.keys as I
from sim.intv import Intervention
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

__author__ = 'Chu-Chang Ku'
__all__ = ['Model']


class ModelIntv(Model):
    def __init__(self, year0=1970):
        Model.__init__(self, year0)
        self.Parent = Model(year0)
        self.Demography = Demography(I)
        self.BgACF = BgACF(I)
        self.VulACF = VulACF(I)

    def __call__(self, t, y, pars, intv=None):
        y = y.reshape((I.N_State_TB, I.N_State_Strata * 2))

        dy = self.Demography.calc_dy(t, y, pars)
        dy += self.Transmission.calc_dy(t, y, pars)
        dy[:, :2] += self.Progression.calc_dy(t, y[:, :2], pars)
        dy[:, :2] += self.Cascade.calc_dy(t, y[:, :2], pars)

        if y[:, 2:].sum() > 0:
            dy[:, 2:] += self.Progression.calc_dy(t, y[:, 2:], pars)
            dy[:, 2:] += self.Cascade.calc_dy(t, y[:, 2:], pars)

        # Baseline ACF
        dy += self.BgACF.calc_dy(t, y, pars, intv)

        # Vulnerability-led ACF
        dy += self.VulACF.calc_dy(t, y, pars, intv)

        return dy.reshape(-1)

    def measure(self, t, y, pars, intv=None):
        y = y.reshape((I.N_State_TB, I.N_State_Strata * 2))

        mea = {'Time': t}
        self.Demography.measure(t, y, pars, mea)

        # baseline metrics
        y_all = y[:, :2] + y[:, 2:]
        self.Transmission.measure(t, y_all, pars, mea)
        self.Progression.measure(t, y_all, pars, mea)
        self.Cascade.measure(t, y_all, pars, mea)

        # ACF metrics
        self.BgACF.measure(t, y, pars, intv, mea)
        self.VulACF.measure(t, y, pars, intv, mea)
        return mea

    @staticmethod
    def dfe(t, y, pars, intv=None):
        ntb = y.reshape((I.N_State_TB, I.N_State_Strata * 2))[I.Infectious].sum()
        return ntb - 0.5

    dfe.terminal = True
    dfe.direction = -1

    def find_baseline(self, p, t_start=2022):
        y0, _, _ = self.Parent.simulate_to_baseline(p, t_start)
        y0, p = self._bridge_y0(y0, p)
        return y0, p

    def _bridge_y0(self, y0, p):
        y0 = y0.reshape((-1, 2))
        y0r = np.zeros((y0.shape[0], 4))
        y0r[:, :2] = y0

        # Update triage
        spec = p['acf_sym_spec'] if 'acf_sym_spec' in p else None
        self._update_triage(y0, p, spec)

        # Baseline ACF
        if 'r_acf_mdu' not in p:
            self._find_r_acf0(y0, p)

        return y0r, p

    @staticmethod
    def _find_r_acf0(y0, p, yield_mdu=430 / 3e6, yield_d2d=104.5 / 3e6):
        pos_sym, pos_cxr, pos_xpert, eligible = p['pos_sym'], p['pos_cxr'], p['pos_xpert'], p['eli']

        p['r_acf_mdu'] = yield_mdu / ((y0 * eligible * pos_cxr * pos_xpert).sum() / y0.sum())
        p['r_acf_d2d'] = yield_d2d / ((y0 * eligible * pos_sym * pos_xpert).sum() / y0.sum())

    @staticmethod
    def _update_triage(y0, p, spec_sym=None):
        sens_cxr, spec_cxr = p['acf_cxr_sens'], p['acf_cxr_spec']
        sens_xpert, spec_xpert = p['acf_xpert_sens'], p['acf_xpert_spec']

        eligible = np.ones_like(y0)
        eligible[I.Tx_DS] = 0
        eligible[I.Tx_DR] = 0
        eligible[I.LTBI_TPT] = 0

        vul = eligible.copy()
        vul[:, 0] = 0

        pos_cxr = np.ones((I.N_State_TB, 1))
        pos_cxr[I.U] = (1 - spec_cxr)
        pos_cxr[I.LTBI] = (1 - spec_cxr)
        pos_cxr[I.Infectious] = sens_cxr

        pos_xpert = np.ones((I.N_State_TB, 1))
        pos_xpert[I.U] = (1 - spec_xpert)
        pos_xpert[I.LTBI] = (1 - spec_xpert)
        pos_xpert[I.Infectious] = sens_xpert

        if spec_sym is None:
            pos = (y0 * eligible * pos_cxr * pos_xpert)
            ppv_mu = pos[I.Infectious].sum() / pos.sum()

            def fn(spec, y):
                pos_sym = np.zeros_like(y)
                pos_sym[I.LTBI + I.Asym] = (1 - spec)
                pos_sym[I.U] = (1 - spec)
                pos_sym[I.Sym + I.ExSym] = 1

                pos = y * eligible * pos_sym * pos_xpert
                ppv_d2d = pos[I.Infectious].sum() / pos.sum()

                return (ppv_d2d - ppv_mu) ** 2

            opt = minimize_scalar(fn, 0.99, args=(y0, ), method='bounded', bounds=(0.5, 1))
            spec_sym = opt.x

        pos_sym = np.ones((I.N_State_TB, 1))
        pos_sym[I.LTBI + I.Asym] = (1 - spec_sym)
        pos_sym[I.U] = (1 - spec_sym)
        pos_sym[I.Sym + I.ExSym] = 1

        pos_vul = np.ones_like(y0)
        pos_vul[:, 1] = 1

        p.update({
            'eli': eligible,
            'eli_vul': vul,
            'pos_sym': pos_sym,
            'pos_vul': pos_vul,
            'pos_cxr': pos_cxr,
            'pos_xpert': pos_xpert,
        })

    def simulate_onward(self, y0, p, intv=None, t_start=2022, t_end=2031.0, dt=0.5):
        t_out = np.linspace(t_start, t_end, int((t_end - t_start) / dt) + 1)
        if intv is not None:
            intv = Intervention.parse_obj(intv)

        if 'sus' not in p:
            p = self.update_parameters(p)

        ys = solve_ivp(self, [t_start, t_end], y0.reshape(-1), args=(p, intv, ), events=self.dfe, dense_output=True)

        if len(ys.t_events[0]) > 0 or not ys.success:
            return None, None, {'succ': False, 'res': 'DFE reached'}

        ms = [self.measure(t, ys.sol(t), p, intv) for t in t_out if t >= t_start]
        ms = pd.DataFrame(ms).set_index('Time')
        return ys, ms, {'succ': True}


if __name__ == '__main__':
    import json
    import matplotlib.pylab as plt

    with open('../../data/test_priors.json', 'r') as f:
        prior = json.load(f)

    m0 = ModelIntv()

    p0 = prior[2]
    p0.update({'beta_ds': 5, 'rr_risk_comorb': 20, 'rr_beta_dr': 1.05, 'p_comorb': 0.3})

    y1, p1 = m0.find_baseline(p0, 2022)

    _, ms0, _ = m0.simulate_onward(y1, p1, intv={'D2D': {'Scale': 0}, 'MDU': {'Scale': 0}})
    # _, ms1, _ = m0.simulate_onward(y1, p1, intv={'D2D': {'Scale': 1}, 'MDU': {'Scale': 1}})
    # _, ms2, _ = m0.simulate_onward(y1, p1, intv={'D2D': {'Scale': 2}, 'MDU': {'Scale': 2}})
    _, ms1, _ = m0.simulate_onward(y1, p1, intv={'VulACF': {'Coverage': 0.15, 'FollowUp': 0.8, 'Duration': 2}})
    _, ms2, _ = m0.simulate_onward(y1, p1, intv={'VulACF': {'Coverage': 0.15, 'FollowUp': 1, 'Duration': 2}})
    # _, ms2, _ = m0.simulate_onward(y1, p1, intv={'PlainACF': {'Coverage': 0.15}})

    print('MDU', ms1.ACF_MDU_Yield[2022.5] * 1e5, 430 / 3e6 * 1e5)
    print('D2D', ms1.ACF_D2D_Yield[2022.5] * 1e5, 104.5 / 3e6 * 1e5)

    print(ms1[['IncR', 'ACF_Vul_Yield']])

    fig, axes = plt.subplots(3, 2)

    ms0.Prev.plot(ax=axes[0, 0])
    ms1.Prev.plot(ax=axes[0, 0])
    ms2.Prev.plot(ax=axes[0, 0])

    axes[0, 0].legend(['I0', 'I1', 'I2'])
    axes[0, 0].set_title('Prevalence')

    ms0.IncR.plot(ax=axes[1, 0])
    ms1.IncR.plot(ax=axes[1, 0])
    ms2.IncR.plot(ax=axes[1, 0])

    axes[1, 0].legend(['I0', 'I1', 'I2'])
    axes[1, 0].set_title('Incidence')

    ms0.MorR.plot(ax=axes[2, 0])
    ms1.MorR.plot(ax=axes[2, 0])
    ms2.MorR.plot(ax=axes[2, 0])

    axes[2, 0].legend(['I0', 'I1', 'I2'])
    axes[2, 0].set_title('Mortality')

    ms0.ACF_MDU_Yield.plot(ax=axes[0, 1])
    ms1.ACF_MDU_Yield.plot(ax=axes[0, 1])
    ms2.ACF_MDU_Yield.plot(ax=axes[0, 1])

    axes[0, 1].legend(['I0', 'I1', 'I2'])
    axes[0, 1].set_title('Yield, MDU')

    ms0.ACF_D2D_Yield.plot(ax=axes[1, 1])
    ms1.ACF_D2D_Yield.plot(ax=axes[1, 1])
    ms2.ACF_D2D_Yield.plot(ax=axes[1, 1])

    axes[1, 1].legend(['I0', 'I1', 'I2'])
    axes[1, 1].set_title('Yield, D2D')

    ms0.ACF_Vul_Yield.plot(ax=axes[2, 1])
    ms1.ACF_Vul_Yield.plot(ax=axes[2, 1])
    ms2.ACF_Vul_Yield.plot(ax=axes[2, 1])

    axes[2, 1].legend(['I0', 'I1', 'I2'])
    axes[2, 1].set_title('Yield, Vul')

    fig.tight_layout()

    plt.show()
