import numpy as np
import pandas as pd
from sim.intv import Demography, ProcACF, ProcAltACF
from sim.dy import Model
import sim.dy.keys as I
from sim.intv import Intervention
from scipy.integrate import solve_ivp

__author__ = 'Chu-Chang Ku'
__all__ = ['Model']


class ModelIntv(Model):
    def __init__(self, year0=1970):
        Model.__init__(self, year0)
        self.Parent = Model(year0)
        self.Demography = Demography(I)
        self.ACF = ProcACF(I)
        self.AltACF = ProcAltACF(I)

    def __call__(self, t, y, pars, intv=None):
        y = y.reshape((I.N_State_TB, I.N_State_Strata * 2))

        dy = self.Demography.calc_dy(t, y, pars)
        dy += self.Transmission.calc_dy(t, y, pars)
        dy[:, :2] += self.Progression.calc_dy(t, y[:, :2], pars)
        dy[:, :2] += self.Cascade.calc_dy(t, y[:, :2], pars)

        if y[:, 2:].sum() > 0:
            dy[:, 2:] += self.Progression.calc_dy(t, y[:, 2:], pars)
            dy[:, 2:] += self.Cascade.calc_dy(t, y[:, 2:], pars)

        # ACF
        dy += self.ACF.calc_dy(t, y, (pars, intv))

        # Vulnerability-led ACF
        dy += self.AltACF.calc_dy(t, y, (pars, intv))

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
        self.ACF.measure(t, y, (pars, intv), mea)
        self.AltACF.measure(t, y, (pars, intv), mea)
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
        spec = p['acf_sym_spec'] if 'acf_sym_spec' in p else 1 - 0.036
        self.update_triage(y0, p, spec)
        return y0r, p

    @staticmethod
    def _find_r_acf0(y0, p, yield_mdu=430 / 3e6, yield_d2d=104.5 / 3e6):
        pos_sym, pos_cxr, pos_xpert, eligible = p['pos_sym'], p['pos_cxr'], p['pos_xpert'], p['eli']

        p['r_acf_mdu'] = yield_mdu / ((y0 * eligible * (1 - (1 - pos_sym) * (1 - pos_cxr)) * pos_xpert).sum() / y0.sum())
        p['r_acf_d2d'] = yield_d2d / ((y0 * eligible * pos_sym * pos_xpert).sum() / y0.sum())
        p['acf_sym_spec'] = p['acf_sym_spec']

    @staticmethod
    def update_triage(y0, p, spec_sym=1 - 0.036):
        sens_cxr, spec_cxr = p['acf_cxr_sens'], p['acf_cxr_spec']
        sens_xpert, spec_xpert = p['acf_xpert_sens'], p['acf_xpert_spec']

        eligible = np.ones((y0.shape[0], 1))
        eligible[I.Tx_DS] = 0
        eligible[I.Tx_DR] = 0
        eligible[I.LTBI_TPT] = 0

        pos_vul = np.ones((I.N_State_TB, 4))
        pos_vul[:, 0] = 0
        pos_vul[:, 0] = 0

        pos_cxr = np.ones((I.N_State_TB, 1))
        pos_cxr[I.U] = (1 - spec_cxr)
        pos_cxr[I.LTBI] = (1 - spec_cxr)
        pos_cxr[I.Infectious] = sens_cxr

        pos_xpert = np.ones((I.N_State_TB, 1))
        pos_xpert[I.U] = (1 - spec_xpert)
        pos_xpert[I.LTBI] = (1 - spec_xpert)
        pos_xpert[I.Infectious] = sens_xpert

        pos_sym = np.zeros((I.N_State_TB, 1))
        pos_sym[I.LTBI + I.Asym] = (1 - spec_sym)
        pos_sym[I.U] = (1 - spec_sym)
        pos_sym[I.Sym + I.ExSym] = 1

        p['alg:Sy'] = {
            'pos': pos_sym * pos_xpert,
            'use_sym': 1, 'use_vul': 0, 'use_vs': 0, 'use_cxr': 0,
            'use_xpert': pos_sym
        }

        p['alg:SyCx'] = {
            'pos': (1 - (1 - pos_sym) * (1 - pos_cxr)) * pos_xpert,
            'use_sym': 1, 'use_vul': 0, 'use_vs': 0, 'use_cxr': 1,
            'use_xpert': (1 - (1 - pos_sym) * (1 - pos_cxr))
        }

        pos_sym = pos_sym.reshape(-1)
        pos_cxr = pos_cxr.reshape(-1)

        pos = np.zeros((I.N_State_TB, 4))
        pos[:, 0] = pos_sym + (1 - pos_sym) * pos_vul[:, 0] * pos_cxr
        pos[:, 1] = pos_sym + (1 - pos_sym) * pos_vul[:, 1] * pos_cxr
        pos[:, 2] = pos_sym + (1 - pos_sym) * pos_cxr
        pos[:, 3] = pos_sym + (1 - pos_sym) * pos_cxr

        cxr = np.zeros((I.N_State_TB, 4))
        cxr[:, 0] = (1 - pos_sym) * pos_vul[:, 0]
        cxr[:, 1] = (1 - pos_sym) * pos_vul[:, 1]
        cxr[:, 2] = (1 - pos_sym)
        cxr[:, 3] = (1 - pos_sym)

        p['alg:VSC'] = {
            'pos': pos * pos_xpert,
            'use_sym': np.array([0, 0, 1, 1]).reshape((1, -1)),
            'use_vul': 0,
            'use_vs': np.array([1, 1, 0, 0]).reshape((1, -1)),
            'use_cxr': cxr,
            'use_xpert': pos_sym.reshape((-1, 1))
        }

        p['eli'] = eligible

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
    p0.update({'beta_ds': 25, 'rr_risk_comorb': 20, 'rr_beta_dr': 1.05, 'p_comorb': 0.15})

    y1, p1 = m0.find_baseline(p0, 2022)

    _, ms0, _ = m0.simulate_onward(y1, p1, intv={'FullACF': {'Coverage': 0}})
    # _, ms1, _ = m0.simulate_onward(y1, p1, intv={'D2D': {'Scale': 1}, 'MDU': {'Scale': 1}})
    # _, ms2, _ = m0.simulate_onward(y1, p1, intv={'D2D': {'Scale': 2}, 'MDU': {'Scale': 2}})
    # _, ms1, _ = m0.simulate_onward(y1, p1, intv={'FullACF': {'Coverage': 0.15, 'FollowUp': 3,
    #                                                          'Duration': 2, 'ScreenAlg': 'VSC'}})
    # _, ms2, _ = m0.simulate_onward(y1, p1, intv={'FullACF': {'Coverage': 0.15, 'FollowUp': 1,
    #                                                          'Duration': 2, 'ScreenAlg': 'VSC'}})

    _, ms1, _ = m0.simulate_onward(y1, p1, intv={'FullACF': {'Coverage': 0.15, 'ScreenAlg': 'Sy'}})
    _, ms2, _ = m0.simulate_onward(y1, p1, intv={'AltACF': {'Coverage': 0.15, 'ScreenAlg': 'Sy'}})

    # print('MDU', ms1.ACF_MDU_Yield[2022.5] * 1e5, 430 / 3e6 * 1e5)
    # print('D2D', ms1.ACF_D2D_Yield[2022.5] * 1e5, 104.5 / 3e6 * 1e5)

    print(ms1[['IncR', 'ACF_Yield']])

    fig, axes = plt.subplots(3, 2)

    ms0.Prev.plot(ax=axes[0, 0])
    ms1.Prev.plot(ax=axes[0, 0])
    ms2.Prev.plot(ax=axes[0, 0])

    axes[0, 0].legend(['I0', 'I1', 'I2'])
    axes[0, 0].set_title('Prevalence')

    ms0.IncR.plot(ax=axes[1, 0])
    ms1.IncR.plot(ax=axes[1, 0])
    ms2.IncR.plot(ax=axes[1, 0])
    print(ms1.IncR)
    print(ms2.IncR)

    axes[1, 0].legend(['I0', 'I1', 'I2'])
    axes[1, 0].set_title('Incidence')

    ms0.MorR.plot(ax=axes[2, 0])
    ms1.MorR.plot(ax=axes[2, 0])
    ms2.MorR.plot(ax=axes[2, 0])

    axes[2, 0].legend(['I0', 'I1', 'I2'])
    axes[2, 0].set_title('Mortality')

    ms0.PrOnTPT.plot(ax=axes[0, 1])
    ms1.PrOnTPT.plot(ax=axes[0, 1])
    ms2.PrOnTPT.plot(ax=axes[0, 1])

    axes[0, 1].legend(['I0', 'I1', 'I2'])
    axes[0, 1].set_title('On Pseudo TPT')

    ms0.ACF_Yield.plot(ax=axes[1, 1])
    ms1.ACF_Yield.plot(ax=axes[1, 1])
    ms2.ACF_Yield.plot(ax=axes[1, 1])

    axes[1, 1].legend(['I0', 'I1', 'I2'])
    axes[1, 1].set_title('Yield, Vul-led')

    ms0.ACF_fu_Yield.plot(ax=axes[2, 1])
    ms1.ACF_fu_Yield.plot(ax=axes[2, 1])
    ms2.ACF_fu_Yield.plot(ax=axes[2, 1])

    axes[2, 1].legend(['I0', 'I1', 'I2'])
    axes[2, 1].set_title('Yield, Follow-up')

    fig.tight_layout()

    plt.show()
