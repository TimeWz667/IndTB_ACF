import numpy as np
import pandas as pd
from sim.dy import Model
import sim.dy.keys as I
from sim.intv import Intervention
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

__author__ = 'Chu-Chang Ku'
__all__ = ['Model']


class ModelIntv(Model):
    def __init__(self, year0=1970, time_out=0):
        Model.__init__(self, year0)
        self.Parent = Model(year0)
        self.TimeOut = time_out

    def __call__(self, t, y, pars, intv=None):
        y = y.reshape((I.N_State_TB, I.N_State_Strata * 2))
        dy = np.zeros_like(y)

        dy[:, :2] = self.Parent(t, y[:, :2], pars).reshape((-1, 2))
        if y[:, 2:].sum() > 0:
            dy[:, 2:] = self.Parent(t, y[:, 2:], pars).reshape((-1, 2))

        # ACF
        pos_sym, pos_cxr, pos_xpert, eligible = pars['pos_sym'], pars['pos_cxr'], pars['pos_xpert'], pars['eligible']
        r_acf_mu, r_acf_d2d, p_dst = pars['r_acf_mu'], pars['r_acf_d2d'], pars['acf_dst_sens']
        if intv is not None:
            r_acf_mu, r_acf_d2d, p_dst = intv.modify_acf_bg(t, r_acf_mu, r_acf_d2d, p_dst)

        # intv
        acf_mu = r_acf_mu * eligible * pos_cxr * pos_xpert * y[:, :2]
        acf_mu_tp_DS, acf_mu_tp_DR = acf_mu[I.Infectious_DS], acf_mu[I.Infectious_DR]
        acf_d2d = r_acf_d2d * eligible * pos_sym * pos_xpert * y[:, :2]
        acf_d2d_tp_DS, acf_d2d_tp_DR = acf_d2d[I.Infectious_DS], acf_d2d[I.Infectious_DR]

        acf_tp_ds, acf_tp_dr = acf_mu_tp_DS + acf_d2d_tp_DS, acf_mu_tp_DR + acf_d2d_tp_DR

        dy[I.Infectious_DS, :2] -= acf_tp_ds
        dy[I.Infectious_DR, :2] -= acf_tp_dr
        dy[I.Txf_Pub_DS, :2] += acf_tp_ds.sum(0)
        dy[I.Txf_Pub_DR, :2] += acf_tp_dr.sum(0) * (1 - p_dst)
        dy[I.Txs_Pub_DR, :2] += acf_tp_dr.sum(0) * p_dst

        # Vulnerability-led ACF
        r_acf_vul = 0
        if intv is not None:
            r_acf_vul = intv.modify_acf_vul(t, r_acf_vul)

            acf_vul = r_acf_vul * eligible * pars['pos_vul'] * pos_xpert * y[:, :2]
            acf_vul_tp_ds, acf_vul_tp_dr = acf_vul[I.Infectious_DS], acf_vul[I.Infectious_DR]

            dy[I.Infectious_DS, :2] -= acf_vul_tp_ds
            dy[I.Infectious_DR, :2] -= acf_vul_tp_dr
            dy[I.Txf_Pub_DS, :2] += acf_vul_tp_ds.sum(0)
            dy[I.Txf_Pub_DR, :2] += acf_vul_tp_dr.sum(0) * (1 - p_dst)
            dy[I.Txs_Pub_DR, :2] += acf_vul_tp_dr.sum(0) * p_dst
            # todo timeout setup

        if self.TimeOut > 0:
            lost = dy[:, 2:] / self.TimeOut
            dy[:, 2:] -= lost
            dy[:, :2] += lost

        return dy.reshape(-1)

    def measure(self, t, y, pars, intv=None):
        y = y.reshape((I.N_State_TB, I.N_State_Strata * 2))

        # ACF metrics
        pos_sym, pos_cxr, pos_xpert, eligible = pars['pos_sym'], pars['pos_cxr'], pars['pos_xpert'], pars['eligible']
        pos_vul = pars['pos_vul']
        r_acf_mu, r_acf_d2d, p_dst = pars['r_acf_mu'], pars['r_acf_d2d'], pars['acf_dst_sens']
        r_acf_vul = 0
        if intv is not None:
            r_acf_mu, r_acf_d2d, p_dst = intv.modify_acf_bg(t, r_acf_mu, r_acf_d2d, p_dst)
            r_acf_vul = intv.modify_acf_vul(t, r_acf_vul)

        acf_mu_reach0 = r_acf_mu * eligible * y[:, :2]
        acf_mu_reach1 = acf_mu_reach0 * pos_cxr
        acf_mu_reach2 = acf_mu_reach1 * pos_xpert

        acf_d2d_reach0 = r_acf_d2d * eligible * y[:, :2]
        acf_d2d_reach1 = acf_d2d_reach0 * pos_sym
        acf_d2d_reach2 = acf_d2d_reach1 * pos_xpert

        acf_vul_reach0 = r_acf_vul * eligible * y[:, :2]
        acf_vul_reach1 = acf_vul_reach0 * pos_vul
        acf_vul_reach2 = acf_vul_reach1 * pos_xpert

        n = y.sum()

        mea_acf = {
            'Reach_ACF_MU1': acf_mu_reach0.sum() / n,
            'Reach_ACF_MU2': acf_mu_reach1.sum() / n,
            'Yield_ACF_MU': acf_mu_reach2.sum() / n,
            'Reach_ACF_D2D1': acf_d2d_reach0.sum() / n,
            'Reach_ACF_D2D2': acf_d2d_reach1.sum() / n,
            'Yield_ACF_D2D': acf_d2d_reach2.sum() / n,
            'Reach_ACF_Vul1': acf_vul_reach0.sum() / n,
            'Reach_ACF_Vul2': acf_vul_reach1.sum() / n,
            'Yield_ACF_Vul': acf_vul_reach2.sum() / n
        }

        # baseline metrics
        y = y[:, :2] + y[:, 2:]
        mea = Model.measure(self, t, y, pars)
        mea.update(mea_acf)
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
        if 'r_acf_mu' not in p:
            self._find_r_acf0(y0, p)

        return y0r, p

    @staticmethod
    def _find_r_acf0(y0, p, yield_mu=104.5 / 3e6, yield_d2d=430 / 3e6):
        pos_sym, pos_cxr, pos_xpert, eligible = p['pos_sym'], p['pos_cxr'], p['pos_xpert'], p['eligible']

        p['r_acf_mu'] = yield_mu / ((y0 * eligible * pos_cxr * pos_xpert).sum() / y0.sum())
        p['r_acf_d2d'] = yield_d2d / ((y0 * eligible * pos_sym * pos_xpert).sum() / y0.sum())

    @staticmethod
    def _update_triage(y0, p, spec_sym=None):
        sens_cxr, spec_cxr = p['acf_cxr_sens'], p['acf_cxr_spec']
        sens_xpert, spec_xpert = p['acf_xpert_sens'], p['acf_xpert_spec']

        eligible = np.ones_like(y0)
        eligible[I.Tx_DS] = 0
        eligible[I.Tx_DR] = 0

        pos_cxr = np.zeros_like(y0)
        pos_cxr[I.U] = (1 - spec_cxr)
        pos_cxr[I.LTBI] = (1 - spec_cxr)
        pos_cxr[I.Infectious] = sens_cxr

        pos_xpert = np.zeros_like(y0)
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

            opt = minimize_scalar(fn, 1, args=(y0, ), method='bounded', bounds=(0.5, 1))
            spec_sym = opt.x

        pos_sym = np.zeros_like(y0)
        pos_sym[I.LTBI + I.Asym] = (1 - spec_sym)
        pos_sym[I.U] = (1 - spec_sym)
        pos_sym[I.Sym + I.ExSym] = 1

        pos_vul = np.zeros_like(y0)
        pos_vul[:, 1] = 1

        p.update({
            'eligible': eligible,
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

    _, ms0, _ = m0.simulate_onward(y1, p1, intv={'D2D': {'Scale': 0}, 'MU': {'Scale': 0}})
    # _, ms1, _ = m0.simulate_onward(y1, p1, intv={'D2D': {'Scale': 1}, 'MU': {'Scale': 1}})
    # _, ms2, _ = m0.simulate_onward(y1, p1, intv={'D2D': {'Scale': 2}, 'MU': {'Scale': 2}})
    _, ms1, _ = m0.simulate_onward(y1, p1, intv={'D2D': {'Scale': 0}, 'MU': {'Scale': 0}, 'VulACF': {'Coverage': 0.2}})
    _, ms2, _ = m0.simulate_onward(y1, p1, intv={'D2D': {'Scale': 0}, 'MU': {'Scale': 0}, 'VulACF': {'Coverage': 0.5}})
    # _, ms2, _ = m.simulate_onward(y1, p0, intv={'ACF': {'Yield': .05, 'HiRisk': False}})
    # _, ms3, _ = m.simulate_onward(y1, p0, intv={'ACF': {'Yield': .03, 'HiRisk': True}})

    print('MU', ms1.Yield_ACF_MU[2022.5] * 1e5, 104.5 / 3e6 * 1e5)
    print('D2D', ms1.Yield_ACF_D2D[2022.5] * 1e5, 430 / 3e6 * 1e5)

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

    ms0.Yield_ACF_MU.plot(ax=axes[0, 1])
    ms1.Yield_ACF_MU.plot(ax=axes[0, 1])
    ms2.Yield_ACF_MU.plot(ax=axes[0, 1])

    axes[0, 1].legend(['I0', 'I1', 'I2'])
    axes[0, 1].set_title('Yield, MU')

    ms0.Yield_ACF_D2D.plot(ax=axes[1, 1])
    ms1.Yield_ACF_D2D.plot(ax=axes[1, 1])
    ms2.Yield_ACF_D2D.plot(ax=axes[1, 1])

    axes[1, 1].legend(['I0', 'I1', 'I2'])
    axes[1, 1].set_title('Yield, D2D')

    ms0.Yield_ACF_Vul.plot(ax=axes[2, 1])
    ms1.Yield_ACF_Vul.plot(ax=axes[2, 1])
    ms2.Yield_ACF_Vul.plot(ax=axes[2, 1])

    axes[2, 1].legend(['I0', 'I1', 'I2'])
    axes[2, 1].set_title('Yield, Vul')


    fig.tight_layout()

    plt.show()
