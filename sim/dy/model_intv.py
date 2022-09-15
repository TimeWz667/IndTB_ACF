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
    def __init__(self, year0=1970):
        Model.__init__(self, year0)
        self.Parent = Model(year0)

    @staticmethod
    def _calc_bg_acf(y,  r_acf, eli, pos_screen, pos_confirm, p_dst, mea=False, label='bg'):
        arrived = r_acf * y
        screened = arrived * eli
        confirmed = screened * pos_screen
        pos = confirmed * pos_confirm
        tp_ds, tp_dr = pos[I.Infectious_DS], pos[I.Infectious_DR]
        tp_dr_fl = tp_dr * (1 - p_dst)
        tp_dr_sl = tp_dr - tp_dr_fl

        if mea:
            n = y.sum()
            return {
                f'ACF_{label}_Footfall': arrived.sum() / n,
                f'ACF_{label}_Screened': screened.sum() / n,
                f'ACF_{label}_Confirmed': confirmed.sum() / n,
                f'ACF_{label}_Yield': pos.sum() / n,
                f'ACF_{label}_TP': (tp_ds + tp_dr).sum() / n,
                f'ACF_{label}_DS_Fl': tp_ds.sum() / n,
                f'ACF_{label}_DR_Fl': tp_dr_fl.sum() / n,
                f'ACF_{label}_DR_Sl': tp_dr_sl.sum() / n,
            }
        else:
            dy = np.zeros_like(y)
            dy[I.Infectious_DS, :2] -= tp_ds
            dy[I.Infectious_DR, :2] -= tp_dr
            dy[I.Txf_Pub_DS, :2] += tp_ds.sum(0)
            dy[I.Txf_Pub_DR, :2] += tp_dr_fl.sum(0)
            dy[I.Txs_Pub_DR, :2] += tp_dr_sl.sum(0)
            return dy

    @staticmethod
    def _calc_vul_acf(y, cov, r_fu, r_lost, eli, pos_screen, pos_confirm, p_dst, mea=False, label='vul'):
        n = y.sum()
        n_target = cov * y.sum()
        r_acf = n_target / (eli * y[:, :2]).sum()
        arrived = r_acf * y[:, :2]
        screened = arrived * eli
        confirmed = screened * pos_screen
        pos = confirmed * pos_confirm
        tp_ds, tp_dr = pos[I.Infectious_DS], pos[I.Infectious_DR]
        tp_dr_fl = tp_dr * (1 - p_dst)
        tp_dr_sl = tp_dr - tp_dr_fl

        fp_tpt = pos[I.LTBI]

        if mea:
            return {
                f'ACF_{label}_Footfall': arrived.sum() / n,
                f'ACF_{label}_Screened': screened.sum() / n,
                f'ACF_{label}_Confirmed': confirmed.sum() / n,
                f'ACF_{label}_Yield': pos.sum() / n,
                f'ACF_{label}_TP': (tp_ds + tp_dr).sum() / n,
                f'ACF_{label}_DS_Fl': tp_ds.sum() / n,
                f'ACF_{label}_DR_Fl': tp_dr_fl.sum() / n,
                f'ACF_{label}_DR_Sl': tp_dr_sl.sum() / n,
            }
        else:
            dy = np.zeros_like(y)
            dy[I.Infectious_DS, :2] -= tp_ds
            dy[I.Infectious_DR, :2] -= tp_dr
            dy[I.Txf_Pub_DS, :2] += tp_ds.sum(0)
            dy[I.Txf_Pub_DR, :2] += tp_dr_fl.sum(0)
            dy[I.Txs_Pub_DR, :2] += tp_dr_sl.sum(0)

            if r_fu > 0:
                lost = r_lost * y[:, 2:]
                dy[:, 2:] -= lost
                dy[:, :2] += lost

            return dy

    def __call__(self, t, y, pars, intv=None):
        y = y.reshape((I.N_State_TB, I.N_State_Strata * 2))
        dy = np.zeros_like(y)

        dy[:, :2] = self.Parent(t, y[:, :2], pars).reshape((-1, 2))
        if y[:, 2:].sum() > 0:
            pars_acf = dict(pars)

            dy[:, 2:] = self.Parent(t, y[:, 2:], pars_acf).reshape((-1, 2))

        # ACF, background
        pos_sym, pos_cxr, pos_xpert, eligible = pars['pos_sym'], pars['pos_cxr'], pars['pos_xpert'], pars['eli']
        r_acf_mu, r_acf_d2d, p_dst = pars['r_acf_mu'], pars['r_acf_d2d'], pars['acf_dst_sens']

        eli_vul = pars['eli_vul']
        r_acf_vul = r_acf_plain = 0

        if intv is not None:
            r_acf_mu, r_acf_d2d, p_dst = intv.modify_acf_bg(t, r_acf_mu, r_acf_d2d, p_dst)
            # MDU
            if r_acf_mu > 0:
                dy[:, :2] += self._calc_bg_acf(y[:, :2], r_acf_mu, eligible, pos_cxr, pos_xpert, p_dst, mea=False)
            # D2D
            if r_acf_d2d > 0:
                dy[:, :2] += self._calc_bg_acf(y[:, :2], r_acf_d2d, eligible, pos_sym, pos_xpert, p_dst, mea=False)

            # Vulnerability-led ACF
            r_acf_vul, r_fu, r_lost = intv.modify_acf_vul(t, r_acf_vul, 0, 0)
            if r_acf_vul > 0:
                dy += self._calc_vul_acf(y, r_acf_vul, r_fu, r_lost, eli_vul, pos_cxr, pos_xpert, p_dst, mea=False)

            # Plain ACF
            r_acf_plain, r_fu, r_lost = intv.modify_acf_plain(t, r_acf_plain, 0, 0)
            if r_acf_plain > 0:
                dy += self._calc_vul_acf(y, r_acf_plain, r_fu, r_lost, eligible, pos_cxr, pos_xpert, p_dst, mea=False)

        return dy.reshape(-1)

    def measure(self, t, y, pars, intv=None):
        y = y.reshape((I.N_State_TB, I.N_State_Strata * 2))

        # ACF metrics
        pos_sym, pos_cxr, pos_xpert, eligible = pars['pos_sym'], pars['pos_cxr'], pars['pos_xpert'], pars['eli']
        r_acf_mu, r_acf_d2d, p_dst = pars['r_acf_mu'], pars['r_acf_d2d'], pars['acf_dst_sens']

        eli_vul = pars['eli_vul']
        r_acf_vul = r_acf_plain = 0

        mea_acf = {}

        if intv is not None:
            r_acf_mu, r_acf_d2d, p_dst = intv.modify_acf_bg(t, r_acf_mu, r_acf_d2d, p_dst)

            # MDU
            mea_acf.update(self._calc_bg_acf(y[:, :2], r_acf_mu, eligible, pos_cxr, pos_xpert, p_dst,
                                             mea=True, label='MDU'))

            # D2D
            mea_acf.update(self._calc_bg_acf(y[:, :2], r_acf_mu, eligible, pos_cxr, pos_xpert, p_dst,
                                             mea=True, label='D2D'))

            # Vulnerability-led ACF
            r_acf_vul, r_fu, r_lost = intv.modify_acf_vul(t, r_acf_vul, 0, 0)
            mea_acf.update(self._calc_vul_acf(y, r_acf_vul, r_fu, r_lost, eli_vul, pos_cxr, pos_xpert, p_dst,
                                              mea=True, label='Vul'))

            # Plain ACF
            r_acf_plain, r_fu, r_lost = intv.modify_acf_plain(t, r_acf_plain, 0, 0)
            mea_acf.update(self._calc_vul_acf(y, r_acf_plain, r_fu, r_lost, eligible, pos_cxr, pos_xpert, p_dst,
                                              mea=True, label='Plain'))

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
        pos_sym, pos_cxr, pos_xpert, eligible = p['pos_sym'], p['pos_cxr'], p['pos_xpert'], p['eli']

        p['r_acf_mu'] = yield_mu / ((y0 * eligible * pos_cxr * pos_xpert).sum() / y0.sum())
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

            opt = minimize_scalar(fn, 0.99, args=(y0, ), method='bounded', bounds=(0.5, 1))
            spec_sym = opt.x

        pos_sym = np.zeros_like(y0)
        pos_sym[I.LTBI + I.Asym] = (1 - spec_sym)
        pos_sym[I.U] = (1 - spec_sym)
        pos_sym[I.Sym + I.ExSym] = 1

        pos_vul = np.zeros_like(y0)
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

    _, ms0, _ = m0.simulate_onward(y1, p1, intv={'D2D': {'Scale': 0}, 'MU': {'Scale': 0}})
    # _, ms1, _ = m0.simulate_onward(y1, p1, intv={'D2D': {'Scale': 1}, 'MU': {'Scale': 1}})
    # _, ms2, _ = m0.simulate_onward(y1, p1, intv={'D2D': {'Scale': 2}, 'MU': {'Scale': 2}})
    _, ms1, _ = m0.simulate_onward(y1, p1, intv={'VulACF': {'Coverage': 0.15}})
    _, ms2, _ = m0.simulate_onward(y1, p1, intv={'PlainACF': {'Coverage': 0.15}})
    # _, ms2, _ = m.simulate_onward(y1, p0, intv={'ACF': {'Yield': .05, 'HiRisk': False}})
    # _, ms3, _ = m.simulate_onward(y1, p0, intv={'ACF': {'Yield': .03, 'HiRisk': True}})

    print('MU', ms1.ACF_MDU_Yield[2022.5] * 1e5, 104.5 / 3e6 * 1e5)
    print('D2D', ms1.ACF_D2D_Yield[2022.5] * 1e5, 430 / 3e6 * 1e5)

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
    axes[0, 1].set_title('Yield, MU')

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
