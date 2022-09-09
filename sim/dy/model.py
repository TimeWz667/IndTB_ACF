import numpy as np
import pandas as pd
from sim.components import Demography, Transmission, Progression, Cascade
from sim.util import simulate
import sim.dy.keys as I
from scipy.integrate import solve_ivp

__author__ = 'Chu-Chang Ku'
__all__ = ['Model']


class Model:
    def __init__(self, year0=1970):
        self.Demography = Demography(I)
        self.Transmission = Transmission(I)
        self.Progression = Progression(I)
        self.Cascade = Cascade(I)

        self.Year0 = year0

    @staticmethod
    def update_parameters(pars):
        pars = dict(pars)

        pars['sus'] = sus = np.zeros((I.N_State_TB, I.N_State_Strata))
        sus[I.U] = 1
        sus[I.SLat] = pars['rr_sus_ltbi']
        sus[I.RLow] = pars['rr_sus_ltbi']
        sus[I.RHigh] = pars['rr_sus_ltbi']
        sus[I.RSt] = pars['rr_sus_ltbi']

        pars['trans_ds'] = trans = np.zeros((I.N_State_TB, I.N_State_Strata))
        trans[I.Infectious_DS] = 1
        trans[I.Asym] *= pars['rr_inf_asym']

        pars['trans_dr'] = trans = np.zeros((I.N_State_TB, I.N_State_Strata))
        trans[I.Infectious_DR] = 1
        trans[I.Asym] *= pars['rr_inf_asym']

        pars['beta_dr'] = pars['beta_ds'] * pars['rr_beta_dr']

        for sector in ['pub', 'pri']:
            pars[f'r_succ_fl_{sector}_ds'] = pars['r_succ_txf']
            p_ltfu = 1 - pars[f'p_succ_txf_{sector}'] - pars['p_die_txf_ds']
            pars[f'r_ltfu_fl_{sector}_ds'] = pars['r_succ_txf'] * p_ltfu / pars[f'p_succ_txf_{sector}']

            pars[f'r_succ_fl_{sector}_dr'] = 0
            pars[f'r_ltfu_fl_{sector}_dr'] = pars['r_succ_txf']

        pars['r_succ_sl_pub_dr'] = pars['r_succ_txs']
        p_ltfu = 1 - pars['p_succ_txs_pub'] - pars['p_die_txs_dr']
        pars['r_ltfu_sl_pub_dr'] = pars['r_succ_txs'] * p_ltfu / pars['p_succ_txs_pub']

        pars['r_die_ds'] = pars['r_succ_txf'] * pars['p_die_txf_ds'] / pars['p_succ_txf_pub']
        pars['r_die_dr'] = pars['r_succ_txs'] * pars['p_die_txs_dr'] / pars['p_succ_txs_pub']

        return pars

    def get_y0(self, pars):
        y0 = np.zeros((I.N_State_TB, I.N_State_Strata))
        p_hi = pars['p_comorb']
        n0 = np.array([1 - p_hi, p_hi])

        y0[I.Sym_DS] = 1e-2 * n0
        y0[I.SLat_DS] = 0.4 * n0
        y0[I.U] = n0 - y0.sum(0)
        return y0

    def __call__(self, t, y, pars):
        y = y.reshape((I.N_State_TB, I.N_State_Strata))

        dy = self.Demography.calc_dy(t, y, pars)
        dy += self.Transmission.calc_dy(t, y, pars)
        dy += self.Progression.calc_dy(t, y, pars)
        dy += self.Cascade.calc_dy(t, y, pars)

        return dy.reshape(-1)

    def measure(self, t, y, pars):
        y = y.reshape((I.N_State_TB, I.N_State_Strata))

        mea = {'Time': t}
        self.Demography.measure(t, y, pars, mea)
        self.Transmission.measure(t, y, pars, mea)
        self.Progression.measure(t, y, pars, mea)
        self.Cascade.measure(t, y, pars, mea)

        return mea

    @staticmethod
    def dfe(t, y, pars):
        ntb = y.reshape((I.N_State_TB, I.N_State_Strata))[I.Infectious].sum()
        return ntb - 0.5

    dfe.terminal = True
    dfe.direction = -1

    def simulate_to_fit(self, p):
        if 'sus' not in p:
            p = self.update_parameters(p)

        t_start = 2020
        y0 = self.get_y0(p).reshape(-1)
        ys0 = solve_ivp(self, [t_start - 300, t_start], y0, args=(p, ), events=self.dfe, method='RK23')
        if len(ys0.t_events[0]) > 0 or not ys0.success:
            return None, None, {'succ': False, 'res': 'DFE reached'}

        return ys0, self.measure(t_start, ys0.y[:, -1], p), {'succ': True}

    def simulate_to_baseline(self, p, t_start=2022):
        if 'sus' not in p:
            p = self.update_parameters(p)

        y0 = self.get_y0(p).reshape(-1)
        ys = solve_ivp(self, [t_start - 300, t_start], y0, args=(p, ), events=self.dfe, method='RK23')
        if len(ys.t_events[0]) > 0 or not ys.success:
            return None, None, {'succ': False, 'res': 'DFE reached'}
        y1 = ys.y[:, -1]
        mea = self.measure(t_start, y1, p)
        return y1, mea, {'succ': True}

    def simulate(self, p):
        if 'sus' not in p:
            p = self.update_parameters(p)

        ys, ms, msg = simulate(self,
                               pars=p,
                               t_warmup=300,
                               t_out=np.linspace(1970, 2020, int(50 * 2) + 1),
                               dfe=self.dfe)
        return ys, ms, msg

    def simulate_onward(self, y0, p, t_end=2030, dt=0.5):
        t_start = 2022
        t_out = np.linspace(t_start, t_end, int((t_end - t_start) / dt) + 1)

        if 'sus' not in p:
            p = self.update_parameters(p)

        ys = solve_ivp(self, [t_start, t_end], y0, args=(p, ), events=self.dfe, dense_output=True)

        if len(ys.t_events[0]) > 0 or not ys.success:
            return None, None, {'succ': False, 'res': 'DFE reached'}

        ms = [self.measure(t, ys.sol(t), p) for t in t_out if t >= t_start]
        ms = pd.DataFrame(ms).set_index('Time')

        return ys, ms, {'succ': True}


if __name__ == '__main__':
    import json
    import matplotlib.pylab as plt

    with open('../../data/test_priors.json', 'r') as f:
        prior = json.load(f)

    m = Model(year0=1970)

    p0 = prior[2]
    p0.update({'beta_ds': 5, 'rr_risk_comorb': 20, 'rr_beta_dr': 1.05, 'p_comorb': 0.3})

    _, ms, _ = m.simulate_to_fit(p0)
    print(ms)

    y1, ms, msg = m.simulate_to_baseline(p0, 2022)

    _, ms1, _ = m.simulate_onward(y1, p0)
    # _, ms2, _ = m.simulate_onward(y1, p0, intv={'ACF': {'Yield': .05, 'HiRisk': False}})
    # _, ms3, _ = m.simulate_onward(y1, p0, intv={'ACF': {'Yield': .03, 'HiRisk': True}})

    fig, axes = plt.subplots(2, 2)

    print(ms1.IncR.tail(5))

    ms1.Pop.plot(ax=axes[0, 0])
    ms1.Pop_RiskLo.plot(ax=axes[0, 0])
    ms1.Pop_RiskHi.plot(ax=axes[0, 0])
    axes[0, 0].legend(['All', 'Lo', 'Hi'])
    axes[0, 0].set_title('Population')

    ms1.LTBI.plot(ax=axes[0, 1])
    ms1.LTBI_RiskLo.plot(ax=axes[0, 1])
    ms1.LTBI_RiskHi.plot(ax=axes[0, 1])
    axes[0, 1].legend(['All', 'Lo', 'Hi'])
    axes[0, 1].set_title('LTBI')

    ms1.IncR.plot(ax=axes[1, 0])
    ms1.IncR_RiskLo.plot(ax=axes[1, 0])
    ms1.IncR_RiskHi.plot(ax=axes[1, 0])
    axes[1, 0].legend(['All', 'Lo', 'Hi'])
    axes[1, 0].set_title('Incidence')

    ms1.Prev.plot(ax=axes[1, 1])
    ms1.Prev_RiskLo.plot(ax=axes[1, 1])
    ms1.Prev_RiskHi.plot(ax=axes[1, 1])
    axes[1, 1].legend(['All', 'Lo', 'Hi'])
    axes[1, 1].set_title('Prevalence')

    fig.tight_layout()

    plt.show()
