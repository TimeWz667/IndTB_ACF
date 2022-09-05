import numpy as np
import pandas as pd
from sim.components import Demography, Transmission, Progression, Cascade
from sim.intv import Intervention
from sim.util import simulate
import sim.dy.keys as I
from scipy.integrate import solve_ivp

__author__ = 'Chu-Chang Ku'
__all__ = ['Model']


class Model:
    def __init__(self, inp, year0=2010):
        self.Inputs = inp

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
        n0 = np.array([1 - p_hi, p_hi]) * self.Inputs['N0']

        y0[I.Sym_DS] = 1e-2 * n0
        y0[I.SLat_DS] = 0.4 * n0
        y0[I.U] = n0 - y0.sum(0)
        return y0

    def __call__(self, t, y, pars, intv=None):
        y = y.reshape((I.N_State_TB, I.N_State_Strata))

        dy = self.Demography.calc_dy(t, y, pars, intv)
        dy += self.Transmission.calc_dy(t, y, pars, intv)
        dy += self.Progression.calc_dy(t, y, pars, intv)
        dy += self.Cascade.calc_dy(t, y, pars, intv)

        return dy.reshape(-1)

    def measure(self, t, y, pars, intv=None):
        y = y.reshape((I.N_State_TB, I.N_State_Strata))

        mea = {'Time': t}
        self.Demography.measure(t, y, pars, intv, mea)
        self.Transmission.measure(t, y, pars, intv, mea)
        self.Progression.measure(t, y, pars, intv, mea)
        self.Cascade.measure(t, y, pars, intv, mea)

        return mea

    @staticmethod
    def dfe(t, y, pars, intv=None):
        ntb = y.reshape((I.N_State_TB, I.N_State_Strata))[I.Infectious].sum()
        return ntb - 0.5

    dfe.terminal = True
    dfe.direction = -1

    def simulate_to_fit(self, p):
        if 'sus' not in p:
            p = self.update_parameters(p)

        t_start = 2020
        y0 = self.get_y0(p).reshape(-1)
        ys0 = solve_ivp(self, [t_start - 300, t_start], y0, args=(p,), events=self.dfe, method='RK23')
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

    def form_non_tb_population(self, t0, y0, p, ppv0, spec):
        mea = dict()
        self.Cascade.measure(t0, y0.reshape((-1, 2)), p, None, mea)
        tp = mea['TP']
        fp = tp * (1 - ppv0) / ppv0
        p['PPV'] = ppv0
        p['spec0'] = spec
        p['NonTB'] = fp / (p['r_cs_s'] * (1 - spec))

    def simulate_onward(self, y0, p, intv=None, t_end=2030, dt=0.5, ppv0=0.85, spec=0.99):
        t_start = 2022
        t_out = np.linspace(t_start, t_end, int((t_end - t_start) / dt) + 1)

        if 'sus' not in p:
            p = self.update_parameters(p)
        if 'PPV' not in p or (p['PPV'] is not ppv0) or 'NonTB' not in p:
            self.form_non_tb_population(t_start, y0, p, ppv0, spec)

        if intv is not None and not isinstance(intv, Intervention):
            intv = Intervention.parse_obj(intv)

        ys = solve_ivp(self, [t_start, t_end], y0, args=(p, intv), events=self.dfe, dense_output=True)

        if len(ys.t_events[0]) > 0 or not ys.success:
            return None, None, {'succ': False, 'res': 'DFE reached'}

        ms = [self.measure(t, ys.sol(t), p, intv) for t in t_out if t >= t_start]
        ms = pd.DataFrame(ms).set_index('Time')

        return ys, ms, {'succ': True}


if __name__ == '__main__':
    import json
    from sim import load_inputs
    import matplotlib.pylab as plt

    with open('../../data/test_priors.json', 'r') as f:
        prior = json.load(f)

    inputs = load_inputs('../../data/pars.json')

    m = Model(inputs, year0=1970)

    p0 = prior[2]
    p0.update({'beta_ds': 15, 'rr_risk_comorb': 20, 'rr_beta_dr': 1.05, 'p_comorb': 0.3})

    _, ms, _ = m.simulate_to_fit(p0)
    print(ms)

    y1, ms, msg = m.simulate_to_baseline(p0, 2022)

    _, ms1, _ = m.simulate_onward(y1, p0)
    _, ms2, _ = m.simulate_onward(y1, p0, intv={'ACF': {'Yield': .05, 'HiRisk': False}})
    _, ms3, _ = m.simulate_onward(y1, p0, intv={'ACF': {'Yield': .03, 'HiRisk': True}})

    fig, axes = plt.subplots(2, 3)

    print(ms1.IncR.tail(5))

    ms1.Pop.plot(ax=axes[0, 0])
    ms1.Pop_RiskLo.plot(ax=axes[0, 0])
    ms1.Pop_RiskHi.plot(ax=axes[0, 0])

    # ms.PropComorb.plot()

    ms1.LTBI.plot(ax=axes[0, 1])
    ms1.LTBI_RiskLo.plot(ax=axes[0, 1])
    ms1.LTBI_RiskHi.plot(ax=axes[0, 1])

    ms1.IncR.plot(ax=axes[0, 2])
    ms2.IncR.plot(ax=axes[0, 2])
    ms3.IncR.plot(ax=axes[0, 2])

    ms1.R_ACF_Reached.plot(ax=axes[1, 2])
    ms2.R_ACF_Reached.plot(ax=axes[1, 2])
    ms3.R_ACF_Reached.plot(ax=axes[1, 2])

    # ms.RR_inf_comorb.plot()
    # ms.RR_inc_comorb.plot()

    # ms.Prev.plot()
    # ms.Prev_RiskLo.plot()
    # ms.Prev_RiskHi.plot()
    # ms2.Prev_RiskLo.plot()
    # ms2.Prev_RiskHi.plot()
    # ms.PrDR_Inc.plot(ax=axes[1, 0])
    ms1.PPV_Acf.plot(ax=axes[1, 0])
    ms2.PPV_Acf.plot(ax=axes[1, 0])
    ms3.PPV_Acf.plot(ax=axes[1, 0])

    # ms.IncR.plot()
    # ms.Prev.plot()
    # ms.MorR.plot()
    # ms.CNR_RiskLo.plot(ax=axes[1, 1])
    # ms1.CNR_RiskLo.plot(ax=axes[1, 1])
    # ms2.CNR_RiskLo.plot(ax=axes[1, 1])

    ms1.TP.plot(ax=axes[1, 1])
    ms2.TP.plot(ax=axes[1, 1])
    ms3.TP.plot(ax=axes[1, 1])

    # ms.PrSp_Asym.plot()
    # ms.PrSp_Sym.plot()
    # ms.PrSym.plot()
    # ms.PrDR_Inc.plot()
    # ms.IncR_RiskHi.plot(ax=axes[1, 0])
    # ms1.IncR_RiskHi.plot(ax=axes[1, 0])
    # ms2.IncR_RiskHi.plot(ax=axes[1, 0])

    # ms.IncR_Remote.plot(ax=axes[1, 0])
    # ms1.IncR_Remote.plot(ax=axes[1, 0])
    # ms2.IncR_Remote.plot(ax=axes[1, 0])

    plt.legend()
    plt.show()

    # print(ms3)
