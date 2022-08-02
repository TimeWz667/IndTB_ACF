import numpy as np
import pandas as pd
from sim.components import Demography, Transmission, Progression, Cascade
from sim.intv import Intervention
from sim.util import simulate, update
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

    def update_parameters(self, pars):
        pars = dict(pars)

        pars.update({
            'r_succ_txf_pub': pars['r_succ_txf'],
            'r_succ_txf_pri': pars['r_succ_txf'],
            'r_ltfu_txf_pub': pars['r_succ_txf'] * (1 - pars['p_succ_txf_pub']) / pars['p_succ_txf_pub'],
            'r_ltfu_txf_pri': pars['r_succ_txf'] * (1 - pars['p_succ_txf_pri']) / pars['p_succ_txf_pri'],
        })

        pars.update({
            'r_succ_txs_pub': pars['r_succ_txs'],
            'r_succ_txs_pri': 0,
            'r_ltfu_txs_pub': pars['r_succ_txs'] * (1 - pars['p_succ_txs_pub']) / pars['p_succ_txs_pub'],
            'r_ltfu_txs_pri': pars['r_succ_txs']
        })

        pars['sus'] = sus = np.zeros((I.N_State_TB, I.N_State_Strata))
        sus[I.U] = 1
        sus[I.SLat] = pars['rr_sus_ltbi']
        sus[I.RLow] = pars['rr_sus_ltbi']
        sus[I.RHigh] = pars['rr_sus_ltbi']
        sus[I.RSt] = pars['rr_sus_ltbi']
        # sus[: I.RiskHi] *= pars['rr_risk_comorb']

        pars['trans_ds'] = trans = np.zeros((I.N_State_TB, I.N_State_Strata))
        trans[I.Infectious_Sn_DS] = pars['rr_inf_sn']
        trans[I.Infectious_Sp_DS] = 1
        trans[I.Asym] *= pars['rr_inf_asym']

        pars['trans_dr'] = trans = np.zeros((I.N_State_TB, I.N_State_Strata))
        trans[I.Infectious_Sn_DR] = pars['rr_inf_sn']
        trans[I.Infectious_Sp_DR] = 1
        trans[I.Asym] *= pars['rr_inf_asym']

        pars['beta_dr'] = pars['beta_ds'] * pars['rr_beta_dr']

        return pars

    def get_y0(self, pars):
        y0 = np.zeros((I.N_State_TB, I.N_State_Strata))
        p_hi = pars['p_comorb']
        n0 = np.array([1 - p_hi, p_hi]) * self.Inputs['N0']

        y0[I.Sym_Sp_DS] = 1e-2 * n0
        y0[I.SLat_DS] = 0.4 * n0
        y0[I.U] = n0 - y0.sum(0)
        return y0

    def collect_calc(self, t, y, pars, intv):
        calc = dict()
        self.Demography(t, y, pars, intv, calc)
        self.Transmission(t, y, pars, intv, calc)
        self.Progression(t, y, pars, intv, calc)
        self.Cascade(t, y, pars, intv, calc)
        return calc

    def __call__(self, t, y, pars, intv=None):
        y = y.reshape((I.N_State_TB, I.N_State_Strata))

        calc = dict()
        self.Demography(t, y, pars, intv, calc)
        self.Transmission(t, y, pars, intv, calc)

        dy = np.zeros_like(y)

        # Infection
        dy -= calc['infection_ds'] + calc['infection_dr']
        dy[I.FLat_DS] += calc['infection_ds'].sum(0)
        dy[I.FLat_DR] += calc['infection_dr'].sum(0)

        # Demography
        dy[I.U, 0] += calc['births']
        dy -= calc['deaths'] + calc['deaths_tb']

        dy[:, 0] -= calc['prog_comorb']
        dy[:, 1] += calc['prog_comorb']

        dy += self.Progression.calc_dy2(t, y, pars, intv)
        dy += self.Cascade.calc_dy2(t, y, pars, intv)

        return dy.reshape(-1)

    def measure(self, t, y, pars, intv=None):
        y = y.reshape((I.N_State_TB, I.N_State_Strata))

        calc = self.collect_calc(t, y, pars, intv)

        mea = {'Time': t}
        self.Demography.measure(t, y, pars, intv, calc, mea)
        self.Transmission.measure(t, y, pars, intv, calc, mea)
        self.Progression.measure(t, y, pars, intv, calc, mea)
        self.Cascade.measure(t, y, pars, intv, calc, mea)

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

        t_start = 2006
        y0 = self.get_y0(p).reshape(-1)
        ys0 = solve_ivp(self, [t_start - 300, t_start], y0, args=(p,), events=self.dfe, method='RK23')
        if len(ys0.t_events[0]) > 0 or not ys0.success:
            return None, None, {'succ': False, 'res': 'DFE reached'}

        y1 = ys0.y[:, -1]
        ys = solve_ivp(self, [t_start, 2018], y1, args=(p,), events=self.dfe, method='RK23', dense_output=True)
        if len(ys.t_events[0]) > 0 or not ys.success:
            return None, None, {'succ': False, 'res': 'DFE reached'}

        mea = [self.measure(t, ys.sol(t), p) for t in [2006, 2012, 2018]]
        return ys, mea, {'succ': True}

    def simulate(self, p):
        if 'sus' not in p:
            p = self.update_parameters(p)

        ys, ms, msg = simulate(self,
                               pars=p,
                               t_warmup=300,
                               t_out=np.linspace(1970, 2020, int(50 * 2) + 1),
                               dfe=self.dfe)
        return ys, ms, msg

    def form_non_tb_population(self, t0, y0, pars, ppv0, spec):
        calc = dict()
        self.Cascade(t0, y0.reshape((I.N_State_TB, 2)), pars, None, calc)
        tp = calc['det_txf_pub_s'] + calc['det_txf_pub_c'] + calc['det_txs_pub_s'] + calc['det_txs_pub_c']
        tp += calc['det_txf_pri_s'] + calc['det_txf_pri_c']
        tp = tp.sum()
        fp = tp * (1 - ppv0) / ppv0

        pars['NonTB'] = fp / (pars['r_cs_s'] * (1 - 0.99)) / y0.sum()
        pars['spec0'] = spec

    def simulate_onward(self, y0, p, intv=None, t_end=2030, dt=0.5, ppv0=0.85, spec=0.99):
        t_start = 2020
        p0 = p
        if 'sus' not in p:
            p = self.update_parameters(p)
        if 'NonTB' not in p:
            self.form_non_tb_population(t_start, y0, p, ppv0, spec)

        if intv is not None and not isinstance(intv, Intervention):
            intv = Intervention.parse_obj(intv)

        n_ts = int((t_end - t_start) / dt) + 1
        ys, ms, msg = update(self, y0, pars=p,
                                  intv=intv,
                                  t_out=np.linspace(t_start, t_end, n_ts),
                                  dfe=self.dfe)
        msg['pars'] = p0
        return ys, ms, msg


if __name__ == '__main__':
    import json
    from sim import load_inputs
    import matplotlib.pylab as plt

    with open('../../data/test_priors.json', 'r') as f:
        prior = json.load(f)

    inputs = load_inputs('../../data/pars.json')

    m = Model(inputs, year0=1970)

    p0 = prior[0]
    p0.update({'beta_ds': 11, 'rr_risk_comorb': 20, 'rr_beta_dr': 1.02, 'p_comorb': 0.3})

    ys, ms, msg = m.simulate(p0)
    ys = ys.y.T[-1]
    _, ms1, _ = m.simulate_onward(ys, p0)
    _, ms2, _ = m.simulate_onward(ys, p0, intv={'ACF': {'Yield': 1, 'Type': 'high', 'HiRisk': False}})
    _, ms3, _ = m.simulate_onward(ys, p0, intv={'ACF': {'Yield': 1, 'Type': 'high', 'HiRisk': True}})

    ms = pd.concat([ms, ms1.iloc[1:]])

    fig, axes = plt.subplots(2, 3)

    ms = ms[ms.index > 2000]

    print(ms.IncR.tail(5))

    ms.Pop.plot(ax=axes[0, 0])
    ms.Pop_RiskLo.plot(ax=axes[0, 0])
    ms.Pop_RiskHi.plot(ax=axes[0, 0])

    # ms.PropComorb.plot()

    ms.LTBI.plot(ax=axes[0, 1])
    ms.LTBI_RiskLo.plot(ax=axes[0, 1])
    ms.LTBI_RiskHi.plot(ax=axes[0, 1])

    ms.IncR.plot(ax=axes[0, 2])
    ms2.IncR.plot(ax=axes[0, 2])
    ms3.IncR.plot(ax=axes[0, 2])

    ms.R_ACF_Reached.plot(ax=axes[1, 2])
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
    ms.PPV_ACF.plot(ax=axes[1, 0])
    ms2.PPV_ACF.plot(ax=axes[1, 0])
    ms3.PPV_ACF.plot(ax=axes[1, 0])

    # ms.IncR.plot()
    # ms.Prev.plot()
    # ms.MorR.plot()
    # ms.CNR_RiskLo.plot(ax=axes[1, 1])
    # ms1.CNR_RiskLo.plot(ax=axes[1, 1])
    # ms2.CNR_RiskLo.plot(ax=axes[1, 1])

    ms.TP.plot(ax=axes[1, 1])
    ms1.TP.plot(ax=axes[1, 1])
    ms2.TP.plot(ax=axes[1, 1])

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

    # ms.IncR_DS.plot()
    # ms.IncR_DR.plot()
    # ms.PrDR_Inc.plot()

    # ms.IncR_DS.plot(ax=axes[1, 0])
    # ms.IncR_DR.plot(ax=axes[1, 0])
    # ms1.IncR_DR.plot()
    # ms2.IncR_DR.plot()

    plt.legend()
    plt.show()

    # print(ms3)
