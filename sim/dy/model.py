import numpy as np
import pandas as pd
from sim.components import Demography, Transmission, Progression, Cascade
from sim.intv import Intervention
from sim.util import simulate, update_intv
import sim.dy.keys as I

__author__ = 'Chu-Chang Ku'
__all__ = ['Model']


class Model:
    def __init__(self, inp, year0=2010):
        self.Inputs = inp

        intv = Intervention()
        self.Demography = Demography(I, intv)
        self.Transmission = Transmission(I, intv)
        self.Progression = Progression(I, intv)
        self.Cascade = Cascade(I, intv)
        self.__intervention = intv

        self.Year0 = year0

    @property
    def Intervention(self):
        return self.__intervention

    @Intervention.setter
    def Intervention(self, intv):
        if isinstance(intv, dict):
            intv = Intervention.parse_obj(intv)

        self.__intervention = intv
        self.Demography.Intervention = intv
        self.Transmission.Intervention = intv
        self.Progression.Intervention = intv
        self.Cascade.Intervention = intv

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

        pars['trans_dr'] = trans = np.zeros((I.N_State_TB, I.N_State_Strata))
        trans[I.Infectious_Sn_DR] = pars['rr_inf_sn']
        trans[I.Infectious_Sp_DR] = 1

        pars['mixing'] = np.ones((I.N_State_Strata, I.N_State_Strata))

        return pars

    def get_y0(self, pars):
        y0 = np.zeros((I.N_State_TB, I.N_State_Strata))
        n0 = self.Inputs['N0'] * np.array([0.95, 0.05])

        y0[I.Sym_Sp_DS] = 1e-2 * n0
        y0[I.SLat_DS] = 0.4 * n0
        y0[I.U] = n0 - y0.sum(0)
        return y0

    def collect_calc(self, t, y, pars):
        #t = max(t, self.Year0)
        calc = dict()
        self.Demography(t, y, pars, calc)
        self.Transmission(t, y, pars, calc)
        self.Progression(t, y, pars, calc)
        self.Cascade(t, y, pars, calc)
        return calc

    def __call__(self, t, y, pars):
        y = y.reshape((I.N_State_TB, I.N_State_Strata))

        calc = self.collect_calc(t, y, pars)

        dy = np.zeros_like(y)
        # Infection
        dy -= calc['infection_ds'] + calc['infection_dr']
        dy[I.FLat_DS] += calc['infection_ds'].sum(0)
        dy[I.FLat_DR] += calc['infection_dr'].sum(0)

        # Incidence
        dy[I.FLat] -= calc['act']
        dy[I.SLat] -= calc['react']
        dy[I.RLow] -= calc['rel_tc']
        dy[I.RHigh] -= calc['rel_td']
        dy[I.RSt] -= calc['rel_st']

        dy[I.Asym] += calc['inc_smr']

        # Progression
        dy[I.FLat] -= calc['stab_fl']
        dy[I.SLat] += calc['stab_fl']
        dy[I.RLow] -= calc['stab_tc']
        dy[I.RHigh] -= calc['stab_td']
        dy[I.RSt] += calc['stab_tc'] + calc['stab_td']
        #
        dy[I.Asym] -= calc['sc_a'] + calc['sym_onset']
        dy[I.Sym] += calc['sym_onset'] - calc['sc_s']
        dy[I.ExSym] -= calc['sc_c']

        sc_asc = calc['sc_a'] + calc['sc_s'] + calc['sc_c']

        dy[I.RHigh_DS] += sc_asc[0] + sc_asc[1]
        dy[I.RHigh_DR] += sc_asc[2] + sc_asc[3]

        # Smear conversion
        con_a, con_s, con_c = calc['convert_a'], calc['convert_s'], calc['convert_c']
        dy[I.Asym_Sn] -= con_a
        dy[I.Asym_Sp] += con_a
        dy[I.Sym_Sn] -= con_s
        dy[I.Sym_Sp] += con_s
        dy[I.ExSym_Sn] -= con_c
        dy[I.ExSym_Sp] += con_c

        # DR development
        develop_dr_pub, develop_dr_pri = calc['develop_dr_pub'], calc['develop_dr_pri']
        dy[[I.Txf_Pub_Sn_DS, I.Txf_Pub_Sp_DS]] -= develop_dr_pub
        dy[[I.Txf_Pub_Sn_DR, I.Txf_Pub_Sp_DR]] += develop_dr_pub
        dy[[I.Txf_Pri_Sn_DS, I.Txf_Pri_Sp_DS]] -= develop_dr_pri
        dy[[I.Txf_Pri_Sn_DR, I.Txf_Pri_Sp_DR]] += develop_dr_pri

        # Dx
        det_1_pub_s, det_1_pub_c = calc['det_txf_pub_s'], calc['det_txf_pub_c']
        det_2_pub_s, det_2_pub_c = calc['det_txs_pub_s'], calc['det_txs_pub_c']
        det_1_pri_s, det_1_pri_c = calc['det_txf_pri_s'], calc['det_txf_pri_c']

        fn_s = calc['fn_pub_s'] + calc['fn_pri_s']

        dy[I.Sym] -= det_1_pub_s + det_2_pub_s + det_1_pri_s + fn_s
        dy[I.ExSym] += fn_s - (det_1_pub_c + det_2_pub_c + det_1_pri_c)
        dy[I.Txf_Pub] += det_1_pub_s + det_1_pub_c
        dy[I.Txf_Pri] += det_1_pri_s + det_1_pri_c
        dy[I.Txs_Pub] += det_2_pub_s + det_2_pub_c
        dy[I.Txs_Pri] += 0

        acf_1_pub_s, acf_1_pub_c = calc['acf_txf_pub_s'], calc['acf_txf_pub_c']
        acf_2_pub_s, acf_2_pub_c = calc['acf_txs_pub_s'], calc['acf_txs_pub_c']
        dy[I.Sym] -= acf_1_pub_s + acf_2_pub_s
        dy[I.ExSym] -= acf_1_pub_c + acf_2_pub_c
        dy[I.Txf_Pub] += acf_1_pub_s + acf_1_pub_c
        dy[I.Txs_Pub] += acf_2_pub_s + acf_2_pub_c

        # Tx
        tc_1_pub, td_1_pub = calc['tx_succ_txf_pub'], calc['tx_ltfu_txf_pub']
        tc_1_pri, td_1_pri = calc['tx_succ_txf_pri'], calc['tx_ltfu_txf_pri']
        tc_2_pub, td_2_pub = calc['tx_succ_txs_pub'], calc['tx_ltfu_txs_pub']
        tc_2_pri, td_2_pri = calc['tx_succ_txs_pri'], calc['tx_ltfu_txs_pri']

        dy[I.Txf_Pub] -= tc_1_pub + td_1_pub
        dy[I.Txf_Pri] -= tc_1_pri + td_1_pri
        dy[I.Txs_Pub] -= tc_2_pub + td_2_pub
        dy[I.Txs_Pri] -= tc_2_pri + td_2_pri

        tc = tc_1_pub + tc_1_pri + tc_2_pub + tc_2_pri
        td = td_1_pub + td_1_pri + td_2_pub + td_2_pri
        dy[I.RLow] += tc[[0, 2]] + tc[[1, 3]]
        dy[I.RHigh] += td[[0, 2]] + td[[1, 3]]

        tx_switch_pub = calc['tx_switch_pub']
        dy[I.Txf_Pub] -= tx_switch_pub
        dy[I.Txs_Pub] += tx_switch_pub

        # Self-clearance
        dy[I.SLat] -= calc['clear_sl']
        dy[I.RSt] -= calc['clear_rst']
        dy[I.U] += (calc['clear_sl'] + calc['clear_rst']).sum(0)

        # Demography
        dy[I.U, 0] += calc['births']
        dy -= calc['deaths'] + calc['deaths_tb']

        dy[:, 0] -= calc['prog_comorb']
        dy[:, 1] += calc['prog_comorb']

        return dy.reshape(-1)

    def measure(self, t, y, pars):
        y = y.reshape((I.N_State_TB, I.N_State_Strata))

        calc = self.collect_calc(t, y, pars)

        mea = {'Time': t}
        self.Demography.measure(t, y, pars, calc, mea)
        self.Transmission.measure(t, y, pars, calc, mea)
        self.Progression.measure(t, y, pars, calc, mea)
        self.Cascade.measure(t, y, pars, calc, mea)

        return mea

    @staticmethod
    def dfe(t, y, pars):
        ntb = y.reshape((I.N_State_TB, I.N_State_Strata))[I.Infectious].sum()
        return ntb - 0.5

    dfe.terminal = True
    dfe.direction = -1

    def simulate(self, p):
        p0 = p
        p = self.update_parameters(p)

        ys, ms, msg = simulate(self,
                               pars=p,
                               t_warmup=300,
                               t_out=np.linspace(1970, 2020, int(50 * 2) + 1),
                               dfe=self.dfe)
        msg['pars'] = p0
        return ys, ms, msg

    def simulate_onward(self, y0, p, intv=None, t_end=2030, dt=0.5):
        if intv is None:
            intv = self.Intervention

        t_start = 2020
        p0 = p
        if 'sus' not in p:
            p = self.update_parameters(p)

        n_ts = int((t_end - t_start) / dt) + 1
        ys, ms, msg = update_intv(self, y0, pars=p,
                                  intv=intv,
                                  t_out=np.linspace(t_start, t_end, n_ts),
                                  dfe=self.dfe)
        msg['pars'] = p0
        return ys, ms, msg


if __name__ == '__main__':
    from sim.dy.prior import get_bn
    from sim import load_inputs
    from sims_pars import sample
    import matplotlib.pylab as plt

    inputs = load_inputs('../../data/pars.json')

    m = Model(inputs, year0=1970)

    sc = get_bn()
    p0 = sample(sc, {'rr_risk_comorb': 20, 'beta_dr': 0, 'r_mdr_tx': 0})

    ys, ms, msg = m.simulate(p0)
    ys = ys.y.T[-1]
    _, ms1, _ = m.simulate_onward(ys, p0)
    _, ms2, _ = m.simulate_onward(ys, p0, intv={'ACFPlain': {'R_ACF': 0.2, 'Type': 'high'}})
    #_, ms2, _ = m.simulate_onward(ys, p0, intv={'ACF': {'Scale': 0.2, 'Type': 'mod'}})

    ms = pd.concat([ms, ms1.iloc[1:]])

    fig, axes = plt.subplots(2, 2)

    # ms = ms[ms.index > 2000]
    # ms.Pop.plot(ax=axes[0, 0])
    # ms.Pop_RiskLo.plot(ax=axes[0, 0])
    # ms.Pop_RiskHi.plot(ax=axes[0, 0])

    ms.PropComorb.plot(ax=axes[0, 0])

    ms.LTBI.plot(ax=axes[0, 1])
    ms.LTBI_RiskLo.plot(ax=axes[0, 1])
    ms.LTBI_RiskHi.plot(ax=axes[0, 1])

    # ms.RR_inf_comorb.plot()
    # ms.RR_inc_comorb.plot()

    # ms.Prev.plot()
    # ms.Prev_RiskLo.plot()
    # ms.Prev_RiskHi.plot()
    # ms2.Prev_RiskLo.plot()
    # ms2.Prev_RiskHi.plot()
    # ms.PrPrev_DR.plot()

    # ms.IncR.plot()
    # ms.Prev.plot()
    # ms.MorR.plot()
    # ms.CNR_RiskLo.plot(ax=axes[1, 1])
    # ms1.CNR_RiskLo.plot(ax=axes[1, 1])
    # ms2.CNR_RiskLo.plot(ax=axes[1, 1])

    ms.CNR.plot(ax=axes[1, 1])
    ms1.CNR.plot(ax=axes[1, 1])
    ms2.CNR.plot(ax=axes[1, 1])

    # ms.PrSp_Asym.plot()
    # ms.PrSp_Sym.plot()
    # ms.PrSym.plot()
    # ms.PrDR_Inc.plot()
    ms.IncR_RiskHi.plot(ax=axes[1, 0])
    ms1.IncR_RiskHi.plot(ax=axes[1, 0])
    ms2.IncR_RiskHi.plot(ax=axes[1, 0])

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
