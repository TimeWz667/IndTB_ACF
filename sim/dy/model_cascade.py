import numpy as np
from sim.dy.model import Model
import sim.dy.keys as I
from sim.util import update_intv

__author__ = 'Chu-Chang Ku'
__all__ = ['ModelCascade']


class ModelCascade(Model):
    def __init__(self, inp, year0=2010, year=2022):
        Model.__init__(self, inp, year0=year0)
        self.Year = year

    def __call__(self, t, y, pars):
        y = y.reshape((I.N_State_TB, I.N_State_Strata))

        calc = self.collect_calc(self.Year, y, pars)

        dy = np.zeros_like(y)
        #
        # dy -= calc['infection_ds'] + calc['infection_dr']
        #
        # dy[I.SLat] += calc['lat']
        #
        # dy[I.SLat] -= calc['react']
        # dy[I.RLow] -= calc['rel_tc']
        # dy[I.RHigh] -= calc['rel_td']
        # dy[I.RSt] -= calc['rel_st']
        # #
        # dy[I.Asym] += calc['inc_smr']

        # Progression
        # dy[I.RLow] -= calc['stab_tc']
        # dy[I.RHigh] -= calc['stab_td']
        # dy[I.RSt] += calc['stab_tc'] + calc['stab_td']
        #
        dy[I.Asym] -= calc['sc_a'] + calc['sym_onset']
        dy[I.Sym] += calc['sym_onset'] - calc['sc_s']
        dy[I.ExSym] -= calc['sc_c']

        sc_asc = calc['sc_a'] + calc['sc_s'] + calc['sc_c']
        sc = np.zeros((2, I.N_State_Strata))
        sc[0] = sc_asc[0] + sc_asc[2]
        sc[1] = sc_asc[1] + sc_asc[3]

        # dy[I.RHigh] += sc[0] + sc[1]

        # Smear convertion
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

        # tc = tc_1_pub + tc_1_pri + tc_2_pub + tc_2_pri
        # td = td_1_pub + td_1_pri + td_2_pub + td_2_pri
        # dy[I.RLow] += tc[[0, 2]] + tc[[1, 3]]
        # dy[I.RHigh] += td[[0, 2]] + td[[1, 3]]

        tx_switch_pub = calc['tx_switch_pub']
        dy[I.Txf_Pub] -= tx_switch_pub
        dy[I.Txs_Pub] += tx_switch_pub


        # # Self-clearance
        # dy[I.SLat] -= calc['clear_sl']
        # dy[I.RSt] -= calc['clear_rst']
        # dy[I.U] += (calc['clear_sl'] + calc['clear_rst'])

        # Demography
        # dy[I.U, 0] += calc['births']
        dy -= calc['deaths'] + calc['deaths_tb']

        # dy[:, 0] -= calc['prog_comorb']
        # dy[:, 1] += calc['prog_comorb']

        # if t <= self.Year0:
        #     ns = y.sum(0, keepdims=True)
        #     ns[ns == 0] = 1e-10
        #     dy -= y / ns * dy.sum(0, keepdims=True)
        # else:
        #     pass
        # if t <= self.Year0:
        #     dy -= y / y.sum() * dy.sum()

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
        return 1

    dfe.terminal = True
    dfe.direction = -1

    def simulate(self, p):
        raise AttributeError

    def simulate_onward(self, y0, p, intv=None, t_end=2030, dt=0.5):
        if intv is None:
            intv = self.Intervention

        t_start = self.Year
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

    m = ModelCascade(inputs, year0=1970, year=2022)

    sc = get_bn()
    p0 = sample(sc, {'rr_risk_comorb': 20})

    y0 = np.zeros_like(m.get_y0(p0))
    y0[I.Asym_Sn_DS, I.RiskHi] = 1
    y0 = y0.reshape(-1)

    _, ms1, _ = m.simulate_onward(y0, p0)
    _, ms2, _ = m.simulate_onward(y0, p0, intv={'ACF': {'Scale': 1, 'Type': 'mod'}})

    # ms = ms[ms.index > 2000]
    # ms.Pop.plot()
    # ms.Pop_RiskLo.plot()
    # ms.Pop_RiskHi.plot()
    ms1.Pop.plot()
    ms2.Pop.plot()

    # ms.PropComorb.plot()

    # ms.LTBI.plot()
    # ms.LTBI_RiskLo.plot()
    # ms.LTBI_RiskHi.plot()

    # ms.RR_inf_comorb.plot()
    # ms.RR_inc_comorb.plot()

    # ms1.Prev_RiskLo.plot()
    # ms1.Prev_RiskHi.plot()
    # ms2.Prev_RiskLo.plot()
    # ms2.Prev_RiskHi.plot()
    # ms.PrPrev_DR.plot()

    # ms.IncR.plot()
    # ms.Prev.plot()
    # ms.MorR.plot()
    # ms.CNR.plot()
    # ms.CNR_Pub.plot()
    # ms.CNR_Pri.plot()
    # ms1.CNR_Pub.plot()
    # ms2.CNR_Pub.plot()

    # ms.PrSp_Asym.plot()
    # ms.PrSp_Sym.plot()
    # ms.PrSym.plot()
    # ms.PrDR_Inc.plot()
    # ms.IncR.plot()
    # ms1.IncR.plot()
    # ms2.IncR.plot()
    # ms.IncR_DS.plot()
    # ms.IncR_DR.plot()
    # ms.PrDR_Inc.plot()
    # ms.IncR_rural.plot()
    # ms.IncR_urban.plot()
    # ms.IncR_slum.plot()

    # ms.IncR_DS.plot()
    # ms.IncR_DR.plot()
    # ms1.IncR_DR.plot()
    # ms2.IncR_DR.plot()

    plt.legend()
    plt.show()
