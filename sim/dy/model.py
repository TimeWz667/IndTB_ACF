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

        pars['sus'] = sus = np.zeros((I.N_State_TB, I.N_State_Strata))
        sus[I.U] = 1
        sus[I.SLat] = pars['rr_sus_slat']
        sus[I.RLow] = pars['rr_sus_rec']
        sus[I.RHigh] = pars['rr_sus_rec']
        sus[I.RSt] = pars['rr_sus_rec']

        pars['trans'] = trans = np.zeros((I.N_State_TB, I.N_State_Strata))
        trans[I.Asym] = pars['rr_inf_asym']
        trans[I.Sym] = 1
        trans[I.ExSym] = pars['rr_inf_cs']

        pars['mixing'] = np.ones((I.N_State_Strata, I.N_State_Strata))

        return pars

    def get_y0(self, pars):
        y0 = np.zeros((I.N_State_TB, I.N_State_Strata))

        n0 = np.array([self.Inputs['N0'], 0])

        y0[I.Sym] = 1e-2 * n0
        y0[I.SLat] = 0.4 * n0
        y0[I.U] = n0 - y0.sum(0)
        return y0

    def collect_calc(self, t, y, pars):
        t = max(t, self.Year0)

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

        dy -= calc['infection']
        dy[I.FLat] += calc['infection'].sum(0)

        dy[I.FLat] -= calc['act']
        dy[I.SLat] -= calc['react']
        dy[I.RLow] -= calc['rel_tc']
        dy[I.RHigh] -= calc['rel_td']
        dy[I.RSt] -= calc['rel_st']

        dy[I.Asym] += calc['inc']

        # Progression
        dy[I.FLat] += - calc['stab_fl']
        dy[I.SLat] += calc['stab_fl']

        dy[I.RLow] -= calc['stab_tc']
        dy[I.RHigh] -= calc['stab_td']
        dy[I.RSt] += calc['stab_tc'] + calc['stab_td']

        dy[I.Asym] -= calc['sc_a'] + calc['sym_onset']
        dy[I.Sym] += calc['sym_onset'] - calc['sc_s']
        dy[I.ExSym] -= calc['sc_c']
        dy[I.RSt] += calc['sc_a'] + calc['sc_s'] + calc['sc_c']

        # Dx
        det_s = calc['det_s']
        fn_s = calc['fn_s']
        det_c = calc['det_c']

        dy[I.Sym] -= det_s + fn_s
        dy[I.ExSym] += fn_s - det_c
        dy[I.Tx] += det_s + det_c

        # Tx
        tc, td = calc['tx_succ_tx'], calc['tx_ltfu_tx']

        dy[I.Tx] -= tc + td
        dy[I.RLow] += tc
        dy[I.RHigh] += td

        # Self-clearance
        dy[I.SLat] -= calc['clear_sl']
        dy[I.RSt] -= calc['clear_rst']
        dy[I.U] += (calc['clear_sl'] + calc['clear_rst'])

        # Demography
        dy[I.U] += calc['births']
        dy -= calc['deaths'] + calc['deaths_tb']

        if t <= self.Year0:
            ns = y.sum(0, keepdims=True)
            ns[ns == 0] = 1e-10
            dy -= y / ns * dy.sum(0, keepdims=True)

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
        ntb = y.reshape((I.N_State_TB, I.N_State_Strata))[I.PTB].sum()
        return ntb - 0.5

    dfe.terminal = True
    dfe.direction = -1

    def simulate(self, p):
        p0 = p
        if 'sus' not in p:
            p = self.update_parameters(p)

        ys, ms, msg = simulate(self,
                               pars=p,
                               t_warmup=300,
                               t_out=np.linspace(2000, 2020, int(20 * 2) + 1),
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

    m = Model(inputs, year0=2000)
    m.CascadeOutput = False

    sc = get_bn()
    p0 = sample(sc)

    ys, ms, msg = m.simulate(p0)
    ys = ys.y.T[-1]
    _, ms1, _ = m.simulate_onward(ys, p0)
    _, ms2, _ = m.simulate_onward(ys, p0, intv={'PPM': {'Scale': 1}, 'UACF': {'Scale': 0.2}})

    ms = pd.concat([ms, ms1.iloc[1:]])
    # ms.Pop.plot()
    # ms.Pop_rural.plot()
    # ms.Pop_urban.plot()
    # ms.Pop_slum.plot()

    # ms.LTBI.plot()
    # ms.LTBI_rural.plot()
    # ms.LTBI_urban.plot()
    # ms.LTBI_slum.plot()

    # ms.Prev.plot()
    # ms.Prev_DS.plot()
    # ms.Prev_DR.plot()
    # ms.PrPrev_DR.plot()

    ms.IncR.plot()
    ms.Prev.plot()
    ms.MorR.plot()
    # ms1.IncR.plot()
    # # ms2.IncR.plot()
    # ms.IncR_DS.plot()
    # ms.IncR_DR.plot()
    # ms.PrDR_Inc.plot()
    # ms.IncR_rural.plot()
    # ms.IncR_urban.plot()
    # ms.IncR_slum.plot()

    ms.CNR.plot()

    plt.legend()
    plt.show()
