import numpy as np
from sim.components.base import Process

__author__ = 'Chu-Chang Ku'
__all__ = ['Cascade']


class Cascade(Process):
    def __call__(self, t, y, pars, calc):
        I = self.Keys

        r_cs_s, r_cs_c = pars['r_cs_s'], pars['r_cs_c']

        r_cs_s, r_cs_c = self.Intervention.modify_cs(t, r_cs_s, r_cs_c)

        p_dx, p_txi = pars['p_dx'], pars['p_txi']

        p_dx = self.Intervention.modify_dx(t, p_dx)

        p_tp = p_dx * p_txi

        r_det_s = r_cs_s * p_tp
        r_fn_s = r_cs_s - r_det_s

        r_det_c = r_cs_c * p_tp
        r_fn_c = r_cs_c - r_det_c

        calc['det_s'] = r_det_s * y[I.Sym]
        calc['fn_s'] = r_fn_s * y[I.Sym]
        calc['det_c'] = r_det_c * y[I.ExSym]
        calc['fn_c'] = r_fn_c * y[I.ExSym]

        r_succ_tx, r_ltfu_tx = pars['r_succ_tx'], pars['r_ltfu_tx']
        calc['tx_succ_tx'] = r_succ_tx * y[I.Tx]
        calc['tx_ltfu_tx'] = r_ltfu_tx * y[I.Tx]

    def measure(self, t, y, pars, calc, mea):
        I = self.Keys

        ns = y.sum(0)
        n = ns.sum()

        det = calc['det_s'] + calc['det_c']

        mea['CNR'] = det.sum() / n

        for i, strata in enumerate(I.Tag_Strata):
            n = max(ns[i], 1e-15)

            mea[f'CNR_{strata}'] = det[i] / n
