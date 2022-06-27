import numpy as np
from sim.components.base import Process

__author__ = 'Chu-Chang Ku'
__all__ = ['Cascade']


class Cascade(Process):
    def __init__(self, keys, intv=None):
        Process.__init__(self, keys, intv)
        self.T0_Rif = 1970
        self.T0_DOTS = 1990
        self.T1_DOTS = 2007

    def __call__(self, t, y, pars, calc):
        I = self.Keys

        r_cs_s, r_cs_c = pars['r_cs_s'], pars['r_cs_c']

        p_entry_pri = pars['p_entry_pri']
        p_entry_pub = 1 - p_entry_pri

        if t < self.T0_DOTS:
            p_entry_pub, p_entry_pri = 0, 1
        elif t < self.T1_DOTS:
            p_entry_pub *= (t - self.T0_DOTS) / (self.T1_DOTS - self.T0_DOTS)
            p_entry_pri = 1 - p_entry_pub

        p_dx_pub, p_dx_pri = pars['p_dx_pub'], pars['p_dx_pri']

        p_dst = pars['p_dst_pcf']
        p_dst_pub = np.array([0, 0, p_dst, p_dst]).reshape((-1, 1))

        p_txf_pub, p_txf_pri = pars['p_txf_pub'], pars['p_txf_pri']
        p_txs_pub, p_txs_pri = pars['p_txs_pub'], pars['p_txs_pri']

        p_txi = p_dx_pub * (1 - p_dst_pub) * p_txf_pub
        calc['det_txf_pub_s'] = r_cs_s * p_entry_pub * p_txi * y[I.Sym]
        calc['det_txf_pub_c'] = r_cs_c * p_entry_pub * p_txi * y[I.ExSym]

        p_txi = p_dx_pub * p_dst_pub * p_txs_pub
        calc['det_txs_pub_s'] = r_cs_s * p_entry_pub * p_txi * y[I.Sym]
        calc['det_txs_pub_c'] = r_cs_c * p_entry_pub * p_txi * y[I.ExSym]

        p_txi = p_dx_pri * p_txf_pri
        calc['det_txf_pri_s'] = r_cs_s * p_entry_pri * p_txi * y[I.Sym]
        calc['det_txf_pri_c'] = r_cs_c * p_entry_pri * p_txi * y[I.ExSym]

        p_txi = p_dx_pub * p_txf_pub
        calc['fn_pub_s'] = r_cs_s * p_entry_pub * (1 - p_txi) * y[I.Sym]
        calc['fn_pub_c'] = r_cs_c * p_entry_pub * (1 - p_txi) * y[I.ExSym]

        p_txi = p_dx_pri * p_txf_pri
        calc['fn_pri_s'] = r_cs_s * p_entry_pri * (1 - p_txi) * y[I.Sym]
        calc['fn_pri_c'] = r_cs_c * p_entry_pri * (1 - p_txi) * y[I.ExSym]

        # ACF

        r_acf, r_acf_det, p_dst = 0, 0, 0
        if t > self.Intervention.T0_Intv:
            r_acf, r_acf_det, p_dst = self.Intervention.modify_acf(t, r_acf, r_acf_det, p_dst, pars)

        p_dst_acf = np.array([0, 0, p_dst, p_dst]).reshape((-1, 1))

        calc['acf_s'] = r_acf * y[I.Sym]
        calc['acf_c'] = r_acf * y[I.ExSym]

        calc['acf_txf_pub_s'] = r_acf_det * (1 - p_dst_acf) * y[I.Sym]
        calc['acf_txf_pub_c'] = r_acf_det * (1 - p_dst_acf) * y[I.ExSym]
        calc['acf_txs_pub_s'] = r_acf_det * p_dst_acf * y[I.Sym]
        calc['acf_txs_pub_c'] = r_acf_det * p_dst_acf * y[I.ExSym]

        # Tx
        r_succ_txf_pub, r_succ_txf_pri = pars['r_succ_txf_pub'], pars['r_succ_txf_pri']
        r_succ_txs_pub, r_succ_txs_pri = pars['r_succ_txs_pub'], pars['r_succ_txs_pri']
        r_ltfu_txf_pub, r_ltfu_txf_pri = pars['r_ltfu_txf_pub'], pars['r_ltfu_txf_pri']
        r_ltfu_txs_pub, r_ltfu_txs_pri = pars['r_ltfu_txs_pub'], pars['r_ltfu_txs_pri']

        if t < self.T0_Rif:
            calc['tx_succ_txf_pub'] = r_succ_txf_pub * y[I.Txf_Pub]
            calc['tx_succ_txf_pri'] = r_succ_txf_pri * y[I.Txf_Pri]
            calc['tx_succ_txs_pub'] = r_succ_txs_pub * y[I.Txs_Pub]
            calc['tx_succ_txs_pri'] = r_succ_txs_pri * y[I.Txs_Pri]

            calc['tx_ltfu_txf_pub'] = r_ltfu_txf_pub * (1 - pars['p_tr_pub']) * y[I.Txf_Pub]
            calc['tx_ltfu_txf_pri'] = r_ltfu_txf_pri * y[I.Txf_Pri]
            calc['tx_ltfu_txs_pub'] = r_ltfu_txs_pub * y[I.Txs_Pub]
            calc['tx_ltfu_txs_pri'] = r_ltfu_txs_pri * y[I.Txs_Pri]

            calc['tx_switch_pub'] = r_ltfu_txf_pub * pars['p_tr_pub'] * y[I.Txf_Pub]
        else:
            calc['tx_succ_txf_pub'] = 0 * y[I.Txf_Pub]
            calc['tx_succ_txf_pri'] = 0 * y[I.Txf_Pri]
            calc['tx_succ_txs_pub'] = 0 * y[I.Txs_Pub]
            calc['tx_succ_txs_pri'] = 0 * y[I.Txs_Pri]

            calc['tx_ltfu_txf_pub'] = (r_succ_txf_pub + r_ltfu_txf_pub) * y[I.Txf_Pub]
            calc['tx_ltfu_txf_pri'] = (r_succ_txf_pri + r_ltfu_txf_pri) * y[I.Txf_Pri]
            calc['tx_ltfu_txs_pub'] = (r_succ_txs_pub + r_ltfu_txs_pub) * y[I.Txs_Pub]
            calc['tx_ltfu_txs_pri'] = (r_succ_txs_pri + r_ltfu_txs_pri) * y[I.Txs_Pri]

            calc['tx_switch_pub'] = 0 * y[I.Txf_Pub]

    def measure(self, t, y, pars, calc, mea):
        I = self.Keys

        ns = y.sum(0)
        n = ns.sum()

        det_pub = calc['det_txf_pub_s'] + calc['det_txf_pub_c'] + calc['det_txs_pub_s'] + calc['det_txs_pub_c']
        det_pri = calc['det_txf_pri_s'] + calc['det_txf_pri_c']
        det_acf = calc['acf_txf_pub_s'] + calc['acf_txf_pub_c'] + calc['acf_txs_pub_s'] + calc['acf_txs_pub_c']

        det = det_pub + det_pri + det_acf

        mea['CNR'] = det.sum() / n
        mea['CNR_Pub'] = det_pub.sum() / n
        mea['CNR_Pri'] = det_pri.sum() / n
        mea['CNR_Acf'] = det_acf.sum() / n
        mea['CNR_DS'] = det[2:].sum() / n
        mea['CNR_DR'] = det[:2].sum() / n

        for i, strata in enumerate(I.Tag_Strata):
            n = max(ns[i], 1e-15)
            mea[f'CNR_{strata}'] = det[:, i].sum() / n
            mea[f'CNR_DS_{strata}'] = det[2:, i].sum() / n
            mea[f'CNR_DR_{strata}'] = det[:2, i].sum() / n
            mea[f'CNR_Acf_{strata}'] = det_acf[:, i].sum() / n

        mea['N_Pub_Detected'] = det_pub.sum()
        mea['N_Pri_Detected'] = det_pri.sum()
        mea['N_ACF_Detected'] = det_acf.sum()
        mea['N_ACF_Reached'] = (calc['acf_s'] + calc['acf_c']).sum()
