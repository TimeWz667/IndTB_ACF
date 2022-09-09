import numpy as np
from sim.components.base import Process
from sim.util import calc_dy, extract_tr

__author__ = 'Chu-Chang Ku'
__all__ = ['Cascade']


class Cascade(Process):
    def __init__(self, keys):
        Process.__init__(self, keys)
        self.T0_Rif = 1970
        self.T0_DOTS = 1990
        self.T1_DOTS = 2007

    def get_trs(self, t, y, pars):
        I = self.Keys

        trs = list()
        # PCF
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
        p_txi = pars['p_txi']

        r_succ_fl_pub_ds = pars['r_succ_fl_pub_ds']
        r_succ_fl_pri_ds = pars['r_succ_fl_pri_ds']
        r_succ_fl_pub_dr = pars['r_succ_fl_pub_dr']
        r_succ_fl_pri_dr = pars['r_succ_fl_pri_dr']
        r_succ_sl_pub_dr = pars['r_succ_sl_pub_dr']

        r_ltfu_fl_pub_ds = pars['r_ltfu_fl_pub_ds']
        r_ltfu_fl_pri_ds = pars['r_ltfu_fl_pri_ds']
        r_ltfu_fl_pub_dr = pars['r_ltfu_fl_pub_dr']
        r_ltfu_fl_pri_dr = pars['r_ltfu_fl_pri_dr']
        r_ltfu_sl_pub_dr = pars['r_ltfu_sl_pub_dr']

        p_tr_pub = pars['p_tr_pub']

        for sym, cs, txf_pub, txf_pri, txs_pub in zip(I.Sym, I.ExSym, I.Txf_Pub, I.Txf_Pri, I.Txs_Pub):
            ds = sym is I.Sym_DS

            # PCF
            if ds:
                trs += [
                    (sym, txf_pub, r_cs_s * p_entry_pub * p_dx_pub * p_txi, 'pcf_pub'),
                    (cs, txf_pub, r_cs_c * p_entry_pub * p_dx_pub * p_txi, 'pcf_pub'),
                    (sym, cs, r_cs_s * p_entry_pub * (1 - p_dx_pub * p_txi), 'fn0'),
                ]
            else:
                trs += [
                    (sym, txf_pub, r_cs_s * p_entry_pub * (1 - p_dst) * p_dx_pub * p_txi, 'pcf_pub'),
                    (sym, txs_pub, r_cs_s * p_entry_pub * p_dst * p_dx_pub * p_txi, 'pcf_pub'),
                    (cs, txf_pub, r_cs_c * p_entry_pub * (1 - p_dst) * p_dx_pub * p_txi, 'pcf_pub'),
                    (cs, txs_pub, r_cs_c * p_entry_pub * p_dst * p_dx_pub * p_txi, 'pcf_pub'),
                    (sym, cs, r_cs_s * p_entry_pub *
                     (1 - (1 - p_dst) * p_dx_pub * p_txi - p_dst * p_dx_pub * p_txi), 'fn0'),
                ]

            trs += [
                (sym, txf_pri, r_cs_s * p_entry_pri * p_dx_pri * p_txi, 'pcf_pri'),
                (cs, txf_pri, r_cs_c * p_entry_pri * p_dx_pri * p_txi, 'pcf_pri'),
                (sym, cs, r_cs_s * p_entry_pri * (1 - p_dx_pri * p_txi), 'fn0'),
            ]

            # Treatment outcome
            if ds:
                rl, rh = I.RLow_DS, I.RHigh_DS

                if t < self.T0_Rif:
                    trs += [
                        (txf_pub, rh, r_ltfu_fl_pub_ds + r_succ_fl_pub_ds, 'tx_ltfu'),
                        (txf_pri, rh, r_ltfu_fl_pri_ds + r_succ_fl_pri_ds, 'tx_ltfu')
                    ]
                else:
                    trs += [
                        (txf_pub, rl, r_succ_fl_pub_ds, 'tx_succ'),
                        (txf_pub, rh, r_ltfu_fl_pub_ds, 'tx_ltfu'),
                        (txf_pri, rl, r_succ_fl_pri_ds, 'tx_succ'),
                        (txf_pri, rh, r_ltfu_fl_pri_ds, 'tx_ltfu'),
                    ]
            else:
                rl, rh = I.RLow_DR, I.RHigh_DR

                if t < self.T0_Rif:
                    trs += [
                        (txf_pub, rh, r_ltfu_fl_pub_dr + r_succ_fl_pub_dr, 'tx_ltfu'),
                        (txs_pub, rh, r_ltfu_sl_pub_dr + r_succ_sl_pub_dr, 'tx_ltfu'),
                        (txf_pri, rh, r_ltfu_fl_pri_dr + r_succ_fl_pri_dr, 'tx_ltfu')
                    ]
                else:
                    trs += [
                        (txf_pub, rl, r_succ_fl_pub_dr, 'tx_succ'),
                        (txf_pub, rh, r_ltfu_fl_pub_dr * (1 - p_tr_pub), 'tx_ltfu'),
                        (txf_pub, txs_pub, r_ltfu_fl_pub_dr * p_tr_pub, 'tx_ltfu'),

                        (txs_pub, rl, r_succ_sl_pub_dr, 'tx_succ'),
                        (txs_pub, rh, r_ltfu_sl_pub_dr, 'tx_ltfu'),
                        (txf_pri, rl, r_succ_fl_pri_dr, 'tx_succ'),
                        (txf_pri, rh, r_ltfu_fl_pri_dr, 'tx_ltfu'),
                    ]

        return trs, trs

    def calc_dy(self, t, y, pars):
        trs_l, trs_h = self.get_trs(t, y, pars)
        dy = np.zeros_like(y)
        dy[:, 0] = calc_dy(y[:, 0], trs_l)
        dy[:, 1] = calc_dy(y[:, 1], trs_h)
        return dy

    def measure(self, t, y, pars, mea):
        I = self.Keys

        trs_l, trs_h = self.get_trs(t, y, pars)

        fil = lambda x: x[3] == 'pcf_pub'
        det_pub = np.array([extract_tr(y[:, 0], trs_l, fil), extract_tr(y[:, 1], trs_h, fil)])

        fil = lambda x: x[3] == 'pcf_pri'
        det_pri = np.array([extract_tr(y[:, 0], trs_l, fil), extract_tr(y[:, 1], trs_h, fil)])

        det = det_pub + det_pri

        fil = lambda x: x[0] in I.Infectious_DR and x[3] == 'pcf_pub'
        det_dr = np.array([extract_tr(y[:, 0], trs_l, fil), extract_tr(y[:, 1], trs_h, fil)])

        ns = y.sum(0)
        n = ns.sum()

        mea['TP'] = det.sum() / n
        mea['TP_Pub'] = det_pub.sum() / n
        mea['TP_Pri'] = det_pri.sum() / n
        mea['TP_Pcf'] = mea['TP_Pub'] + mea['TP_Pri']
        mea['N_Pub_Detected'] = det_pub.sum()
        mea['N_Pri_Detected'] = det_pri.sum()
        mea['PrDR_CNR'] = det_dr.sum() / max(det_pub.sum(), 1e-10)
