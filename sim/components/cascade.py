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

    def __call__(self, t, y, pars, intv, calc):
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
        r_acf, r_acf_tp, r_acf_fp, p_dst = 0, 0, 0, 0
        n_nontb = pars['NonTB'] * y.sum(0) if 'NonTB' in pars else 0
        n_tb = y[I.Sym].sum(0) + y[I.ExSym].sum(0)
        n = y.sum()
        if intv is not None and t > intv.T0_Intv:
            r_acf, r_acf_tp, r_acf_fp, p_dst = \
                intv.modify_acf(t, r_acf, r_acf_tp, r_acf_fp, p_dst, pars, n_tb / n, n_nontb / n)

        p_dst_acf = np.array([0, 0, p_dst, p_dst]).reshape((-1, 1))

        calc['r_acf'] = r_acf
        calc['r_acf_fp'] = r_acf_fp

        calc['acf_s'] = r_acf * y[I.Sym]
        calc['acf_c'] = r_acf * y[I.ExSym]

        calc['acf_txf_pub_s'] = r_acf_tp * (1 - p_dst_acf) * y[I.Sym]
        calc['acf_txf_pub_c'] = r_acf_tp * (1 - p_dst_acf) * y[I.ExSym]
        calc['acf_txs_pub_s'] = r_acf_tp * p_dst_acf * y[I.Sym]
        calc['acf_txs_pub_c'] = r_acf_tp * p_dst_acf * y[I.ExSym]

        # ACF non-TB population
        n_nontb = pars['NonTB'] * y.sum(0, keepdims=True) if 'NonTB' in pars else np.zeros(2)
        calc['acf_nontb'] = n_nontb * r_acf
        calc['acf_tx_fp'] = n_nontb * r_acf_fp

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

    def get_trs(self, t, y, pars, intv):
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
        p_txf_pub, p_txf_pri = pars['p_txf_pub'], pars['p_txf_pri']
        p_txs_pub, p_txs_pri = pars['p_txs_pub'], pars['p_txs_pri']

        r_succ_txf = pars['r_succ_txf']
        r_ltfu_txf_pub = pars['r_succ_txf'] * (1 - pars['p_succ_txf_pub']) / pars['p_succ_txf_pub']
        r_ltfu_txf_pri = pars['r_succ_txf'] * (1 - pars['p_succ_txf_pri']) / pars['p_succ_txf_pri']
        r_succ_txs = pars['r_succ_txs']
        r_ltfu_txs_pub = pars['r_succ_txs'] * (1 - pars['p_succ_txs_pub']) / pars['p_succ_txs_pub']
        p_tr_pub = pars['p_tr_pub']

        for sym, cs, txf_pub, txf_pri, txs_pub in zip(I.Sym, I.ExSym, I.Txf_Pub, I.Txf_Pri, I.Txs_Pub):
            ds = sym in I.Sym_DS

            # PCF
            if ds:
                trs += [
                    (sym, txf_pub, r_cs_s * p_entry_pub * p_dx_pub * p_txf_pub, 'pcf_pub'),
                    (cs, txf_pub, r_cs_c * p_entry_pub * p_dx_pub * p_txf_pub, 'pcf_pub'),
                    (sym, cs, r_cs_s * p_entry_pub * (1 - p_dx_pub * p_txf_pub), 'fn0'),
                ]
            else:
                trs += [
                    (sym, txf_pub, r_cs_s * p_entry_pub * (1 - p_dst) * p_dx_pub * p_txf_pub, 'pcf_pub'),
                    (sym, txs_pub, r_cs_s * p_entry_pub * p_dst * p_dx_pub * p_txs_pub, 'pcf_pub'),
                    (cs, txf_pub, r_cs_c * p_entry_pub * (1 - p_dst) * p_dx_pub * p_txf_pub, 'pcf_pub'),
                    (cs, txs_pub, r_cs_c * p_entry_pub * p_dst * p_dx_pub * p_txs_pub, 'pcf_pub'),
                    (sym, cs, r_cs_s * p_entry_pub *
                     (1 - (1 - p_dst) * p_dx_pub * p_txf_pub - p_dst * p_dx_pub * p_txs_pub), 'fn0'),
                ]

            trs += [
                (sym, txf_pri, r_cs_s * p_entry_pri * p_dx_pri * p_txf_pri, 'pcf_pri'),
                (cs, txf_pri, r_cs_c * p_entry_pri * p_dx_pri * p_txf_pri, 'pcf_pri'),
                (sym, cs, r_cs_s * p_entry_pri * (1 - p_dx_pri * p_txf_pri), 'fn0'),
            ]

            # Treatment outcome
            if ds:
                rl, rh = I.RLow_DS, I.RHigh_DS
            else:
                rl, rh = I.RLow_DR, I.RHigh_DR

            if t < self.T0_Rif:
                trs += [
                    (txf_pub, rh, r_ltfu_txf_pub + r_succ_txf, 'tx_ltfu'),
                    (txs_pub, rh, r_ltfu_txs_pub + r_succ_txs, 'tx_ltfu'),
                    (txf_pri, rh, r_ltfu_txf_pri + r_succ_txf, 'tx_ltfu')
                ]
            else:
                trs += [
                    (txf_pub, rl, r_succ_txf, 'tx_succ'),
                    (txf_pub, rh, r_ltfu_txf_pub * (1 - p_tr_pub), 'tx_ltfu'),
                    (txf_pub, txs_pub, r_ltfu_txf_pub * p_tr_pub, 'tx_ltfu'),

                    (txs_pub, rl, r_succ_txs, 'tx_succ'),
                    (txs_pub, rh, r_ltfu_txs_pub, 'tx_ltfu'),
                    (txf_pri, rl, r_succ_txf, 'tx_succ'),
                    (txf_pri, rh, r_ltfu_txf_pri, 'tx_ltfu'),
                ]

        if intv is not None and t >= intv.T0_Intv:
            acf_l, acf_h, rates = self.trs_acf(t, y, pars, intv)
            trs_l = trs + acf_l
            trs_h = trs + acf_h
        else:
            trs_l = trs_h = trs
            rates = dict()

        return trs_l, trs_h, rates

    def trs_acf(self, t, y, pars, intv):
        if intv is None or t < intv.T0_Intv:
            return list(), list()

        I = self.Keys

        # ACF
        r_acf0, r_acf_tp, r_acf_fp, p_dst = 0, 0, 0, 0
        n_nontb = pars['NonTB'] * y.sum(0) if 'NonTB' in pars else 0
        n_tb = y[I.Sym].sum(0) + y[I.ExSym].sum(0)
        n = y.sum()
        p_tb, p_nontb = n_tb / n, n_nontb / n

        r_acf0, r_acf_tp, r_acf_fp, p_dst = intv.modify_acf(t, r_acf0, r_acf_tp, r_acf_fp, p_dst,
                                                            pars, p_tb, p_nontb)

        trs_l = [
            (I.Sym_Sn_DS, I.Txf_Pub_Sn_DS, r_acf_tp[0, 0], 'acf'),
            (I.Sym_Sp_DS, I.Txf_Pub_Sp_DS, r_acf_tp[1, 0], 'acf'),
            (I.Sym_Sn_DR, I.Txf_Pub_Sn_DR, r_acf_tp[0, 0] * (1 - p_dst), 'acf'),
            (I.Sym_Sp_DR, I.Txf_Pub_Sp_DR, r_acf_tp[1, 0] * (1 - p_dst), 'acf'),
            (I.Sym_Sn_DR, I.Txs_Pub_Sn_DR, r_acf_tp[0, 0] * p_dst, 'acf'),
            (I.Sym_Sp_DR, I.Txs_Pub_Sp_DR, r_acf_tp[1, 0] * p_dst, 'acf'),

            (I.ExSym_Sn_DS, I.Txf_Pub_Sn_DS, r_acf_tp[0, 0], 'acf'),
            (I.ExSym_Sp_DS, I.Txf_Pub_Sp_DS, r_acf_tp[1, 0], 'acf'),
            (I.ExSym_Sn_DR, I.Txf_Pub_Sn_DR, r_acf_tp[0, 0] * (1 - p_dst), 'acf'),
            (I.ExSym_Sp_DR, I.Txf_Pub_Sp_DR, r_acf_tp[1, 0] * (1 - p_dst), 'acf'),
            (I.ExSym_Sn_DR, I.Txs_Pub_Sn_DR, r_acf_tp[0, 0] * p_dst, 'acf'),
            (I.ExSym_Sp_DR, I.Txs_Pub_Sp_DR, r_acf_tp[1, 0] * p_dst, 'acf'),
        ]
        trs_h = [
            (I.Sym_Sn_DS, I.Txf_Pub_Sn_DS, r_acf_tp[0, 1], 'acf'),
            (I.Sym_Sp_DS, I.Txf_Pub_Sp_DS, r_acf_tp[1, 1], 'acf'),
            (I.Sym_Sn_DR, I.Txf_Pub_Sn_DR, r_acf_tp[0, 1] * (1 - p_dst), 'acf'),
            (I.Sym_Sp_DR, I.Txf_Pub_Sp_DR, r_acf_tp[1, 1] * (1 - p_dst), 'acf'),
            (I.Sym_Sn_DR, I.Txs_Pub_Sn_DR, r_acf_tp[0, 1] * p_dst, 'acf'),
            (I.Sym_Sp_DR, I.Txs_Pub_Sp_DR, r_acf_tp[1, 1] * p_dst, 'acf'),

            (I.ExSym_Sn_DS, I.Txf_Pub_Sn_DS, r_acf_tp[0, 1], 'acf'),
            (I.ExSym_Sp_DS, I.Txf_Pub_Sp_DS, r_acf_tp[1, 1], 'acf'),
            (I.ExSym_Sn_DR, I.Txf_Pub_Sn_DR, r_acf_tp[0, 1] * (1 - p_dst), 'acf'),
            (I.ExSym_Sp_DR, I.Txf_Pub_Sp_DR, r_acf_tp[1, 1] * (1 - p_dst), 'acf'),
            (I.ExSym_Sn_DR, I.Txs_Pub_Sn_DR, r_acf_tp[0, 1] * p_dst, 'acf'),
            (I.ExSym_Sp_DR, I.Txs_Pub_Sp_DR, r_acf_tp[1, 1] * p_dst, 'acf'),
        ]

        return trs_l, trs_h, {'r_acf0': r_acf0, 'r_acf_tp': r_acf_tp, 'r_acf_fp': r_acf_fp, 'p_dst': p_dst}

    def calc_dy(self, t, y, pars, intv):
        trs_l, trs_h, _ = self.get_trs(t, y, pars, intv)
        dy = np.zeros_like(y)
        dy[:, 0] = calc_dy(y[:, 0], trs_l)
        dy[:, 1] = calc_dy(y[:, 1], trs_h)
        return dy

    def measure(self, t, y, pars, intv, mea):
        I = self.Keys

        trs_l, trs_h, rates = self.get_trs(t, y, pars, intv)

        fil = lambda x: x[3] == 'pcf_pub'
        det_pub = np.array([extract_tr(y[:, 0], trs_l, fil), extract_tr(y[:, 1], trs_h, fil)])

        fil = lambda x: x[3] == 'pcf_pri'
        det_pri = np.array([extract_tr(y[:, 0], trs_l, fil), extract_tr(y[:, 1], trs_h, fil)])

        fil = lambda x: x[3] == 'acf'
        det_acf = np.array([extract_tr(y[:, 0], trs_l, fil), extract_tr(y[:, 1], trs_h, fil)])

        det = det_pub + det_pri + det_acf

        ns = y.sum(0)
        n = ns.sum()

        mea['TP'] = det.sum() / n
        mea['TP_Pub'] = det_pub.sum() / n
        mea['TP_Pri'] = det_pri.sum() / n
        mea['TP_Pcf'] = mea['TP_Pub'] + mea['TP_Pri']
        mea['TP_Acf'] = det_acf.sum() / n
        mea['N_Pub_Detected'] = det_pub.sum()
        mea['N_Pri_Detected'] = det_pri.sum()

        if 'r_acf_fp' in rates and 'NonTB' in pars:
            n_tb = y[I.Sym].sum(0) + y[I.ExSym].sum(0)
            n_nontb = pars['NonTB'] * y.sum(0)
            fp_pcf = n_nontb * pars['r_cs_s'] * (1 - pars['spec0'])
            fp_acf = n_nontb * rates['r_acf_fp']
            mea['FP_Pcf'] = fp_pcf.sum() / n
            mea['FP_Acf'] = fp_acf.sum() / n
            mea['FP'] = mea['FP_Pcf'] + mea['FP_Acf']

            mea['PPV_Pcf'] = mea['TP_Pcf'] / (mea['TP_Pcf'] + mea['FP_Pcf'] + 1e-10)
            mea['PPV_Acf'] = mea['TP_Acf'] / (mea['TP_Acf'] + mea['FP_Acf'] + 1e-10)

            mea['N_ACF_TB_Reached'] = (n_tb * rates['r_acf0']).sum()
            mea['N_ACF_NonTB_Reached'] = (n_nontb * rates['r_acf0']).sum()
            mea['N_ACF_Reached'] = mea['N_ACF_TB_Reached'] + mea['N_ACF_NonTB_Reached']

        else:
            mea['FP'] = mea['FP_Pcf'] = mea['FP_Acf'] = 0
            mea['PPV'] = mea['PPV_Pcf'] = mea['PPV_Acf'] = 1
            mea['N_ACF_Reached'] = mea['N_ACF_TB_Reached'] = mea['N_ACF_NonTB_Reached'] = 0

        mea['R_ACF_Reached'] = mea['N_ACF_Reached'] / n
