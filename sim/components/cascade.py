import numpy as np
from sim.components.base import Process
from sim.util import calc_dy

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

    def calc_dy(self, t, y, pars, intv):
        I = self.Keys
        calc = dict()
        self(t, y, pars, intv, calc)
        dy = np.zeros_like(y)

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

        return dy

    def calc_dy2(self, t, y, pars, intv):
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

        trs += [
            (I.Sym_Sn_DS, I.Txf_Pub_Sn_DS, r_cs_s * p_entry_pub * p_dx_pub * p_txf_pub),
            (I.Sym_Sp_DS, I.Txf_Pub_Sp_DS, r_cs_s * p_entry_pub * p_dx_pub * p_txf_pub),
            (I.Sym_Sn_DR, I.Txf_Pub_Sn_DR, r_cs_s * p_entry_pub * p_dx_pub * (1 - p_dst) * p_txf_pub),
            (I.Sym_Sp_DR, I.Txf_Pub_Sp_DR, r_cs_s * p_entry_pub * p_dx_pub * (1 - p_dst) * p_txf_pub),
            (I.Sym_Sn_DR, I.Txs_Pub_Sn_DR, r_cs_s * p_entry_pub * p_dx_pub * p_dst * p_txs_pub),
            (I.Sym_Sp_DR, I.Txs_Pub_Sp_DR, r_cs_s * p_entry_pub * p_dx_pub * p_dst * p_txs_pub),

            (I.Sym_Sn_DS, I.Txf_Pri_Sn_DS, r_cs_s * p_entry_pri * p_dx_pri * p_txf_pri),
            (I.Sym_Sp_DS, I.Txf_Pri_Sp_DS, r_cs_s * p_entry_pri * p_dx_pri * p_txf_pri),
            (I.Sym_Sn_DR, I.Txf_Pri_Sn_DR, r_cs_s * p_entry_pri * p_dx_pri * p_txf_pri),
            (I.Sym_Sp_DR, I.Txf_Pri_Sp_DR, r_cs_s * p_entry_pri * p_dx_pri * p_txf_pri),

            (I.ExSym_Sn_DS, I.Txf_Pub_Sn_DS, r_cs_c * p_entry_pub * p_dx_pub * p_txf_pub),
            (I.ExSym_Sp_DS, I.Txf_Pub_Sp_DS, r_cs_c * p_entry_pub * p_dx_pub * p_txf_pub),
            (I.ExSym_Sn_DR, I.Txf_Pub_Sn_DR, r_cs_c * p_entry_pub * p_dx_pub * (1 - p_dst) * p_txf_pub),
            (I.ExSym_Sp_DR, I.Txf_Pub_Sp_DR, r_cs_c * p_entry_pub * p_dx_pub * (1 - p_dst) * p_txf_pub),
            (I.ExSym_Sn_DR, I.Txs_Pub_Sn_DR, r_cs_c * p_entry_pub * p_dx_pub * p_dst * p_txs_pub),
            (I.ExSym_Sp_DR, I.Txs_Pub_Sp_DR, r_cs_c * p_entry_pub * p_dx_pub * p_dst * p_txs_pub),

            (I.ExSym_Sn_DS, I.Txf_Pri_Sn_DS, r_cs_c * p_entry_pri * p_dx_pri * p_txf_pri),
            (I.ExSym_Sp_DS, I.Txf_Pri_Sp_DS, r_cs_c * p_entry_pri * p_dx_pri * p_txf_pri),
            (I.ExSym_Sn_DR, I.Txf_Pri_Sn_DR, r_cs_c * p_entry_pri * p_dx_pri * p_txf_pri),
            (I.ExSym_Sp_DR, I.Txf_Pri_Sp_DR, r_cs_c * p_entry_pri * p_dx_pri * p_txf_pri),
        ]

        trs += [
            (I.Sym_Sn_DS, I.ExSym_Sn_DS, r_cs_s * p_entry_pub * (1 - p_dx_pub * p_txf_pub)),
            (I.Sym_Sp_DS, I.ExSym_Sp_DS, r_cs_s * p_entry_pub * (1 - p_dx_pub * p_txf_pub)),
            (I.Sym_Sn_DR, I.ExSym_Sn_DR, r_cs_s * p_entry_pub * (1 - p_dst) * (1 - p_dx_pub * p_txf_pub)),
            (I.Sym_Sp_DR, I.ExSym_Sp_DR, r_cs_s * p_entry_pub * (1 - p_dst) * (1 - p_dx_pub * p_txf_pub)),
            (I.Sym_Sn_DR, I.ExSym_Sn_DR, r_cs_s * p_entry_pub * p_dst * (1 - p_dx_pub * p_txs_pub)),
            (I.Sym_Sp_DR, I.ExSym_Sp_DR, r_cs_s * p_entry_pub * p_dst * (1 - p_dx_pub * p_txs_pub)),

            (I.Sym_Sn_DS, I.ExSym_Sn_DS, r_cs_s * p_entry_pri * (1 - p_dx_pri * p_txf_pri)),
            (I.Sym_Sp_DS, I.ExSym_Sp_DS, r_cs_s * p_entry_pri * (1 - p_dx_pri * p_txf_pri)),
            (I.Sym_Sn_DR, I.ExSym_Sn_DR, r_cs_s * p_entry_pri * (1 - p_dx_pri * p_txf_pri)),
            (I.Sym_Sp_DR, I.ExSym_Sp_DR, r_cs_s * p_entry_pri * (1 - p_dx_pri * p_txf_pri))
        ]

        r_succ_txf = pars['r_succ_txf']
        r_ltfu_txf_pub = pars['r_succ_txf'] * (1 - pars['p_succ_txf_pub']) / pars['p_succ_txf_pub']
        r_ltfu_txf_pri = pars['r_succ_txf'] * (1 - pars['p_succ_txf_pri']) / pars['p_succ_txf_pri']
        r_succ_txs = pars['r_succ_txs']
        r_ltfu_txs_pub = pars['r_succ_txs'] * (1 - pars['p_succ_txs_pub']) / pars['p_succ_txs_pub']

        if t < self.T0_Rif:
            trs += [
                (I.Txf_Pub_Sn_DS, I.RLow_DS, r_succ_txf),
                (I.Txf_Pub_Sp_DS, I.RLow_DS, r_succ_txf),
                (I.Txf_Pub_Sn_DR, I.RLow_DR, r_succ_txf),
                (I.Txf_Pub_Sp_DR, I.RLow_DR, r_succ_txf),
                (I.Txs_Pub_Sn_DR, I.RLow_DS, r_succ_txs),
                (I.Txs_Pub_Sp_DR, I.RLow_DS, r_succ_txs),

                (I.Txf_Pub_Sn_DS, I.RHigh_DS, r_ltfu_txf_pub * (1 - pars['p_tr_pub'])),
                (I.Txf_Pub_Sp_DS, I.RHigh_DS, r_ltfu_txf_pub * (1 - pars['p_tr_pub'])),
                (I.Txf_Pub_Sn_DR, I.RHigh_DR, r_ltfu_txf_pub * (1 - pars['p_tr_pub'])),
                (I.Txf_Pub_Sp_DR, I.RHigh_DR, r_ltfu_txf_pub * (1 - pars['p_tr_pub'])),
                (I.Txs_Pub_Sn_DR, I.RHigh_DR, r_ltfu_txs_pub),
                (I.Txs_Pub_Sp_DR, I.RHigh_DR, r_ltfu_txs_pub),

                (I.Txf_Pri_Sn_DS, I.RLow_DS, r_succ_txf),
                (I.Txf_Pri_Sp_DS, I.RLow_DS, r_succ_txf),
                (I.Txf_Pri_Sn_DR, I.RLow_DR, r_succ_txf),
                (I.Txf_Pri_Sp_DR, I.RLow_DR, r_succ_txf),

                (I.Txf_Pri_Sn_DS, I.RHigh_DS, r_ltfu_txf_pri),
                (I.Txf_Pri_Sp_DS, I.RHigh_DS, r_ltfu_txf_pri),
                (I.Txf_Pri_Sn_DR, I.RHigh_DR, r_ltfu_txf_pri),
                (I.Txf_Pri_Sp_DR, I.RHigh_DR, r_ltfu_txf_pri),

                (I.Txf_Pub_Sn_DS, I.Txs_Pub_Sn_DS, r_ltfu_txf_pub * pars['p_tr_pub']),
                (I.Txf_Pub_Sp_DS, I.Txs_Pub_Sp_DS, r_ltfu_txf_pub * pars['p_tr_pub']),
                (I.Txf_Pub_Sn_DR, I.Txs_Pub_Sn_DR, r_ltfu_txf_pub * pars['p_tr_pub']),
                (I.Txf_Pub_Sp_DR, I.Txs_Pub_Sp_DR, r_ltfu_txf_pub * pars['p_tr_pub']),
            ]
        else:
            trs += [
                (I.Txf_Pub_Sn_DS, I.RHigh_DS, r_ltfu_txf_pub + r_succ_txf),
                (I.Txf_Pub_Sp_DS, I.RHigh_DS, r_ltfu_txf_pub + r_succ_txf),
                (I.Txf_Pub_Sn_DR, I.RHigh_DR, r_ltfu_txf_pub + r_succ_txf),
                (I.Txf_Pub_Sp_DR, I.RHigh_DR, r_ltfu_txf_pub + r_succ_txf),
                (I.Txs_Pub_Sn_DR, I.RHigh_DR, r_ltfu_txs_pub + r_succ_txs),
                (I.Txs_Pub_Sp_DR, I.RHigh_DR, r_ltfu_txs_pub + r_succ_txs),

                (I.Txf_Pri_Sn_DS, I.RHigh_DS, r_ltfu_txf_pri + r_succ_txf),
                (I.Txf_Pri_Sp_DS, I.RHigh_DS, r_ltfu_txf_pri + r_succ_txf),
                (I.Txf_Pri_Sn_DR, I.RHigh_DR, r_ltfu_txf_pri + r_succ_txf),
                (I.Txf_Pri_Sp_DR, I.RHigh_DR, r_ltfu_txf_pri + r_succ_txf),
            ]

        if intv is not None and t > intv.T0_Intv:
            acf_l, acf_h = self.trs_acf(t, y, pars, intv)
            trs_l = trs + acf_l
            trs_h = trs + acf_h
        else:
            trs_l = trs_h = trs

        dy = np.zeros_like(y)
        dy[:, 0] = calc_dy(y[:, 0],
                           frs=np.array([fr for fr, _, _ in trs_l], int),
                           tos=np.array([to for _, to, _ in trs_l], int),
                           rates=np.array([rate for _, _, rate in trs_l])
                           )

        dy[:, 1] = calc_dy(y[:, 1],
                           frs=np.array([fr for fr, _, _ in trs_h], int),
                           tos=np.array([to for _, to, _ in trs_h], int),
                           rates=np.array([rate for _, _, rate in trs_h])
                           )
        return dy

    def trs_acf(self, t, y, pars, intv):
        I = self.Keys

        # ACF
        r_acf, r_acf_tp, r_acf_fp, p_dst = 0, 0, 0, 0
        n_nontb = pars['NonTB'] * y.sum(0) if 'NonTB' in pars else 0
        n_tb = y[I.Sym].sum(0) + y[I.ExSym].sum(0)
        n = y.sum()
        if intv is not None and t > intv.T0_Intv:
            r_acf, r_acf_tp, r_acf_fp, p_dst = \
                intv.modify_acf(t, r_acf, r_acf_tp, r_acf_fp, p_dst, pars, n_tb / n, n_nontb / n)

        trs_l = [
            (I.Sym_Sn_DS, I.Txf_Pub_Sn_DS, r_acf_tp[0, 0]),
            (I.Sym_Sp_DS, I.Txf_Pub_Sp_DS, r_acf_tp[1, 0]),
            (I.Sym_Sn_DR, I.Txf_Pub_Sn_DR, r_acf_tp[2, 0] * (1 - p_dst)),
            (I.Sym_Sp_DR, I.Txf_Pub_Sp_DR, r_acf_tp[3, 0] * (1 - p_dst)),
            (I.Sym_Sn_DR, I.Txs_Pub_Sn_DR, r_acf_tp[2, 0] * p_dst),
            (I.Sym_Sp_DR, I.Txs_Pub_Sp_DR, r_acf_tp[3, 0] * p_dst),

            (I.ExSym_Sn_DS, I.Txf_Pub_Sn_DS, r_acf_tp[0, 0]),
            (I.ExSym_Sp_DS, I.Txf_Pub_Sp_DS, r_acf_tp[1, 0]),
            (I.ExSym_Sn_DR, I.Txf_Pub_Sn_DR, r_acf_tp[2, 0] * (1 - p_dst)),
            (I.ExSym_Sp_DR, I.Txf_Pub_Sp_DR, r_acf_tp[3, 0] * (1 - p_dst)),
            (I.ExSym_Sn_DR, I.Txs_Pub_Sn_DR, r_acf_tp[2, 0] * p_dst),
            (I.ExSym_Sp_DR, I.Txs_Pub_Sp_DR, r_acf_tp[3, 0] * p_dst),
        ]
        trs_h = [
            (I.Sym_Sn_DS, I.Txf_Pub_Sn_DS, r_acf_tp[0, 1]),
            (I.Sym_Sp_DS, I.Txf_Pub_Sp_DS, r_acf_tp[1, 1]),
            (I.Sym_Sn_DR, I.Txf_Pub_Sn_DR, r_acf_tp[2, 1] * (1 - p_dst)),
            (I.Sym_Sp_DR, I.Txf_Pub_Sp_DR, r_acf_tp[3, 1] * (1 - p_dst)),
            (I.Sym_Sn_DR, I.Txs_Pub_Sn_DR, r_acf_tp[2, 1] * p_dst),
            (I.Sym_Sp_DR, I.Txs_Pub_Sp_DR, r_acf_tp[3, 1] * p_dst),

            (I.ExSym_Sn_DS, I.Txf_Pub_Sn_DS, r_acf_tp[0, 1]),
            (I.ExSym_Sp_DS, I.Txf_Pub_Sp_DS, r_acf_tp[1, 1]),
            (I.ExSym_Sn_DR, I.Txf_Pub_Sn_DR, r_acf_tp[2, 1] * (1 - p_dst)),
            (I.ExSym_Sp_DR, I.Txf_Pub_Sp_DR, r_acf_tp[3, 1] * (1 - p_dst)),
            (I.ExSym_Sn_DR, I.Txs_Pub_Sn_DR, r_acf_tp[2, 1] * p_dst),
            (I.ExSym_Sp_DR, I.Txs_Pub_Sp_DR, r_acf_tp[3, 1] * p_dst),
        ]

        return trs_l, trs_h

    def measure(self, t, y, pars, intv, calc, mea):
        I = self.Keys

        ns = y.sum(0)
        n = ns.sum()

        det_pub = calc['det_txf_pub_s'] + calc['det_txf_pub_c'] + calc['det_txs_pub_s'] + calc['det_txs_pub_c']
        det_pri = calc['det_txf_pri_s'] + calc['det_txf_pri_c']
        det_acf = calc['acf_txf_pub_s'] + calc['acf_txf_pub_c'] + calc['acf_txs_pub_s'] + calc['acf_txs_pub_c']

        det = det_pub + det_pri + det_acf

        n_nontb = pars['NonTB'] * y.sum(0, keepdims=True) if 'NonTB' in pars else np.zeros(2)
        calc['acf_nontb'] = n_nontb * calc['r_acf']
        calc['acf_tx_fp'] = n_nontb * calc['r_acf_fp']

        if 'NonTB' in pars:
            n_nontb = pars['NonTB'] * n
            fp = n_nontb * pars['r_cs_s'] * (1 - pars['spec0'])
        else:
            fp = 0

        tp = det.sum()
        mea['TP'] = det.sum() / n
        mea['TP_Pub'] = det_pub.sum() / n
        mea['TP_Pri'] = det_pri.sum() / n
        mea['TP_Acf'] = det_acf.sum() / n
        mea['TP_DS'] = det[2:].sum() / n
        mea['TP_DR'] = det[:2].sum() / n

        fp_acf = calc['acf_tx_fp']

        mea['FP_PCF'] = fp / n
        mea['FP_ACF'] = fp_acf / n
        mea['PPV'] = tp / (fp + tp + 1e-10)
        mea['PPV_ACF'] = det_acf.sum() / (fp_acf.sum() + det_acf.sum() + 1e-10)

        mea['N_Pub_Detected'] = det_pub.sum()
        mea['N_Pri_Detected'] = det_pri.sum()
        mea['N_ACF_Detected'] = det_acf.sum()
        mea['N_ACF_TB_Reached'] = (calc['acf_s'] + calc['acf_c']).sum()
        mea['N_ACF_NonTB_Reached'] = np.array(calc['acf_nontb']).sum()
        mea['N_ACF_Reached'] = mea['N_ACF_TB_Reached'] + mea['N_ACF_NonTB_Reached']
        mea['R_ACF_Reached'] = mea['N_ACF_Reached'] / n

        for i, strata in enumerate(I.Tag_Strata):
            n = max(ns[i], 1e-15)
            mea[f'TP_{strata}'] = det[:, i].sum() / n
            mea[f'TP_DS_{strata}'] = det[2:, i].sum() / n
            mea[f'TP_DR_{strata}'] = det[:2, i].sum() / n
            mea[f'TP_ACF_{strata}'] = det_acf[:, i].sum() / n


