from sim.components.base import Process
from sim.util import calc_dy, extract_tr
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['Progression']


class Progression(Process):
    def __init__(self, keys):
        Process.__init__(self, keys)

    def get_trs(self, t, y, pars):
        I = self.Keys

        pr0 = pars['p_primary']
        r_act = pars['r_stab'] * pr0 / (1 - pr0)
        r_react = pars['r_react']
        r_rel_st, r_rel_tc, r_rel_td = pars['r_relapse_st'], pars['r_relapse_tc'], pars['r_relapse_td']

        rr_comorb = pars['rr_risk_comorb']

        inc_l = [
            (I.FLat_DS, I.Asym_DS, r_act, 'inc_recent'),
            (I.FLat_DR, I.Asym_DR, r_act, 'inc_recent'),

            (I.SLat_DS, I.Asym_DS, r_react, 'inc_remote'),
            (I.SLat_DR, I.Asym_DR, r_react, 'inc_remote'),

            (I.RLow_DS, I.Asym_DS, r_rel_tc, 'inc_remote'),
            (I.RLow_DR, I.Asym_DR, r_rel_tc, 'inc_remote'),

            (I.RHigh_DS, I.Asym_DS, r_rel_td, 'inc_remote'),
            (I.RHigh_DR, I.Asym_DR, r_rel_td, 'inc_remote'),

            (I.RSt_DS, I.Asym_DS, r_rel_st, 'inc_remote'),
            (I.RSt_DR, I.Asym_DR, r_rel_st, 'inc_remote'),

            (I.FLatTPT_DS, I.Asym_DS, 0.3 * r_act, 'inc_recent'),
            (I.FLatTPT_DR, I.Asym_DR, 0.3 * r_act, 'inc_recent'),

            (I.SLatTPT_DS, I.Asym_DS, 0.3 * r_react, 'inc_remote'),
            (I.SLatTPT_DR, I.Asym_DR, 0.3 * r_react, 'inc_remote'),

            (I.RLowTPT_DS, I.Asym_DS, 0.3 * r_rel_tc, 'inc_remote'),
            (I.RLowTPT_DR, I.Asym_DR, 0.3 * r_rel_tc, 'inc_remote'),

            (I.RHighTPT_DS, I.Asym_DS, 0.3 * r_rel_td, 'inc_remote'),
            (I.RHighTPT_DR, I.Asym_DR, 0.3 * r_rel_td, 'inc_remote'),

            (I.RStTPT_DS, I.Asym_DS, 0.3 * r_rel_st, 'inc_remote'),
            (I.RStTPT_DR, I.Asym_DR, 0.3 * r_rel_st, 'inc_remote'),
        ]

        inc_h = list(inc_l)
        for i in range(8):
            fr, to, rate, desc = inc_h[i]
            inc_h[i] = fr, to, rate * rr_comorb, desc

        prog = list()

        r_stab = pars['r_stab']
        for fr, to in zip(I.FLat, I.SLat):
            prog.append((fr, to, r_stab, 'stab'))
        for fr, to in zip(I.RLow, I.RSt):
            prog.append((fr, to, r_stab, 'stab'))
        for fr, to in zip(I.RHigh, I.RSt):
            prog.append((fr, to, r_stab, 'stab'))

        for fr, to in zip(I.FLatTPT, I.SLatTPT):
            prog.append((fr, to, r_stab, 'stab'))
        for fr, to in zip(I.RLowTPT, I.RStTPT):
            prog.append((fr, to, r_stab, 'stab'))
        for fr, to in zip(I.RHighTPT, I.RStTPT):
            prog.append((fr, to, r_stab, 'stab'))

        r_sc = pars['r_sc']
        for fr in I.Infectious_DS:
            prog.append((fr, I.RHigh_DS, r_sc, 'self_cure'))

        for fr in I.Infectious_DR:
            prog.append((fr, I.RHigh_DR, r_sc, 'self_cure'))

        r_clear = pars['r_clear']
        for fr in I.SLat + I.RSt + I.SLatTPT + I.RStTPT:
            prog.append((fr, I.U, r_clear, 'self_clear'))

        for fr, to in zip(I.Asym, I.Sym):
            prog.append((fr, to, pars['r_onset'], 'onset'))

        if t > 1980:
            r_mdr_tx = pars['r_mdr_tx']
            prog += [
                (I.Txf_Pub_DS, I.Txf_Pub_DR, r_mdr_tx, 'dev_dr'),
                (I.Txf_Pri_DS, I.Txf_Pri_DR, r_mdr_tx, 'dev_dr')
            ]

        trs_l = inc_l + prog
        trs_h = inc_h + prog
        return trs_l, trs_h

    def calc_dy(self, t, y, pars):
        trs_l, trs_h = self.get_trs(t, y, pars)
        dy = np.zeros_like(y)
        dy[:, 0] = calc_dy(y[:, 0], trs_l)
        dy[:, 1] = calc_dy(y[:, 1], trs_h)
        return dy

    def measure(self, t, y, pars, mea):
        I = self.Keys
        trs_l, trs_h = self.get_trs(t, y, pars)

        inc_recent = np.array([
            extract_tr(y[:, 0], trs_l, lambda x: x[3] == 'inc_recent'),
            extract_tr(y[:, 1], trs_h, lambda x: x[3] == 'inc_recent')
        ])
        inc_remote = np.array([
            extract_tr(y[:, 0], trs_l, lambda x: x[3] == 'inc_remote'),
            extract_tr(y[:, 1], trs_h, lambda x: x[3] == 'inc_remote')
        ])

        inc_dr = np.array([
            extract_tr(y[:, 0], trs_l, lambda x: x[3].startswith('inc') and x[1] is I.Asym_DR),
            extract_tr(y[:, 1], trs_h, lambda x: x[3].startswith('inc') and x[1] is I.Asym_DR)
        ])

        ns = y.sum(0)
        n = ns.sum()
        inc = inc_recent + inc_remote

        mea['IncR'] = inc.sum() / n
        mea['IncR_Recent'] = inc_recent.sum() / n
        mea['IncR_Remote'] = inc_remote.sum() / n
        mea['Recent'] = inc_recent.sum() / max(inc.sum(), 1e-10)
        mea['IncR_DR'] = inc_dr.sum() / n
        mea['IncR_DS'] = mea['IncR'] - mea['IncR_DR']
        mea['PrDR_Inc'] = mea['IncR_DR'] / max(mea['IncR'], 1e-10)

        prev_a = y[I.Asym].sum(0)
        prev_s = y[I.Sym].sum(0)
        prev_c = y[I.ExSym].sum(0)

        prev = prev_a + prev_s + prev_c
        mea['Prev'] = prev.sum() / n
        mea['N_UT'] = prev.sum()
        mea['PrAsym'] = prev_a.sum() / prev.sum()
        mea['PrPreCS'] = prev_s.sum() / prev.sum()
        mea['PrExCS'] = prev_c.sum() / prev.sum()

        # mea['PrSp_Asym'] = y[I.Asym][[1, 3]].sum() / prev_a.sum()
        # mea['PrSp_PreCS'] = y[I.Sym][[1, 3]].sum() / max(prev_s.sum(), 1e-10)
        # mea['PrSp_ExCS'] = y[I.ExSym][[1, 3]].sum() / max(prev_c.sum(), 1e-10)

        ltbi = y[I.LTBI].sum(0)
        mea['LTBI'] = ltbi.sum() / n

        for i, strata in enumerate(I.Tag_Strata):
            n = max(ns[i], 1e-15)
            mea[f'IncR_{strata}'] = inc.sum() / n
            mea[f'IncR_DS_{strata}'] = inc[i].sum() / n
            mea[f'IncR_DR_{strata}'] = inc[i].sum() / n
            mea[f'Recent_{strata}'] = inc_recent[i].sum() / max(inc[i].sum(), 1e-15)

            mea[f'Prev_{strata}'] = prev[i] / n
            mea[f'LTBI_{strata}'] = ltbi[i] / n

        mea['RR_inc_comorb'] = mea['IncR_RiskHi'] / max(mea['IncR_RiskLo'], 1e-10)
        p1, p0 = mea['Prev_RiskHi'], mea['Prev_RiskLo']
        if p1 <= 0 or p0 <= 0:
            mea['OR_prev_comorb'] = 0
        else:
            mea['OR_prev_comorb'] = (p1 / (1 - p1)) / (p0 / (1 - p0))
