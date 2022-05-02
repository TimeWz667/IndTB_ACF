from sim.components.base import Process
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['Progression']


class Progression(Process):
    def __init__(self, keys, intv=None):
        Process.__init__(self, keys)
        self.Intervention = intv

    def __call__(self, t, y, pars, calc):
        I = self.Keys

        r_react = pars['r_react']
        r_rel_st, r_rel_tc, r_rel_td = pars['r_relapse_st'], pars['r_relapse_tc'], pars['r_relapse_td']

        p_sp0 = pars['p_sp0']
        sn_sp = np.array([1 - p_sp0, p_sp0, 1 - p_sp0, p_sp0]).reshape((-1, 1))

        risk_comorb = np.ones(y[I.SLat].shape)
        risk_comorb[:, I.RiskHi] *= pars['rr_risk_comorb']
        calc['react'] = risk_comorb * r_react * y[I.SLat]
        calc['rel_tc'] = risk_comorb * r_rel_tc * y[I.RLow]
        calc['rel_td'] = risk_comorb * r_rel_td * y[I.RHigh]
        calc['rel_st'] = risk_comorb * r_rel_st * y[I.RSt]

        calc['act_smr'] = sn_sp * np.repeat(calc['act'], 2, axis=0)
        calc['react_smr'] = sn_sp * np.repeat(calc['react'], 2, axis=0)
        calc['rel_tc_smr'] = sn_sp * np.repeat(calc['rel_tc'], 2, axis=0)
        calc['rel_td_smr'] = sn_sp * np.repeat(calc['rel_td'], 2, axis=0)
        calc['rel_st_smr'] = sn_sp * np.repeat(calc['rel_st'], 2, axis=0)

        calc['inc_recent'] = calc['act']
        calc['inc_remote'] = calc['react'] + calc['rel_tc'] + calc['rel_td'] + calc['rel_st']
        calc['inc'] = calc['inc_recent'] + calc['inc_remote']
        calc['inc_smr'] = calc['act_smr'] + calc['react_smr'] + calc['rel_tc_smr'] + calc['rel_td_smr'] + calc['rel_st_smr']

        r_stab = pars['r_stab']
        calc['stab_tc'] = r_stab * y[I.RLow]
        calc['stab_td'] = r_stab * y[I.RHigh]

        calc['sc_a'] = pars['r_sc'] * y[I.Asym]
        calc['sc_s'] = pars['r_sc'] * y[I.Sym]
        calc['sc_c'] = pars['r_sc'] * y[I.ExSym]

        # calc['clear_sl'] = pars['r_clear'] * y[I.SLat]
        # calc['clear_rst'] = pars['r_clear'] * y[I.RSt]

        r_onset = np.array([
            pars['r_onset_sn'], pars['r_onset_sp'],
            pars['r_onset_sn'], pars['r_onset_sp']
        ]).reshape((-1, 1))

        calc['sym_onset'] = r_onset * y[I.Asym]

        r_convert_a = pars['r_convert_a']
        r_convert_s = pars['r_convert_s']

        calc['convert_a'] = r_convert_a * y[I.Asym_Sn]
        calc['convert_s'] = r_convert_s * y[I.Sym_Sn]
        calc['convert_c'] = r_convert_s * y[I.ExSym_Sn]

        r_mdr_tx = pars['r_mdr_tx']
        calc['develop_dr_pub'] = r_mdr_tx * y[[I.Txf_Pub_Sn_DS, I.Txf_Pub_Sp_DS]]
        calc['develop_dr_pri'] = r_mdr_tx * y[[I.Txf_Pri_Sn_DS, I.Txf_Pri_Sp_DS]]

    def measure(self, t, y, pars, calc, mea):
        I = self.Keys
        ns = y.sum(0)
        n = ns.sum()
        inc = calc['inc']

        mea['IncR'] = inc.sum() / n
        mea['Recent'] = calc['inc_recent'].sum() / inc.sum()
        mea['IncR_DR'] = inc[1].sum() / n
        mea['IncR_DS'] = inc[0].sum() / n
        mea['PrDR_Inc'] = mea['IncR_DR'] / mea['IncR']

        prev_a = y[I.Asym].sum(0)
        prev_s = y[I.Sym].sum(0)
        prev_c = y[I.ExSym].sum(0)
        prev_sym = prev_s + prev_c
        prev = prev_a + prev_s + prev_c
        mea['Prev'] = prev.sum() / n
        mea['PrSym'] = prev_sym.sum() / prev.sum()

        mea['PrSp_Asym'] = y[I.Asym][[1, 3]].sum() / prev_a.sum()
        mea['PrSp_Sym'] = (y[I.Sym] + y[I.ExSym])[[1, 3]].sum() / prev_sym.sum()

        ltbi = y[I.LTBI].sum(0)
        mea['LTBI'] = ltbi.sum() / n

        for i, strata in enumerate(I.Tag_Strata):
            n = max(ns[i], 1e-15)
            mea[f'IncR_{strata}'] = inc.sum() / n
            mea[f'IncR_DS_{strata}'] = inc[i].sum() / n
            mea[f'IncR_DR_{strata}'] = inc[i].sum() / n
            mea[f'Recent_{strata}'] = calc['inc_recent'][i].sum() / max(inc[i].sum(), 1e-15)

            mea[f'Prev_{strata}'] = prev[i] / n
            mea[f'LTBI_{strata}'] = ltbi[i] / n

        mea['RR_inc_comorb'] = mea['IncR_RiskHi'] / mea['IncR_RiskLo']
        p1, p0 = mea['Prev_RiskHi'], mea['Prev_RiskLo']
        mea['OR_prev_comorb'] = (p1 / (1 - p1)) / (p0 / (1 - p0))
