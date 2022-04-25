from sim.components.base import Process

__author__ = 'Chu-Chang Ku'
__all__ = ['Progression']


class Progression(Process):
    def __init__(self, keys, intv=None):
        Process.__init__(self, keys)
        self.Intervention = intv

    def __call__(self, t, y, pars, calc):
        I = self.Keys

        r_act, r_react, r_rel = pars['r_act'], pars['r_react'], pars['r_relapse']
        r_lat = pars['r_lat']

        r_sym = pars['r_onset']

        calc['act'] = r_act * y[I.FLat]
        calc['react'] = r_react * y[I.SLat]
        calc['rel_tc'] = pars['r_relapse_tc'] * y[I.RLow]
        calc['rel_td'] = pars['r_relapse_td'] * y[I.RHigh]
        calc['rel_st'] = r_rel * y[I.RSt]

        calc['inc_recent'] = calc['act']
        calc['inc_remote'] = calc['react'] + calc['rel_tc'] + calc['rel_td'] + calc['rel_st']
        calc['inc'] = calc['inc_recent'] + calc['inc_remote']

        calc['stab_fl'] = r_lat * y[I.FLat]
        calc['stab_tc'] = pars['r_stab'] * y[I.RLow]
        calc['stab_td'] = pars['r_stab'] * y[I.RHigh]

        calc['sc_a'] = pars['r_sc'] * y[I.Asym]
        calc['sc_s'] = pars['r_sc'] * y[I.Sym]
        calc['sc_c'] = pars['r_sc'] * y[I.ExSym]

        calc['clear_sl'] = pars['r_clear'] * y[I.SLat]
        calc['clear_rst'] = pars['r_clear'] * y[I.RSt]

        calc['sym_onset'] = r_sym * y[I.Asym]

    def measure(self, t, y, pars, calc, mea):
        I = self.Keys
        ns = y.sum(0)
        n = ns.sum()
        inc = calc['inc']

        mea['IncR'] = inc.sum() / n
        mea['Recent'] = calc['inc_recent'].sum() / inc.sum()

        for i, strata in enumerate(I.Tag_Strata):
            n = max(ns[i], 1e-15)
            mea[f'IncR_{strata}'] = inc.sum() / n
            mea[f'IncR_DS_{strata}'] = inc[i].sum() / n
            mea[f'IncR_DR_{strata}'] = inc[i].sum() / n
            mea[f'Recent_{strata}'] = calc['inc_recent'][i].sum() / max(inc[i].sum(), 1e-15)
