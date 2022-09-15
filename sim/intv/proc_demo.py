import numpy as np
from sim.components.base import Process

__author__ = 'Chu-Chang Ku'
__all__ = ['Demography']


class Demography(Process):
    def __call__(self, t, y, pars, calc):
        I = self.Keys

        dr = np.ones_like(y) * pars['r_die']
        gr = pars['r_growth'] if t > 1970 else 0

        dr_tb = np.zeros_like(y)
        dr_tb[I.Asym] = pars['r_die_ut']
        dr_tb[I.Sym] = pars['r_die_ut']
        dr_tb[I.ExSym] = pars['r_die_ut']
        dr_tb[I.Tx_DS] = pars['r_die_ds']
        dr_tb[I.Tx_DR] = pars['r_die_dr']

        n = y.sum(0)
        calc['n'] = n
        calc['deaths'] = dr * y
        calc['deaths_tb'] = dr_tb * y
        calc['deaths_all'] = death_all = calc['deaths'] + calc['deaths_tb']
        calc['dr'] = dr
        calc['dr_tb'] = dr_tb
        calc['gr'] = gr

        calc['births'] = (calc['deaths'] + calc['deaths_tb']).sum() + gr * n.sum()

        r_die_crude = death_all[:, [1, 3]].sum() / y[:, [1, 3]].sum()

        if t > 1970:
            r_comorb = (pars['p_comorb'] / (1 - pars['p_comorb'])) * (gr + r_die_crude)
        else:
            r_comorb = (pars['p_comorb'] / (1 - pars['p_comorb'])) * r_die_crude
        calc['prog_comorb'] = r_comorb * y[:, [0, 2]]

    def calc_dy(self, t, y, pars):
        I = self.Keys

        calc = dict()
        self(t, y, pars, calc)
        dy = np.zeros_like(y)
        dy[I.U, 0] += calc['births']
        dy -= calc['deaths_all']

        dy[:, [0, 2]] -= calc['prog_comorb']
        dy[:, [1, 3]] += calc['prog_comorb']
        return dy

    def measure(self, t, y, pars, mea):
        I = self.Keys
        calc = dict()
        self(t, y, pars, calc)

        ns = y.sum(0)
        mor = calc['deaths'][I.Infectious].sum(0)
        mor_tb = calc['deaths_tb'][I.Infectious].sum(0)

        n = ns.sum()
        mea['Pop'] = n
        mea['MorR'] = (mor + mor_tb).sum() / n
        mea['PropComorb'] = ns[[1, 3]].sum() / n

        for i, strata in enumerate(I.Tag_Strata):
            n = max(ns[[i, i + 2]].sum(), 1e-15)
            mea[f'Pop_{strata}'] = n
            mea[f'MorR_{strata}'] = (mor + mor_tb)[[i, i + 2]].sum() / n