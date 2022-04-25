import numpy as np
from numba import njit
from sim.components.base import Process

__author__ = 'Chu-Chang Ku'
__all__ = ['Transmission']


@njit
def infection(sus, trans, mixing, y):
    n_tb, n_area = y.shape
    infectious = np.zeros(n_area)

    for t_f in range(n_tb):
        for a_f in range(n_area):
            if trans[t_f, a_f] > 0:
                infectious[a_f] += trans[t_f, a_f] * y[t_f, a_f]

    n = y.sum()

    infection = np.zeros_like(sus)
    for t_t in range(n_tb):
        for a_t in range(n_area):
            for a_f in range(n_area):
                infection[t_t, a_t] += infectious[a_f] * mixing[a_f, a_t] * sus[t_t, a_t] * y[t_t, a_t]

    return infection / n


class Transmission(Process):
    def __call__(self, t, y, pars, calc):
        y0_decline, y0_baseline = pars['y0_decline'], pars['y0_baseline']

        adr = pars['adr'] + pars['adr_adj']

        if t < y0_decline:
            t = y0_decline

        adj = np.exp(-adr * (t - y0_baseline))

        calc['infection'] = infection(
            sus=pars['sus'],
            trans=pars['trans'],
            y=y,
            mixing=pars['mixing']
        ) * adj * pars['beta']

    def measure(self, t, y, pars, calc, mea):
        I = self.Keys

        inf = calc['infection'].sum(0)

        ltbi = y[I.LTBI].sum(0)
        prev = y[[I.Asym, I.Sym, I.ExSym]].sum(0)

        ns = y.sum(0)
        n = ns.sum()

        mea['Prev'] = prev.sum() / n
        mea['ARTI'] = inf.sum() / n
        mea['LTBI'] = ltbi.sum() / n

        for i, strata in enumerate(I.Tag_Strata):
            n = max(ns[i], 1e-15)

            mea[f'Prev_{strata}'] = prev[i] / n
            mea[f'ARTI_{strata}'] = inf[i] / n
            mea[f'LTBI_{strata}'] = ltbi[i] / n
