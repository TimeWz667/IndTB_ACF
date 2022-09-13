import numpy as np
from sim.components.base import Process

__author__ = 'Chu-Chang Ku'
__all__ = ['Transmission']


def infection(sus, trans, y):
    foi = (trans * y).sum() / y.sum()
    return (sus * y) * foi


class Transmission(Process):
    def __call__(self, t, y, pars, calc):
        calc['infection_ds'] = infection(sus=pars['sus'], trans=pars['trans_ds'], y=y) * pars['beta_ds']
        calc['infection_dr'] = infection(sus=pars['sus'], trans=pars['trans_dr'], y=y) * pars['beta_dr']

    def calc_dy(self, t, y, pars):
        I = self.Keys
        calc = dict()
        self(t, y, pars, calc)

        dy = np.zeros_like(y)

        # Infection
        dy -= calc['infection_ds'] + calc['infection_dr']
        dy[I.FLat_DS] += calc['infection_ds'].sum(0)
        dy[I.FLat_DR] += calc['infection_dr'].sum(0)

        return dy

    def measure(self, t, y, pars, mea):
        I = self.Keys
        calc = dict()
        self(t, y, pars, calc)

        inf_ds = calc['infection_ds'].sum(0)
        inf_dr = calc['infection_dr'].sum(0)
        inf = inf_ds + inf_dr

        ltbi = y[I.LTBI].sum(0)

        ns = y.sum(0)
        n = ns.sum()

        mea['ARTI'] = inf.sum() / n
        mea['LTBI'] = ltbi.sum() / n

        for i, strata in enumerate(I.Tag_Strata):
            n = max(ns[i], 1e-15)

            mea[f'ARTI_{strata}'] = inf[i] / n
            mea[f'LTBI_{strata}'] = ltbi[i] / n

        if mea['ARTI_RiskLo'] <= 0:
            mea['RR_inf_comorb'] = 0
        else:
            mea['RR_inf_comorb'] = mea['ARTI_RiskHi'] / mea['ARTI_RiskLo']
