import numpy as np
from sim.components.base import Process
from collections import defaultdict

__all__ = ['ProcACF', 'ProcAltACF']


def _mea_acf(n, calc, prefix=''):
    arrived = calc[f'{prefix}arrived'].sum()
    eligible = calc[f'{prefix}eligible'].sum()
    pos = calc[f'{prefix}pos'].sum()

    fl = calc[f'{prefix}tp_ds'].sum() + calc[f'{prefix}tp_dr_fl'].sum()
    sl = calc['fu_tp_dr_sl'].sum()

    return {
        f'ACF_{prefix}Footfall': arrived / n,
        f'ACF_{prefix}Screened': eligible / n,
        f'ACF_{prefix}Yield': pos / n,
        f'ACF_{prefix}TP': (fl + sl) / n,
        f'ACF_{prefix}Fl': fl / n,
        f'ACF_{prefix}Sl': sl / n,
        f'ACF_{prefix}TPT': (pos - fl - sl) / n
    }


class ProcACF(Process):
    def _calc(self, t, y, pars):
        I = self.Keys

        pars, intv = pars
        n = y.sum()
        eli = pars['eli']
        p_eli = (eli * y).sum() / n

        calc = defaultdict(lambda: np.zeros((y.shape[0], 1)))

        r_acf, r_loss, r_fu, alg = intv.find_rates_main(t, p_eli)
        calc['r_acf'], calc['r_loss'] = r_acf, r_loss

        if r_acf > 0:
            alg = intv.FullACF.ScreenAlg
            alg = pars[f'alg:{alg}']
            p_dst = pars['acf_dst_sens']

            calc['arrived'] = arrived = r_acf * y
            calc['eligible'] = eligible = arrived * eli
            calc['pos'] = pos = eligible * alg['pos']
            calc['neg'] = eligible * (1 - alg['pos'])
            calc['lost'] = r_loss * y[:, 2:]

            tp_ds, tp_dr = pos[I.Infectious_DS], pos[I.Infectious_DR]
            calc['tp_ds'] = tp_ds
            calc['tp_dr'] = tp_dr
            calc['tp_dr_fl'] = tp_dr * (1 - p_dst)
            calc['tp_dr_sl'] = tp_dr * p_dst

            for item in ['sym', 'vul', 'vs', 'cxr', 'xpert']:
                calc[f'use_{item}'] = (alg[f'use_{item}'] * eligible).sum()

            # Follow-up period
            if r_loss > 0:
                alg = pars['alg:Sy']

                calc['fu_arrived'] = arrived = np.array([0, 0, r_fu, r_fu]).reshape((1, 4)) * y
                calc['fu_eligible'] = eligible = arrived * eli
                calc['fu_pos'] = eligible * alg['pos']
                calc['fu_neg'] = eligible * (1 - alg['pos'])

                tp_ds, tp_dr = pos[I.Infectious_DS], pos[I.Infectious_DR]
                calc['fu_tp_ds'] = tp_ds
                calc['fu_tp_dr'] = tp_dr
                calc['fu_tp_dr_fl'] = tp_dr * (1 - p_dst)
                calc['fu_tp_dr_sl'] = tp_dr * p_dst

                for item in ['sym', 'vul', 'vs', 'cxr', 'xpert']:
                    calc[f'use_{item}'] += (alg[f'use_{item}'] * eligible).sum()
        return calc

    def calc_dy(self, t, y, pars):
        I = self.Keys

        calc = self._calc(t, y, pars)

        tp_ds, tp_dr = calc['tp_ds'], calc['tp_dr']
        tp_dr_fl, tp_dr_sl = calc['tp_dr_fl'], calc['tp_dr_sl']

        dy = np.zeros_like(y)

        if calc['r_acf'] <= 0:
            return dy

        dy[I.Infectious_DS] -= tp_ds
        dy[I.Infectious_DR] -= tp_dr
        dy[I.Txf_Pub_DS] += tp_ds.sum(0)
        dy[I.Txf_Pub_DR] += tp_dr_fl.sum(0)
        dy[I.Txs_Pub_DR] += tp_dr_sl.sum(0)

        pos = calc['pos']
        fp = pos[I.LTBI0]
        dy[I.LTBI0] -= fp
        dy[I.LTBI_TPT] += fp
        fp = pos[I.U]
        dy[I.U] -= fp
        dy[I.UTPT] += fp

        complete = y[I.LTBI_TPT] * 2
        dy[I.LTBI_TPT] -= complete
        dy[I.LTBI0, :2] += complete[:, :2] + complete[:, 2:]

        complete = y[I.UTPT] * 2
        dy[I.UTPT] -= complete
        dy[I.U, :2] += complete[:2] + complete[2:]

        # to follow up list
        if calc['r_loss'] <= 0:
            return dy

        tp_ds, tp_dr = calc['fu_tp_ds'], calc['fu_tp_dr']
        tp_dr_fl, tp_dr_sl = calc['fu_tp_dr_fl'], calc['fu_tp_dr_sl']

        dy[I.Infectious_DS] -= tp_ds
        dy[I.Infectious_DR] -= tp_dr
        dy[I.Txf_Pub_DS] += tp_ds.sum(0)
        dy[I.Txf_Pub_DR] += tp_dr_fl.sum(0)
        dy[I.Txs_Pub_DR] += tp_dr_sl.sum(0)

        pos = calc['fu_pos']
        fp = pos[I.LTBI0]
        dy[I.LTBI0] -= fp
        dy[I.LTBI_TPT] += fp
        fp = pos[I.U]
        dy[I.U] -= fp
        dy[I.UTPT] += fp

        neg, lost = calc['neg'], calc['lost']
        dy[:, 1] -= neg[:, 1]
        dy[:, 3] += neg[:, 1]

        dy[:, 2:] -= lost
        dy[:, :2] += lost
        return dy

    def measure(self, t, y, pars, mea):
        n = y.sum()
        calc = self._calc(t, y, pars)
        mea.update(_mea_acf(n, calc))
        mea.update(_mea_acf(n, calc, prefix='fu_'))

        for item in ['sym', 'vul', 'vs', 'cxr', 'xpert']:
            mea[f'ACF_N_{item}'] = calc[f'use_{item}'].sum()

        n_tpt = y[self.Keys.UTPT].sum() + y[self.Keys.LTBI_TPT].sum()
        mea['PrOnTPT'] = n_tpt / n


class ProcAltACF(Process):
    def _calc(self, t, y, pars):
        I = self.Keys

        pars, intv = pars
        n = y.sum()
        eli = pars['eli']
        p_eli = (eli * y).sum() / n

        calc = defaultdict(lambda: np.zeros((y.shape[0], 1)))

        r_acf, alg = intv.find_rates_alt(t, p_eli)
        calc['r_acf'] = r_acf

        if r_acf > 0:
            alg = intv.AltACF.ScreenAlg
            alg = pars[f'alg:{alg}']
            p_dst = pars['acf_dst_sens']

            calc['arrived'] = arrived = r_acf * y
            calc['eligible'] = eligible = arrived * eli
            calc['pos'] = pos = eligible * alg['pos']

            tp_ds, tp_dr = pos[I.Infectious_DS], pos[I.Infectious_DR]
            calc['tp_ds'] = tp_ds
            calc['tp_dr'] = tp_dr
            calc['tp_dr_fl'] = tp_dr * (1 - p_dst)
            calc['tp_dr_sl'] = tp_dr * p_dst

            for item in ['sym', 'vul', 'vs', 'cxr', 'xpert']:
                calc[f'use_{item}'] = (alg[f'use_{item}'] * eligible).sum()

        return calc

    def calc_dy(self, t, y, pars):
        I = self.Keys

        calc = self._calc(t, y, pars)
        tp_ds, tp_dr = calc['tp_ds'], calc['tp_dr']
        tp_dr_fl, tp_dr_sl = calc['tp_dr_fl'], calc['tp_dr_sl']

        dy = np.zeros_like(y)

        if calc['r_acf'] <= 0:
            return dy

        dy[I.Infectious_DS] -= tp_ds
        dy[I.Infectious_DR] -= tp_dr
        dy[I.Txf_Pub_DS] += tp_ds.sum(0)
        dy[I.Txf_Pub_DR] += tp_dr_fl.sum(0)
        dy[I.Txs_Pub_DR] += tp_dr_sl.sum(0)

        pos = calc['pos']
        fp = pos[I.LTBI0]
        dy[I.LTBI0] -= fp
        dy[I.LTBI_TPT] += fp
        fp = pos[I.U]
        dy[I.U] -= fp
        dy[I.UTPT] += fp

        complete = y[I.LTBI_TPT] * 2
        dy[I.LTBI_TPT] -= complete
        dy[I.LTBI0, :2] += complete[:, :2] + complete[:, 2:]

        complete = y[I.UTPT] * 2
        dy[I.UTPT] -= complete
        dy[I.U, :2] += complete[:2] + complete[2:]

        return dy

    def measure(self, t, y, pars, mea):
        n = y.sum()
        calc = self._calc(t, y, pars)
        mea.update(_mea_acf(n, calc, 'Alt'))

        for item in ['sym', 'vul', 'vs', 'cxr', 'xpert']:
            mea[f'AltACF_N_{item}'] = calc[f'use_{item}'].sum()
