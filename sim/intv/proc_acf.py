import numpy as np


__all__ = ['BgACF', 'VulACF']


def _mea_acf(n, arrived, screened, confirmed, pos, tp_ds, tp_dr_fl, tp_dr_sl, label):
    return {
        f'ACF_{label}_Footfall': arrived.sum() / n,
        f'ACF_{label}_Screened': screened.sum() / n,
        f'ACF_{label}_Confirmed': confirmed.sum() / n,
        f'ACF_{label}_Yield': pos.sum() / n,
        f'ACF_{label}_TP': (tp_ds + tp_dr_fl + tp_dr_sl).sum() / n,
        f'ACF_{label}_DS_Fl': tp_ds.sum() / n,
        f'ACF_{label}_DR_Fl': tp_dr_fl.sum() / n,
        f'ACF_{label}_DR_Sl': tp_dr_sl.sum() / n,
    }


class BgACF:
    def __init__(self, keys):
        self.Keys = keys

    def _calc_bg_acf(self, y, r_acf, eli, pos_screen, pos_confirm, p_dst, mea=False, label='bg'):
        I = self.Keys
        arrived = r_acf * y
        screened = arrived * eli
        confirmed = screened * pos_screen
        pos = confirmed * pos_confirm
        tp_ds, tp_dr = pos[I.Infectious_DS], pos[I.Infectious_DR]
        tp_dr_fl = tp_dr * (1 - p_dst)
        tp_dr_sl = tp_dr - tp_dr_fl

        if mea:
            n = y.sum()
            return _mea_acf(n, arrived, screened, confirmed, pos, tp_ds, tp_dr_fl, tp_dr_sl, label)
        else:
            dy = np.zeros_like(y)
            dy[I.Infectious_DS, :2] -= tp_ds
            dy[I.Infectious_DR, :2] -= tp_dr
            dy[I.Txf_Pub_DS, :2] += tp_ds.sum(0)
            dy[I.Txf_Pub_DR, :2] += tp_dr_fl.sum(0)
            dy[I.Txs_Pub_DR, :2] += tp_dr_sl.sum(0)
            return dy

    def calc_dy(self, t, y, pars, intv):
        pos_sym, pos_cxr, pos_xpert, eligible = pars['pos_sym'], pars['pos_cxr'], pars['pos_xpert'], pars['eli']
        r_acf_mdu, r_acf_d2d, p_dst = pars['r_acf_mdu'], pars['r_acf_d2d'], pars['acf_dst_sens']
        r_acf_mdu, r_acf_d2d, p_dst = intv.modify_acf_bg(t, r_acf_mdu, r_acf_d2d, p_dst)

        dy = np.zeros_like(y)
        # MDU
        if r_acf_mdu > 0:
            dy[:, :2] += self._calc_bg_acf(y[:, :2], r_acf_mdu, eligible, pos_cxr, pos_xpert, p_dst, mea=False)
        # D2D
        if r_acf_d2d > 0:
            dy[:, :2] += self._calc_bg_acf(y[:, :2], r_acf_d2d, eligible, pos_sym, pos_xpert, p_dst, mea=False)

        return dy

    def measure(self, t, y, pars, intv, mea):
        pos_sym, pos_cxr, pos_xpert, eligible = pars['pos_sym'], pars['pos_cxr'], pars['pos_xpert'], pars['eli']
        r_acf_mdu, r_acf_d2d, p_dst = pars['r_acf_mdu'], pars['r_acf_d2d'], pars['acf_dst_sens']
        r_acf_mdu, r_acf_d2d, p_dst = intv.modify_acf_bg(t, r_acf_mdu, r_acf_d2d, p_dst)

        # MDU
        mea.update(self._calc_bg_acf(y[:, :2], r_acf_mdu, eligible, pos_cxr, pos_xpert, p_dst, mea=True, label='MDU'))

        # D2D
        mea.update(self._calc_bg_acf(y[:, :2], r_acf_d2d, eligible, pos_sym, pos_xpert, p_dst, mea=True, label='D2D'))

        l, lt = y[self.Keys.LTBI0].sum(), y[self.Keys.LTBI_TPT].sum()
        mea['PrOnPseudoTPT'] = lt / (l + lt)


class VulACF:
    def __init__(self, keys):
        self.Keys = keys
        self.Recap = True

    def _calc_vul_acf(self, y, cov, r_fu, r_lost, eli, pos_screen, pos_confirm, p_dst, mea=False, label='vul'):
        I = self.Keys
        n = y.sum()

        n_target = cov * y.sum()
        eli = np.hstack([eli, eli])
        r_acf = n_target / (eli * y).sum()
        r_acf = min(24, r_acf)
        arrived = r_acf * y
        screened = arrived * eli
        confirmed = screened * pos_screen
        pos = confirmed * pos_confirm
        neg = screened - pos
        tp_ds, tp_dr = pos[I.Infectious_DS], pos[I.Infectious_DR]
        tp_dr_fl = tp_dr * (1 - p_dst)
        tp_dr_sl = tp_dr - tp_dr_fl

        if mea:
            m = _mea_acf(n, arrived, screened, confirmed, pos, tp_ds, tp_dr_fl, tp_dr_sl, label)
        else:
            dy = np.zeros_like(y)
            dy[I.Infectious_DS] -= tp_ds
            dy[I.Infectious_DR] -= tp_dr
            dy[I.Txf_Pub_DS] += tp_ds.sum(0)
            dy[I.Txf_Pub_DR] += tp_dr_fl.sum(0)
            dy[I.Txs_Pub_DR] += tp_dr_sl.sum(0)

            dy[:, :2] -= neg[:, :2]
            dy[:, 2:] += neg[:, :2]

            fp = pos[I.LTBI0]
            dy[I.LTBI0] -= fp
            dy[I.LTBI_TPT] += fp

            complete = y[I.LTBI_TPT] * 2

            dy[I.LTBI_TPT] -= complete
            dy[I.LTBI0, :2] += complete[:, :2] + complete[:, 2:]

        if r_lost <= 0:
            if mea:
                return m
            else:
                return dy

        # Follow-up
        arrived = r_fu * y[:, 2:]
        screened = arrived
        confirmed = screened * pos_screen
        pos = confirmed * pos_confirm
        tp_ds, tp_dr = pos[I.Infectious_DS], pos[I.Infectious_DR]
        tp_dr_fl = tp_dr * (1 - p_dst)
        tp_dr_sl = tp_dr - tp_dr_fl

        if mea:
            m.update(_mea_acf(n, arrived, screened, confirmed, pos, tp_ds, tp_dr_fl, tp_dr_sl, f'{label}fu'))
        else:
            dy[I.Infectious_DS, 2:] -= tp_ds
            dy[I.Infectious_DR, 2:] -= tp_dr
            dy[I.Txf_Pub_DS, 2:] += tp_ds.sum(0)
            dy[I.Txf_Pub_DR, 2:] += tp_dr_fl.sum(0)
            dy[I.Txs_Pub_DR, 2:] += tp_dr_sl.sum(0)

            lost = r_lost * y[:, 2:]
            dy[:, 2:] -= lost
            dy[:, :2] += lost

        if mea:
            return m
        else:
            return dy

    def calc_dy(self, t, y, pars, intv):
        pos_cxr, pos_xpert, p_dst = pars['pos_cxr'], pars['pos_xpert'], pars['acf_dst_sens']
        eligible, eli_vul = pars['eli'], pars['eli_vul']

        dy = np.zeros_like(y)

        # Vulnerability-led ACF
        r_acf_vul, r_fu_v, r_lost_v = intv.modify_acf_vul(t, 0, 0, 0)
        if r_acf_vul > 0:
            dy += self._calc_vul_acf(y, r_acf_vul, r_fu_v, r_lost_v, eli_vul, pos_cxr, pos_xpert, p_dst, mea=False)

        # Plain ACF
        r_acf_plain, r_fu_p, r_lost_p = intv.modify_acf_plain(t, 0, 0, 0)
        if r_acf_plain > 0:
            dy += self._calc_vul_acf(y, r_acf_plain, r_fu_p, r_lost_p, eligible, pos_cxr, pos_xpert, p_dst, mea=False)
        return dy

    def measure(self, t, y, pars, intv, mea):
        pos_cxr, pos_xpert, p_dst = pars['pos_cxr'], pars['pos_xpert'], pars['acf_dst_sens']
        eligible, eli_vul = pars['eli'], pars['eli_vul']

        r_acf_vul, r_fu, r_lost = intv.modify_acf_vul(t, 0, 0, 0)
        mea.update(self._calc_vul_acf(y, r_acf_vul, r_fu, r_lost, eli_vul, pos_cxr, pos_xpert, p_dst,
                                      mea=True, label='Vul'))

        # Plain ACF
        r_acf_plain, r_fu, r_lost = intv.modify_acf_plain(t, 0, 0, 0)
        mea.update(self._calc_vul_acf(y, r_acf_plain, r_fu, r_lost, eligible, pos_cxr, pos_xpert, p_dst,
                                      mea=True, label='Plain'))
