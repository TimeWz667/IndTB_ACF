import numpy as np


__all__ = ['BgACF', 'VulACF']


def _mea_acf(n, arrived, screened, confirmed, pos, tp_ds, tp_dr_fl, tp_dr_sl, label, n_sym=0, n_vul=0, n_cxr=0, n_xpert=0):
    return {
        f'ACF_{label}_Footfall': arrived.sum() / n,
        f'ACF_{label}_Screened': screened.sum() / n,
        f'ACF_{label}_Confirmed': confirmed.sum() / n,
        f'ACF_{label}_Sym': n_sym / n,
        f'ACF_{label}_Vul': n_vul / n,
        f'ACF_{label}_CXR': n_cxr / n,
        f'ACF_{label}_Xpert': n_xpert / n,
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

    def _calc_mdu(self, y, r_acf, p_dst, pars, mea=False):
        pos_sym, pos_cxr, pos_xpert, eli = pars['pos_sym'], pars['pos_cxr'], pars['pos_xpert'], pars['eli']

        I = self.Keys
        arrived = r_acf * y
        screened = arrived * eli
        confirmed = screened * (1 - (1 - pos_cxr) * (1 - pos_sym))
        pos = confirmed * pos_xpert
        tp_ds, tp_dr = pos[I.Infectious_DS], pos[I.Infectious_DR]
        tp_dr_fl = tp_dr * (1 - p_dst)
        tp_dr_sl = tp_dr - tp_dr_fl

        if mea:
            n = y.sum()
            return _mea_acf(n, arrived, screened, confirmed, pos, tp_ds, tp_dr_fl, tp_dr_sl, 'MDU',
                            n_sym=screened.sum(), n_cxr=screened.sum(), n_xpert=confirmed.sum())
        else:
            dy = np.zeros_like(y)
            dy[I.Infectious_DS, :2] -= tp_ds
            dy[I.Infectious_DR, :2] -= tp_dr
            dy[I.Txf_Pub_DS, :2] += tp_ds.sum(0)
            dy[I.Txf_Pub_DR, :2] += tp_dr_fl.sum(0)
            dy[I.Txs_Pub_DR, :2] += tp_dr_sl.sum(0)
            return dy

    def _calc_d2d(self, y, r_acf, p_dst, pars, mea=False):
        pos_sym, pos_xpert, eli = pars['pos_sym'], pars['pos_xpert'], pars['eli']

        I = self.Keys
        arrived = r_acf * y
        screened = arrived * eli
        confirmed = screened * pos_sym
        pos = confirmed * pos_xpert
        tp_ds, tp_dr = pos[I.Infectious_DS], pos[I.Infectious_DR]
        tp_dr_fl = tp_dr * (1 - p_dst)
        tp_dr_sl = tp_dr - tp_dr_fl

        if mea:
            n = y.sum()
            return _mea_acf(n, arrived, screened, confirmed, pos, tp_ds, tp_dr_fl, tp_dr_sl, 'D2D',
                            n_sym=screened.sum(), n_xpert=confirmed.sum())
        else:
            dy = np.zeros_like(y)
            dy[I.Infectious_DS, :2] -= tp_ds
            dy[I.Infectious_DR, :2] -= tp_dr
            dy[I.Txf_Pub_DS, :2] += tp_ds.sum(0)
            dy[I.Txf_Pub_DR, :2] += tp_dr_fl.sum(0)
            dy[I.Txs_Pub_DR, :2] += tp_dr_sl.sum(0)
            return dy

    def calc_dy(self, t, y, pars, intv):
        r_acf_mdu, r_acf_d2d, p_dst = pars['r_acf_mdu'], pars['r_acf_d2d'], pars['acf_dst_sens']
        r_acf_mdu, r_acf_d2d, p_dst = intv.modify_acf_bg(t, r_acf_mdu, r_acf_d2d, p_dst)

        dy = np.zeros_like(y)
        # MDU
        if r_acf_mdu > 0:
            dy[:, :2] += self._calc_mdu(y[:, :2], r_acf_mdu, p_dst, pars, mea=False)
        # D2D
        if r_acf_d2d > 0:
            dy[:, :2] += self._calc_d2d(y[:, :2], r_acf_d2d, p_dst, pars, mea=False)
        return dy

    def measure(self, t, y, pars, intv, mea):
        r_acf_mdu, r_acf_d2d, p_dst = pars['r_acf_mdu'], pars['r_acf_d2d'], pars['acf_dst_sens']
        r_acf_mdu, r_acf_d2d, p_dst = intv.modify_acf_bg(t, r_acf_mdu, r_acf_d2d, p_dst)

        # MDU
        mea.update(self._calc_mdu(y[:, :2], r_acf_mdu, p_dst, pars, mea=True))

        # D2D
        mea.update(self._calc_d2d(y[:, :2], r_acf_d2d, p_dst, pars, mea=True))

        l, lt = y[self.Keys.LTBI0].sum(), y[self.Keys.LTBI_TPT].sum()
        mea['PrOnPseudoTPT'] = lt / (l + lt)


class VulACF:
    def __init__(self, keys):
        self.Keys = keys
        self.Recap = True

    def _calc_vul_acf(self, y, r_acf0, r_fu, r_lost, pars, mea=False):
        I = self.Keys
        pos_cxr, pos_xpert, p_dst = pars['pos_cxr'], pars['pos_xpert'], pars['acf_dst_sens']
        pos_vul, pos_sym, eli = pars['pos_vul'], pars['pos_sym'], pars['eli']

        n = y.sum()

        n_target = r_acf0 * y.sum()
        r_acf = min(24, n_target / (eli * y).sum())

        arrived = r_acf * y
        screened0 = arrived * eli
        screened = screened0 * (1 - (1 - pos_sym) * (1 - pos_vul))
        confirmed = screened * pos_cxr
        pos = confirmed * pos_xpert
        neg = screened - pos
        tp_ds, tp_dr = pos[I.Infectious_DS], pos[I.Infectious_DR]
        tp_dr_fl = tp_dr * (1 - p_dst)
        tp_dr_sl = tp_dr - tp_dr_fl

        if mea:
            m = _mea_acf(n, arrived, screened, confirmed, pos, tp_ds, tp_dr_fl, tp_dr_sl, 'Vul',
                         n_vul=screened0.sum(), n_sym=screened0.sum(), n_cxr=screened.sum(), n_xpert=confirmed.sum())
        else:
            dy = np.zeros_like(y)
            dy[I.Infectious_DS] -= tp_ds
            dy[I.Infectious_DR] -= tp_dr
            dy[I.Txf_Pub_DS] += tp_ds.sum(0)
            dy[I.Txf_Pub_DR] += tp_dr_fl.sum(0)
            dy[I.Txs_Pub_DR] += tp_dr_sl.sum(0)

            # to follow up list
            if r_lost > 0:
                dy[:, 1] -= neg[:, 1]
                dy[:, 3] += neg[:, 1]

                lost = r_lost * y[:, 2:]
                dy[:, 2:] -= lost
                dy[:, :2] += lost

                complete = y[I.LTBI_TPT] * 2
                dy[I.LTBI_TPT] -= complete
                dy[I.LTBI0, :2] += complete[:, :2] + complete[:, 2:]

            fp = pos[I.LTBI0]
            dy[I.LTBI0] -= fp
            dy[I.LTBI_TPT] += fp

        if r_lost <= 0:
            if mea:
                return m
            else:
                return dy

        # Follow-up
        arrived = r_fu * y[:, 3]
        screened = arrived
        confirmed = screened * pos_sym
        pos = confirmed * pos_xpert
        tp_ds, tp_dr = pos[I.Infectious_DS], pos[I.Infectious_DR]
        tp_dr_fl = tp_dr * (1 - p_dst)
        tp_dr_sl = tp_dr - tp_dr_fl

        if mea:
            m.update(_mea_acf(n, arrived, screened, confirmed, pos, tp_ds, tp_dr_fl, tp_dr_sl, 'VulFu'))
        else:
            dy[I.Infectious_DS, 3] -= tp_ds
            dy[I.Infectious_DR, 3] -= tp_dr
            dy[I.Txf_Pub_DS, 3] += tp_ds.sum()
            dy[I.Txf_Pub_DR, 3] += tp_dr_fl.sum()
            dy[I.Txs_Pub_DR, 3] += tp_dr_sl.sum()

            fp = pos[I.LTBI0]
            dy[I.LTBI0, 3] -= fp
            dy[I.LTBI_TPT, 3] += fp
        if mea:
            return m
        else:
            return dy

    def _calc_plain_acf(self, y, r_acf0, pars, mea=False):
        I = self.Keys
        pos_cxr, pos_xpert, p_dst = pars['pos_cxr'], pars['pos_xpert'], pars['acf_dst_sens']
        eli = pars['eli']

        n = y.sum()

        n_target = r_acf0 * y.sum()
        r_acf = min(24, n_target / (eli * y).sum())

        arrived = r_acf * y
        screened = arrived * eli
        confirmed = screened * pos_cxr
        pos = confirmed * pos_xpert

        tp_ds, tp_dr = pos[I.Infectious_DS], pos[I.Infectious_DR]
        tp_dr_fl = tp_dr * (1 - p_dst)
        tp_dr_sl = tp_dr - tp_dr_fl

        if mea:
            m = _mea_acf(n, arrived, screened, confirmed, pos, tp_ds, tp_dr_fl, tp_dr_sl, 'Plain',
                         n_cxr=screened.sum(), n_xpert=confirmed.sum())
        else:
            dy = np.zeros_like(y)
            dy[I.Infectious_DS] -= tp_ds
            dy[I.Infectious_DR] -= tp_dr
            dy[I.Txf_Pub_DS] += tp_ds.sum(0)
            dy[I.Txf_Pub_DR] += tp_dr_fl.sum(0)
            dy[I.Txs_Pub_DR] += tp_dr_sl.sum(0)

            fp = pos[I.LTBI0]
            dy[I.LTBI0] -= fp
            dy[I.LTBI_TPT] += fp

        if mea:
            return m
        else:
            return dy

    def calc_dy(self, t, y, pars, intv):
        dy = np.zeros_like(y)
        # Vulnerability-led ACF
        r_acf_vul, r_fu_v, r_lost_v = intv.modify_acf_vul(t, 0, 0, 0)
        if r_acf_vul > 0:
            dy += self._calc_vul_acf(y, r_acf_vul, r_fu_v, r_lost_v, pars, mea=False)

        # Plain ACF
        r_acf_plain = intv.modify_acf_plain(t, 0)
        if r_acf_plain > 0:
            dy += self._calc_plain_acf(y, r_acf_plain, pars, mea=False)
        return dy

    def measure(self, t, y, pars, intv, mea):
        r_acf_vul, r_fu, r_lost = intv.modify_acf_vul(t, 0, 0, 0)
        mea.update(self._calc_vul_acf(y, r_acf_vul, r_fu, r_lost, pars, mea=True))

        # Plain ACF
        r_acf_plain = intv.modify_acf_plain(t, 0)
        mea.update(self._calc_plain_acf(y, r_acf_plain, pars, mea=True))
