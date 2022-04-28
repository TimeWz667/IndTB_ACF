from pydantic import BaseModel
from pydantic.types import confloat
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['Intervention']


def scale_up(t, t0, t1):
    if t < t0:
        return 0
    if t > t1:
        return 1
    return (t - t0) / (t1 - t0)


class ACF(BaseModel):
    Scale: confloat(ge=0, le=1) = 0
    Type: str = 'mod'
    SensDx: confloat(ge=0, le=1) = 0
    PrDST: confloat(ge=0, le=1) = 0
    SensScreen: confloat(ge=0, le=1) = 0


class Intervention(BaseModel):
    ACF: ACF = ACF()
    T0_Intv: float = 2022
    T1_Intv: float = 2025

    def modify_acf(self, t, r_acf, p_dst, pars):
        if t > self.T0_Intv and self.ACF.Scale > 0:
            wt = scale_up(t, self.T0_Intv, self.T1_Intv) * self.ACF.Scale * pars['sens_acf_screen']

            if self.ACF.Type == 'mod':
                sens = np.array([
                    pars['sens_acf_sn_mod'], pars['sens_acf_sp_mod'],
                    pars['sens_acf_sn_mod'], pars['sens_acf_sp_mod']
                ]).reshape((-1, 1))
                p_dst = pars['p_dst_acf_mod']
            else:
                sens = np.array([
                    pars['sens_acf_sn_high'], pars['sens_acf_sp_high'],
                    pars['sens_acf_sn_high'], pars['sens_acf_sp_high']
                ]).reshape((-1, 1))
                p_dst = pars['p_dst_acf_high']
            p_dst = wt * p_dst
            r_acf = wt * sens

        return r_acf, p_dst


if __name__ == '__main__':
    intv_list = {
        'ACF': {'Scale': 0.5},
    }

    intv = Intervention.parse_obj(intv_list)
    print(intv.json())


