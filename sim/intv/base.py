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


class ACFPlain(BaseModel):
    R_ACF: confloat(ge=0) = 0
    Type: str = 'mod'
    Focus: bool = True

#
# class ACF(BaseModel):
#     Yield: confloat(ge=0, le=200) = 0
#     Type: str = 'mod'
#     Asym: bool = False
#     HiRisk: bool = False


class Intervention(BaseModel):
    ACF: ACF = ACF()
    ACFPlain: ACFPlain = ACFPlain()
    T0_Intv: float = 2020
    T1_Intv: float = 2023

    def modify_acf(self, t, r_acf0, r_acf, p_dst, pars):
        if self.ACFPlain.R_ACF > 0:
            r_acf0 = self.ACFPlain.R_ACF
            type = self.ACFPlain.Type

            sens = np.array([
                pars[f'sens_acf_sn_{type}'], pars[f'sens_acf_sp_{type}'],
                pars[f'sens_acf_sn_{type}'], pars[f'sens_acf_sp_{type}']
            ]).reshape((-1, 1))
            p_dst = pars[f'p_dst_acf_{type}']

            p_dst = p_dst
            r_acf = r_acf0 * sens

            if np.any(r_acf > 0):
                if self.ACFPlain.Focus:
                    r_acf_low = np.zeros((4, 1))
                else:
                    r_acf_low = r_acf
                r_acf = np.concatenate([r_acf_low, r_acf], axis=1)

            return r_acf0, r_acf, p_dst

        if t > self.T0_Intv and self.ACF.Scale > 0:
            wt = scale_up(t, self.T0_Intv, self.T1_Intv) * self.ACF.Scale
            r_acf0 = wt
            type = self.ACF.Type

            sens = np.array([
                pars[f'sens_acf_sn_{type}'], pars[f'sens_acf_sp_{type}'],
                pars[f'sens_acf_sn_{type}'], pars[f'sens_acf_sp_{type}']
            ]).reshape((-1, 1))
            p_dst = pars[f'p_dst_acf_{type}']

            p_dst = wt * p_dst
            r_acf = wt * sens #* pars['sens_acf_screen']

        return r_acf0, r_acf, p_dst


if __name__ == '__main__':
    intv_list = {
        'ACF': {'Scale': 0.5},
    }

    intv = Intervention.parse_obj(intv_list)
    print(intv.json())


