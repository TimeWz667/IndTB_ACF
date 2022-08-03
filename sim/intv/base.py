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


# class ACF(BaseModel):
#     Scale: confloat(ge=0, le=1) = 0
#     Type: str = 'mod'
#     SensDx: confloat(ge=0, le=1) = 0
#     PrDST: confloat(ge=0, le=1) = 0
#     SensScreen: confloat(ge=0, le=1) = 0

#
# class ACFPlain(BaseModel):
#     R_ACF: confloat(ge=0) = 0
#     Type: str = 'mod'
#     Focus: bool = True


class ACF(BaseModel):
    Yield: confloat(ge=0, le=2) = 0
    Type: str = 'mod'
    HiRisk: bool = False


class Intervention(BaseModel):
    ACF: ACF = ACF()
    # ACFPlain: ACFPlain = ACFPlain()
    T0_Intv: float = 2022
    T1_Intv: float = 2025

    def modify_acf(self, t, r_acf, r_acf_tp, r_acf_fp, p_dst, pars, p_tb, p_nontb):
        if t >= self.T0_Intv and self.ACF.Yield > 0:
            wt = scale_up(t, self.T0_Intv, self.T1_Intv)
            n2detect = self.ACF.Yield * wt

            if self.ACF.HiRisk:
                detectable = p_tb[1] + p_nontb[1]
                r_acf = n2detect / detectable * np.array([0, 1])
            else:
                detectable = p_tb.sum() + p_nontb.sum()
                r_acf = n2detect / detectable * np.ones(2)

            type = self.ACF.Type

            sens = np.array([pars[f'sens_acf_sn_{type}'], pars[f'sens_acf_sp_{type}']])

            spec = pars[f'spec_acf_{type}']
            p_dst = pars[f'p_dst_acf_{type}']

            p_dst = wt * p_dst
            r_acf_tp = r_acf.reshape((-1, 2)) * sens.reshape((2, -1))
            r_acf_tp[r_acf_tp > 20] = 20

            r_acf_fp = r_acf * (1 - spec)

        return r_acf, r_acf_tp, r_acf_fp, p_dst


if __name__ == '__main__':
    intv_list = {
        'ACF': {'Scale': 0.5},
    }

    intv = Intervention.parse_obj(intv_list)
    print(intv.json())


