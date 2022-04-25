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


class ImpDx(BaseModel):
    Dx: confloat(ge=0, le=1) = 0


class CS(BaseModel):
    Scale: confloat(ge=0, le=1) = 0


class Intervention(BaseModel):
    ImpDx: ImpDx = ImpDx()
    CS: CS = CS()
    T0_Intv: float = 2022
    T1_Intv: float = 2025

    def modify_dx(self, t, p_dx):
        if t > self.T0_Intv and self.ImpDx.Dx > 0:
            wt = scale_up(t, self.T0_Intv, self.T1_Intv)
            p_dx1 = self.ImpDx.Dx

            p_dx = p_dx + (p_dx1 - p_dx) * wt
        return p_dx

    def modify_cs(self, t, r_cs, r_rcs):
        if t > self.T0_Intv and self.CS.Scale > 0:
            wt = scale_up(t, self.T0_Intv, self.T1_Intv)
            r_cs1, r_rcs1 = r_cs.copy(), r_rcs.copy()

            r_cs1 /= 1 - self.CS.Scale
            r_rcs1 /= 1 - self.CS.Scale

            r_cs = r_cs + (r_cs1 - r_cs) * wt
            r_rcs = r_rcs + (r_rcs1 - r_rcs) * wt
        return r_cs, r_rcs


if __name__ == '__main__':
    intv_list = {
        'ImpDx': {'Dx': 0.90},
        'CS': {'Scale': 0.5},
    }

    intv = Intervention.parse_obj(intv_list)
    print(intv.json())


