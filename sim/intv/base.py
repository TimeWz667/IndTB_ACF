from pydantic import BaseModel
from pydantic.types import confloat
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['Intervention']


class VulACF(BaseModel):
    Coverage: confloat(ge=0, le=2) = 0


class MU(BaseModel):
    Scale: confloat(ge=0) = 1


class D2D(BaseModel):
    Scale: confloat(ge=0) = 1


class PlainACF(BaseModel):
    Coverage: confloat(ge=0) = 0


class Intervention(BaseModel):
    VulACF = VulACF()
    MU = MU()
    D2D = D2D()
    PlainACF = PlainACF()
    T0_Intv = 2022

    def modify_acf_bg(self, t, r_acf_mu, r_acf_d2d, p_dst):
        if t >= self.T0_Intv:
            r_acf_mu *= self.MU.Scale
            r_acf_d2d *= self.D2D.Scale

        return r_acf_mu, r_acf_d2d, p_dst

    def modify_acf_vul(self, t, r_acf):
        if t >= 2023:
            r_acf = self.VulACF.Coverage

        return r_acf


if __name__ == '__main__':
    intv_list = {
        'ACF': {'Scale': 0.5},
    }

    intv = Intervention.parse_obj(intv_list)
    print(intv.json())
