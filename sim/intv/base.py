from pydantic import BaseModel
from pydantic.types import confloat
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['Intervention']


class MU(BaseModel):
    Scale: confloat(ge=0) = 0


class D2D(BaseModel):
    Scale: confloat(ge=0) = 0


class PlainACF(BaseModel):
    Coverage: confloat(ge=0) = 0
    FollowUp: confloat(ge=0) = 0
    Duration: confloat(ge=0) = 0


class VulACF(BaseModel):
    Coverage: confloat(ge=0) = 0
    FollowUp: confloat(ge=0) = 0
    Duration: confloat(ge=0) = 0


class Intervention(BaseModel):
    VulACF = VulACF()
    MU = MU()
    D2D = D2D()
    PlainACF = PlainACF()
    T0_Bg = 2022
    T0_Vul = 2023

    def modify_acf_bg(self, t, r_acf_mu, r_acf_d2d, p_dst):
        if t >= self.T0_Bg:
            r_acf_mu *= self.MU.Scale
            r_acf_d2d *= self.D2D.Scale
        return r_acf_mu, r_acf_d2d, p_dst

    def modify_acf_vul(self, t, cov, r_fu, r_lost):
        if t >= self.T0_Vul:
            cov = self.VulACF.Coverage
            if self.VulACF.Duration > 0:
                r_fu = 1 / max(self.VulACF.FollowUp, self.VulACF.Duration)
                r_lost = 1 / self.VulACF.Duration
        return cov, r_fu, r_lost

    def modify_acf_plain(self, t, cov, r_fu, r_lost):
        if t >= self.T0_Vul:
            cov = self.PlainACF.Coverage
            if self.PlainACF.Duration > 0:
                r_fu = 1 / max(self.PlainACF.FollowUp, self.PlainACF.Duration)
                r_lost = 1 / self.PlainACF.Duration
        return cov, r_fu, r_lost


if __name__ == '__main__':
    intv_list = {
        'VulACF': {'Scale': 0.5},
    }

    intv = Intervention.parse_obj(intv_list)
    print(intv.json())
