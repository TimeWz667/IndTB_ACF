from pydantic import BaseModel
from pydantic.types import confloat
import numpy as np

__author__ = 'Chu-Chang Ku'
__all__ = ['Intervention']


class FullACF(BaseModel):
    Coverage: confloat(ge=0) = 0
    ScreenAlg: str = 'VSC'  # Screening algorithm
    FollowUp: confloat(ge=0) = 0  # Number of follow-up per year
    Duration: confloat(ge=0) = 0  # Follow-up period
    Year0: confloat(ge=2022) = 2023


class AltACF(BaseModel):
    Coverage: confloat(ge=0) = 0
    ScreenAlg: str = 'Sy'  # Screening algorithm
    Year0: confloat(ge=2022) = 2023


class Intervention(BaseModel):
    FullACF = FullACF()
    AltACF = AltACF()

    def find_rates_main(self, t, p_eli):
        alg = self.FullACF.ScreenAlg

        if t >= self.FullACF.Year0:
            r_acf = self.FullACF.Coverage  # / p_eli
            r_acf = min(r_acf, 26)
            if self.FullACF.Duration <= 0:
                r_loss, r_fu = 0, 0
            else:
                r_loss = 1 / self.FullACF.Duration
                if self.FullACF.FollowUp <= 0:
                    r_fu = 0
                else:
                    r_fu = self.FullACF.FollowUp
            return r_acf, r_loss, r_fu, alg
        return 0, 0, 0, alg

    def find_rates_alt(self, t, p_eli):
        alg = self.AltACF.ScreenAlg

        if t >= self.AltACF.Year0:
            r_acf = self.AltACF.Coverage  # / p_eli
            r_acf = min(r_acf, 26)
            return r_acf, alg
        return 0, alg


if __name__ == '__main__':
    intv_list = {
        'FullACF': {'Coverage': 0.5},
        'AltACF': {'Coverage': 0.2}
    }

    intv = Intervention.parse_obj(intv_list)
    print(intv.json())
