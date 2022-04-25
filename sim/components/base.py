__author__ = 'Chu-Chang Ku'
__all__ = ['Process']


class Process:
    def __init__(self, keys, intv=None):
        self.Keys = keys
        self.Intervention = intv

    def __call__(self, t, y, pars, calc):
        pass

    def measure(self, t, y, pars, calc, mea):
        pass
