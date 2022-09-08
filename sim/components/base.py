from abc import ABCMeta, abstractmethod

__author__ = 'Chu-Chang Ku'
__all__ = ['Process']


class Process(metaclass=ABCMeta):
    def __init__(self, keys):
        self.Keys = keys

    # def __call__(self, t, y, pars, intv, calc):
    #     pass

    @abstractmethod
    def calc_dy(self, t, y, pars, intv):
        pass

    @abstractmethod
    def measure(self, t, y, pars, intv, mea):
        pass
