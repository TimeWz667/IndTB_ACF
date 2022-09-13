from abc import ABCMeta, abstractmethod

__author__ = 'Chu-Chang Ku'
__all__ = ['Process']


class Process(metaclass=ABCMeta):
    def __init__(self, keys):
        self.Keys = keys

    @abstractmethod
    def calc_dy(self, t, y, pars):
        pass

    @abstractmethod
    def measure(self, t, y, pars, mea):
        pass
