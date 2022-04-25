from sims_pars.fitting.base import AbsObjectiveSimBased
from itertools import groupby
import pandas as pd
from scipy.stats import norm
import numpy as np
import json
import time


__author__ = 'Chu-Chang Ku'
__all__ = ['Objective']


class Objective(AbsObjectiveSimBased):
    def __init__(self, model, bn, filepath_targets, exo=None):
        AbsObjectiveSimBased.__init__(self, bn, exo=exo)
        self.Model = model
        self.Data = dict()

        targets = json.load(open(filepath_targets, 'r'))
        for gp, ds in targets.items():
            self.Data[gp] = dict()
            for k, vs in groupby(ds, key=lambda x: x['Index']):
                self.Data[gp][k] = pd.Series({v['Year']: v['M'] for v in vs})

        self.Map = [
            ('IncR', self.Data['All']['Inc'], 1e-3),
            ('MorR', self.Data['All']['Mor'], 1e-4),
        ]

    def simulate(self, pars):
        time.sleep(0.001)
        return self.Model.simulate(pars)

    def link_likelihood(self, sim):
        ys, ms, msg = sim
        if not msg['succ'] or not ys.success:
            return - np.inf

        li = 0
        for sim, data, scale in self.Map:
            try:
                li += ((ms[sim] - data).dropna() / scale).pow(2).sum()
            except KeyError:
                pass

        return -float(li)


if __name__ == '__main__':
    from sim.inputs import load_inputs
    from sim.dy import Model, get_bn

    inputs = load_inputs('../data/pars.json')

    to_fit = Objective(model=Model(inputs),
                       bn=get_bn('../prior'),
                       filepath_targets='../data/targets.json')

    pars = to_fit.sample_prior()
    print(pars)
    sim = to_fit.simulate(pars)

    print(sim[1].head())

    li = to_fit.link_likelihood(sim)
    print(li)
