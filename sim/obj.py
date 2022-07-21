from sims_pars.fitting.base import AbsObjectiveSimBased
from itertools import groupby
import pandas as pd
import numpy as np
import json
import time


__author__ = 'Chu-Chang Ku'
__all__ = ['Objective']


class Objective(AbsObjectiveSimBased):
    def __init__(self, model, bn, filepath_targets, or_prev=1, exo=None):
        AbsObjectiveSimBased.__init__(self, bn, exo=exo)
        self.Model = model
        self.Data = dict()

        targets = json.load(open(filepath_targets, 'r'))
        if or_prev is not None:
            targets['All'].append({
              "Year": 2020,
              "Index": "OR_prev_comorb",
              "Tag": "All",
              "M": or_prev,
              "L": or_prev - 0.1,
              "U": or_prev + 0.1
            })
        for gp, ds in targets.items():
            self.Data[gp] = dict()
            for k, vs in groupby(ds, key=lambda x: x['Index']):
                self.Data[gp][k] = pd.Series({v['Year']: v['M'] for v in vs})

        self.Map = list()

        for ent in targets['All']:
            scale = abs(ent['U'] - ent['L']) / 4
            self.Map.append((ent['Index'], self.Data['All'][ent['Index']], scale))

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

    to_fit = Objective(model=Model(inputs, year0=1970),
                       bn=get_bn('../prior'),
                       filepath_targets='../data/targets.json')

    pars = to_fit.sample_prior()
    print(pars)
    sim = to_fit.simulate(pars)

    print(sim[1].head())

    li = to_fit.link_likelihood(sim)
    print(li)
