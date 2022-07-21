import pymc as pm
from abc import ABCMeta, abstractmethod

__all__ = ['DataModel', 'post_to_particles', 'prior_to_particles']


class DataModel(metaclass=ABCMeta):
    def __init__(self, obs, eps):
        self.Model = pm.Model()

        with self.Model as dm:
            ps = self.define_prior(dm)

        self.FreeParameters = [p.name for p in ps]

        def to_fit(rng, *args, size=None):
            pars = {k: v for k, v in zip(self.FreeParameters, args)}
            return self.simulate(pars)

        with self.Model:
            pm.Simulator("sim", to_fit, params=ps, epsilon=eps, observed=obs)

    @abstractmethod
    def define_prior(self, dm):
        pass

    @abstractmethod
    def simulate(self, pars):
        pass


def post_to_particles(post):
    po = post.posterior.stack(samples=("draw", "chain"))

    variables = list(po.variables.keys())
    variables = [v for v in variables if v != 'samples']
    vs = {k: po[k].to_numpy() for k in variables}
    n_samples = len(vs[variables[0]])
    pts = [{k: v[i] for k, v in vs.items()} for i in range(n_samples)]
    return pts


def prior_to_particles(prior):
    po = prior.prior.stack(samples=("draw", "chain"))

    variables = list(po.variables.keys())
    variables = [v for v in variables if v != 'samples']
    vs = {k: po[k].to_numpy() for k in variables}
    n_samples = len(vs[variables[0]])
    pts = [{k: v[i] for k, v in vs.items()} for i in range(n_samples)]
    return pts
