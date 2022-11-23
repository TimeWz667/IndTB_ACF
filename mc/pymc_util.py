import pymc as pm

__all__ = ['DataModel', 'post_to_particles', 'prior_to_particles']


class DataModel:
    DefaultExo = {
        'r_die': 1 / 70,
        'r_growth': 0.0167,
        'r_die_ut': 0.127,

        'p_primary': 0.14,

        'r_stab': 0.5,
        'r_sc': 0.206,

        'p_txi': 0.95,
        'r_succ_txf': 2,
        'r_succ_txs': 0.5,
        'p_succ_txf_pub': 0.83, 'p_succ_txf_pri': 0.83,
        'p_succ_txs_pub': 0.7, 'p_succ_txs_pri': 0,
        'p_die_txf_ds': 0.05, 'p_die_txs_dr': 0.14,

        'r_mdr_tx': 0.05,
        'acf_dst_sens': 1
    }

    def __init__(self, obs, eps, simulator, exo=None):
        self.Obs = obs
        self.Eps = eps
        self.Exo = dict(DataModel.DefaultExo)
        self.Simulator = simulator
        if exo is not None:
            self.Exo.update(exo)

        with pm.Model() as dm:
            ps = self.define_prior(dm)

        self.FreeParameters = [p.name for p in ps]

    @staticmethod
    def define_prior(dm):
        with dm:
            p_comorb = pm.Uniform('p_comorb', 0, 0.5)
            rr_risk_comorb = pm.Uniform('rr_risk_comorb', 1, 30)
            beta_ds = pm.Uniform('beta_ds', 1, 40)
            rr_beta_dr = pm.Uniform('rr_beta_dr', 0.9, 1.1)
            rr_inf_asym = pm.Uniform('rr_inf_asym', 0, 1)

            rr_sus_ltbi = pm.Uniform('rr_sus_ltbi', 0.25, 0.75)
            r_react = pm.Uniform('r_react', 0.0005, 0.0015)
            r_relapse_td = pm.Triangular('r_relapse_td', lower=0.105, c=0.14, upper=0.175)
            r_relapse_tc = pm.Triangular('r_relapse_tc', lower=0.024, c=0.032, upper=0.04)
            r_relapse_st = pm.Triangular('r_relapse_st', lower=0.0011, c=0.0019, upper=0.002)

            r_onset = pm.Uniform('r_onset', lower=1, upper=6)

            r_clear = pm.Uniform('r_clear', 0.02, 0.04)

            r_cs_s = pm.Uniform('r_cs_s', lower=1, upper=15)
            r_cs_c = pm.Uniform('r_cs_c', lower=1, upper=15)
            p_entry_pri = pm.Triangular('p_entry_pri', lower=0.4, c=0.5, upper=0.6)
            p_dx_pub = pm.Triangular('p_dx_pub', lower=0.81, c=0.83, upper=0.85)
            p_dx_pri = pm.Uniform('p_dx_pri', lower=0.5, upper=0.8)
            p_dst_pcf = pm.Triangular('p_dst_pcf', lower=0.08, c=0.12, upper=0.2)

            p_tr_pub = pm.Triangular('p_tr_pub', lower=0.85, c=0.88, upper=0.92)

            acf_cxr_sens = pm.Triangular('acf_cxr_sens', lower=0.95, c=0.98, upper=1)
            acf_cxr_spec = pm.Triangular('acf_cxr_spec', lower=0.72, c=0.75, upper=0.79)
            acf_xpert_sens = pm.Triangular('acf_xpert_sens', lower=0.9, c=0.92, upper=0.94)
            acf_xpert_spec = pm.Triangular('acf_xpert_spec', lower=0.9, c=0.99, upper=1)
            acf_dst_sens = pm.Triangular('acf_dst_sens', lower=0.9, c=0.95, upper=0.97)

            return (p_comorb, rr_risk_comorb, beta_ds, rr_beta_dr, rr_inf_asym,
                    rr_sus_ltbi,
                    r_react, r_relapse_td, r_relapse_tc, r_relapse_st,
                    r_onset, r_clear,
                    r_cs_s, r_cs_c, p_entry_pri, p_dx_pub, p_dx_pri, p_dst_pcf, p_tr_pub,
                    acf_cxr_sens, acf_cxr_spec, acf_xpert_sens, acf_xpert_spec, acf_dst_sens
                    )

    def serve(self, pars):
        pars = dict(pars)
        pars.update(self.Exo)
        return pars

    def simulate(self, pars):
        pars = self.serve(pars)
        return self.Simulator(pars)

    def build_model(self):
        dm = pm.Model()
        ps = self.define_prior(dm)

        def to_fit(rng, *args, size=None):
            pars = {k: v for k, v in zip(self.FreeParameters, args)}
            return self.simulate(pars)

        with dm:
            pm.Simulator("sim", to_fit, params=ps, epsilon=self.Eps, observed=self.Obs)
        return dm


def _to_particles(po):
    variables = list(po.variables.keys())
    variables = [v for v in variables if v not in ['samples', 'draw', 'chain']]
    vs = {k: po[k].to_numpy() for k in variables}
    n_samples = len(vs[variables[0]])
    pts = [{k: v[i] for k, v in vs.items()} for i in range(n_samples)]
    return pts


def post_to_particles(post):
    po = post.posterior.stack(samples=("draw", "chain"))
    return _to_particles(po)


def prior_to_particles(prior):
    po = prior.prior.stack(samples=("draw", "chain"))
    return _to_particles(po)
