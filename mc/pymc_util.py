import pymc as pm

__all__ = ['DataModel', 'post_to_particles', 'prior_to_particles']


class DataModel:
    DefaultExo = {
        'r_die': 0.015151515151515152,
        'r_growth': 0.0175,
        'p_primary': 0.14,
        'r_react': 0.001,
        'r_stab': 0.5,
        'p_sp0': 0.6,
        'r_sc': 0.2,
        'r_die_ut': 0.2,
        'p_entry_pri': 0.5,
        'p_dx_pub': 0.83,
        'p_dx_pri': 0.7,
        'p_txf_pub': 0.88, 'p_txf_pri': 0.7,
        'p_txs_pub': 0.88, 'p_txs_pri': 0,
        'p_dst_pcf': 0.12,
        'sens_acf_sp_high': 1, 'sens_acf_sn_high': 0.7, 'spec_acf_high': 0.99, 'p_dst_acf_high': 0.95,
        'sens_acf_sp_mod': 1, 'sens_acf_sn_mod': 0, 'spec_acf_mod': 0.98, 'p_dst_acf_mod': 0.12,
        'sens_acf_screen': 0.7, 'spec_acf_screen': 0.61,
        'r_acf_screen': 365,
        'r_acf_dx': 52,
        'r_succ_txf': 2,
        'r_succ_txs': 0.5,
        'p_succ_txf_pub': 0.85, 'p_succ_txf_pri': 0.6,
        'p_succ_txs_pub': 0.46, 'p_succ_txs_pri': 0,
        'p_tr_pub': 0.88,
        'r_mdr_tx': 0.1,
        'r_succ_txf_pub': 2, 'r_succ_txf_pri': 2,
        'r_succ_txs_pub': 0.5, 'r_succ_txs_pri': 0,
        'r_ltfu_txf_pub': 0.3529411764705883, 'r_ltfu_txf_pri': 1.3333333333333335,
        'r_ltfu_txs_pub': 0.5869565217391305, 'r_ltfu_txs_pri': 0.5,
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

    def define_prior(self, dm):
        with dm:
            p_comorb = pm.Uniform('p_comorb', 0, 0.5)
            rr_risk_comorb = pm.Uniform('rr_risk_comorb', 1, 30)
            beta_ds = pm.Uniform('beta_ds', 1, 15)
            rr_beta_dr = pm.Uniform('rr_beta_dr', 0.9, 1.1)
            rr_inf_asym = pm.Uniform('rr_inf_asym', 0, 1)
            rr_inf_sn = pm.Uniform('rr_inf_sn', 0.1, 0.3)
            red_sus = pm.Uniform('red_sus', 0.15, 0.25)
            rr_sus_ltbi = pm.Uniform('rr_sus_ltbi', 0.15, 0.25)
            r_relapse_td = pm.Triangular('r_relapse_td', lower=0.105, c=0.14, upper=0.175)
            r_relapse_tc = pm.Triangular('r_relapse_tc', lower=0.024, c=0.032, upper=0.04)
            r_relapse_st = pm.Triangular('r_relapse_st', lower=0.0011, c=0.0019, upper=0.002)

            r_onset_sp = pm.Triangular('r_onset_sp', lower=1.02, c=1.24, upper=1.65)
            r_onset_sn = pm.Triangular('r_onset_sn', lower=1.9, c=2.37, upper=3.05)

            r_convert_a = pm.Triangular('r_convert_a', lower=0.62, c=0.63, upper=0.64)
            r_convert_s = pm.Triangular('r_convert_s', lower=0.4, c=0.71, upper=1.04)

            r_clear = pm.Uniform('r_clear', 0.02, 0.04)

            r_cs_s = pm.Triangular('r_cs_s', lower=0.57, upper=0.91, c=0.73)
            r_cs_c = pm.Triangular('r_cs_c', lower=9, upper=15, c=12)

            return (p_comorb, rr_risk_comorb, beta_ds, rr_beta_dr, rr_inf_asym, rr_inf_sn,
                    red_sus, rr_sus_ltbi,
                    r_relapse_td, r_relapse_tc, r_relapse_st,
                    r_onset_sp, r_onset_sn, r_convert_a, r_convert_s,
                    r_clear, r_cs_s, r_cs_c)

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


def post_to_particles(post):
    po = post.posterior.stack(samples=("draw", "chain"))

    variables = list(po.variables.keys())
    variables = [v for v in variables if v != 'samples']
    vs = {k: po[k].to_numpy() for k in variables}
    n_samples = len(vs[variables[0]])
    pts = [{k: v[i] for k, v in vs.items() if k not in ['chain', 'draw']} for i in range(n_samples)]
    return pts


def prior_to_particles(prior):
    po = prior.prior.stack(samples=("draw", "chain"))

    variables = list(po.variables.keys())
    variables = [v for v in variables if v != 'samples']
    vs = {k: po[k].to_numpy() for k in variables}
    n_samples = len(vs[variables[0]])
    pts = [{k: v[i] for k, v in vs.items() if k not in ['chain', 'draw']} for i in range(n_samples)]
    return pts
