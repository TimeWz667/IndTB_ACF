import os
import re
from sim.intv import Intervention


ds = [d for d in os.listdir("../out") if d.startswith('dy_')]

# print(ds)
# for path in ds:

intvs = dict()

intvs['Baseline'] = {}
intvs['High, Cov=100%'] = {'ACF': {'Type': 'high', 'Scale': 1}}
intvs['High, Cov=50%'] = {'ACF': {'Type': 'high', 'Scale': 0.5}}
intvs['High, Cov=20%'] = {'ACF': {'Type': 'high', 'Scale': 0.2}}
intvs['Mod, Cov=100%'] = {'ACF': {'Type': 'mod', 'Scale': 1}}
intvs['Mod, Cov=50%'] = {'ACF': {'Type': 'mod', 'Scale': 0.5}}
intvs['Mod, Cov=20%'] = {'ACF': {'Type': 'mod', 'Scale': 0.2}}


def fn_post(y0, m, intv):
    p = y0['pars']
    _, ms, _ = m.simulate_onward(y0=y0['y0'], p=p, intv=intv)
    return ms


if __name__ == '__main__':
    from joblib import Parallel, delayed
    from sim.util import bind_results
    import pickle as pkl
    from sim import load_inputs
    from sim.dy import Model

    inputs = load_inputs('../data/pars.json')
    m = Model(inputs, year0=1970)

    for path in ds:
        out_path = f'../out/{path}'

        # Posterior run
        y0s = pkl.load(open(f'{out_path}/y0_national.pkl', 'rb'))

        mss = list()
        for scenario, intv in intvs.items():
            print(scenario)
            with Parallel(n_jobs=3, verbose=8) as parallel:
                mss0 = parallel(delayed(fn_post)(y0, m, intv) for y0 in y0s)

            mss.append(bind_results(mss0, keys=True, Scenario=scenario))

        mss = bind_results(mss, keys=False)
        mss.to_csv(f'{out_path}/Runs_Intv.csv')
