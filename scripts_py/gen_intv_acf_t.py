import os
import numpy as np

ds = [d for d in os.listdir("../out") if d.startswith('dy_sc')]

for path in ds:
    print(ds)

intvs = dict()

intvs['Baseline'] = {}

for r_acf in np.linspace(0.5, 3, 6):
    intvs[f'High, Universal, R_ACF={r_acf:.2f}'] = {'ACFPlain': {'Type': 'high', 'R_ACF': r_acf, 'Focus': False}}
    intvs[f'High, Focus, R_ACF={r_acf:.2f}'] = {'ACFPlain': {'Type': 'high', 'R_ACF': r_acf, 'Focus': True}}
    #intvs[f'Mod, R_ACF={r_acf:.2f}'] = {'ACFPlain': {'Type': 'mod', 'R_ACF': r_acf}}


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
        mss.to_csv(f'{out_path}/Runs_ACF_t.csv')
