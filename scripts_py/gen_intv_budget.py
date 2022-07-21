import numpy as np


intvs = dict()

intvs['Baseline'] = {}

for budget in np.concatenate([np.linspace(0, 0.002, 6)[:-1], np.linspace(0.002, 0.05, 13)]):
    intvs[f'High, Budget={budget}, Risk'] = {'ACF': {'Type': 'high', 'Yield': budget, 'HiRisk': True}}
    intvs[f'High, Budget={budget}, Universal'] = {'ACF': {'Type': 'high', 'Yield': budget, 'HiRisk': False}}
        # intvs[f'Mod, Budget={budget}, Risk'] = {'ACF': {'Type': 'mod', 'Yield': budget, 'HiRisk': True}}
        # intvs[f'Mod, Budget={budget}, Universal'] = {'ACF': {'Type': 'mod', 'Yield': budget, 'HiRisk': False}}


def fn_post(y0, m, intv):
    p = y0['pars']
    _, ms, _ = m.simulate_onward(y0=np.array(y0['y0']), p=p, intv=intv)
    return ms


if __name__ == '__main__':
    from joblib import Parallel, delayed
    from sim.util import bind_results
    import pickle as pkl
    from sim import load_inputs
    from sim.dy import Model

    inputs = load_inputs('../data/pars.json')
    m = Model(inputs, year0=1970)

    out_path = '../out/dy'

    # Posterior run
    y0s = pkl.load(open(f'{out_path}/y0_national.pkl', 'rb'))

    mss = list()
    for scenario, intv in intvs.items():
        print(scenario)
        with Parallel(n_jobs=3, verbose=8) as parallel:
            mss0 = parallel(delayed(fn_post)(y0, m, intv) for y0 in y0s)

        mss.append(bind_results(mss0, keys=True, Scenario=scenario))

    mss = bind_results(mss, keys=False)
    mss.to_csv(f'{out_path}/Runs_IntvE.csv')
