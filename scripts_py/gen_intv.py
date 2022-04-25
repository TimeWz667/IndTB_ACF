
intv_all = {
    'ImpDx': {'Dx': 0.90},
    'CS': {'Scale': 0.5},
}


intvs = {'(0)Baseline': dict()}
intv = dict()
for i, (k, v) in enumerate(intv_all.items(), 1):
    intv[k] = v
    if i == 0:
        intvs[f'({i})={k}'] = dict(intv)
    else:
        intvs[f'({i})=({i - 1}) + {k}'] = dict(intv)
    intvs[k] = {k: v}


def fn_post(y0, m, intv):
    _, ms, _ = m.simulate_onward(y0=y0['y0'], p=y0['pars'], intv=intv)
    return ms


if __name__ == '__main__':
    from joblib import Parallel, delayed
    from sim.util import bind_results
    import pickle as pkl
    from sim import load_inputs
    from sim.dy import Model

    inputs = load_inputs('../data/pars.json')
    m = Model(inputs, year0=2010)

    out_path = '../out/dy'

    # Posterior run
    y0s = pkl.load(open(f'{out_path}/y0_national.pkl', 'rb'))
    # y0s = y0s[:min(len(y0s), 200)]

    mss = list()
    for scenario, intv in intvs.items():
        print(scenario)
        with Parallel(n_jobs=3, verbose=8) as parallel:
            mss0 = parallel(delayed(fn_post)(y0, m, intv) for y0 in y0s)

        mss.append(bind_results(mss0, keys=True, Scenario=scenario))

    mss = bind_results(mss, keys=False)
    mss.to_csv(f'{out_path}/Runs_Intv.csv')
