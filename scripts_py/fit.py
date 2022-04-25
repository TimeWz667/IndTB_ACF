from sim import load_inputs, Objective
from sim.dy import Model, get_bn
from sim.util import simulate


inputs = load_inputs('../data/pars.json')
m = Model(inputs, year0=2010)
bn = get_bn('../prior')
to_fit = Objective(bn=bn,
                   model=m,
                   filepath_targets='../data/Targets.json')


def fn_post(p, m):
    ys, ms, msg = m.simulate(p)
    return ys, ms, msg['pars'] if 'pars' in msg else dict()


if __name__ == '__main__':
    from sims_pars.fitting.abcsmc import ApproxBayesComSMC
    import os
    from joblib import Parallel, delayed
    from sim.util import save_results
    import pandas as pd
    import pickle as pkl

    out_path = '../out/dy'

    smc = ApproxBayesComSMC(max_round=30, n_collect=150, n_core=4, verbose=8)
    smc.fit(to_fit)

    post = smc.Collector
    print(smc.Monitor.Trajectories)

    os.makedirs(out_path, exist_ok=True)

    smc.Monitor.save_trajectories(f'{out_path}/Trace.csv')
    post.save_to_csv(f'{out_path}/Post.csv')

    # Posterior run
    post = [dict(row) for _, row in pd.read_csv(f'{out_path}/Post.csv').iterrows()]

    with Parallel(n_jobs=4, verbose=8) as parallel:
        rss = parallel(delayed(fn_post)(pars, m) for pars in post)

    rss = [rs for rs in rss if rs[1] is not None]

    yss = [ys for ys, _, _ in rss]
    mss = [ms for _, ms, _ in rss]
    pss = [ps for _, _, ps in rss]

    save_results(mss, f'{out_path}/Runs_Post.csv')
    pd.DataFrame(pss).to_csv(f'{out_path}/Post_full.csv')

    res = [{'y0': list(ys.y[:, -1]), 'pars': pars} for pars, ys in zip(pss, yss)]

    with open(f'{out_path}/y0_national.pkl', 'wb') as f:
        pkl.dump(res, f)
