from sim.inputs import load_inputs
from sim.dy import Model
from bn.prior import get_bn
from bn.obj import Objective


inputs = load_inputs('../data/pars.json')
m = Model(inputs, year0=1970)
bn = get_bn('../bn/prior')


def fn_post(p, m):
    ys, ms, msg = m.simulate(p)
    return ys, ms, msg['pars'] if 'pars' in msg else dict()


if __name__ == '__main__':
    from sims_pars.fitting import ApproxBayesComSMC
    import os
    from joblib import Parallel, delayed
    from sim.util import bind_results
    import pandas as pd
    import pickle as pkl

    smc = ApproxBayesComSMC(max_round=40, n_collect=1000, n_core=5, verbose=8)

    out_path = f'../out/dy'
    os.makedirs(out_path, exist_ok=True)

    to_fit = Objective(bn=bn, model=m, filepath_targets='../data/Targets.json')

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

    mss = bind_results(mss)
    mss.to_csv(f'{out_path}/Runs_Post.csv')

    pd.DataFrame(pss).to_csv(f'{out_path}/Post_full.csv')
