from sim import load_inputs
from sim.dy import Model


inputs = load_inputs('../data/pars.json')
m = Model(inputs, year0=1970)


def fn_post(p, m):
    ys, ms, msg = m.simulate(p)
    if ms is not None:
        print(ms.PrDR_Inc[2020])
    return ys, ms, msg['pars'] if 'pars' in msg else dict()


if __name__ == '__main__':
    from joblib import Parallel, delayed
    from sim.util import save_results
    import pandas as pd
    import pickle as pkl

    out_path = '../out/DR'

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

    # with open(f'{out_path}/y0_national.pkl', 'wb') as f:
    #     pkl.dump(res, f)
