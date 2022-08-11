from sim import load_inputs
from sim.dy import Model


inputs = load_inputs('../data/pars.json')
m = Model(inputs, year0=1970)


def fn_post(p, m: Model):
    ys, ms, _ = m.simulate_to_baseline(p, 2022)
    return ys, ms


if __name__ == '__main__':
    from joblib import Parallel, delayed
    import pandas as pd
    import pickle as pkl
    import os

    ds = os.listdir('../out')
    ds = [d for d in ds if d.startswith('dy_')]
    print(ds)

    for d in ds:
        print(d)
        out_path = f'../out/{d}'

        # Posterior run
        post = [dict(row) for _, row in pd.read_csv(f'{out_path}/Post.csv').iterrows()]

        with Parallel(n_jobs=4, verbose=8) as parallel:
            rss = parallel(delayed(fn_post)(pars, m) for pars in post)

        rss = [rs for rs in rss if rs[1] is not None]

        yss = [ys for ys, _ in rss]
        mss = [ms for _, ms in rss]

        res = [{'y0': ys, 'pars': pars} for pars, ys in zip(post, yss)]

        with open(f'{out_path}/y0s.pkl', 'wb') as f:
            pkl.dump(res, f)
