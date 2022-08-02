import pandas as pd
from scipy.integrate import solve_ivp
import numpy as np
from numba import njit

__author__ = 'Chu-Chang Ku'
__all__ = ['calc_dy', 'warmup', 'update', 'simulate', 'bind_results']


@njit
def calc_dy(y, frs, tos, rates):
    dy = np.zeros_like(y)

    for i in range(len(frs)):
        fr, to, rate = frs[i], tos[i], rates[i]

        n_tr = rate * y[fr]
        dy[fr] -= n_tr
        dy[to] += n_tr

    return dy


def warmup(model, pars, t_warmup, t_start, dfe=None):
    y0 = model.get_y0(pars).reshape(-1)

    ys0 = solve_ivp(model, [t_start - t_warmup, t_start], y0, args=(pars,), events=dfe,
                    dense_output=True, method='RK23')

    if len(ys0.t_events[0]) > 0 or not ys0.success:
        return None, None, {'succ': False, 'res': 'DFE reached'}
    else:
        return ys0, None, {'succ': True}


def update(model, ys0, pars, t_out, dfe=None, intv=None):
    t_start, t_end = min(t_out), max(t_out)
    ys0 = np.array(ys0)
    ys = solve_ivp(model, [t_start, t_end], ys0, args=(pars, intv), events=dfe, dense_output=True)

    if len(ys.t_events[0]) > 0 or not ys.success:
        return None, None, {'succ': False, 'res': 'DFE reached'}

    ms = [model.measure(t, ys.sol(t), pars) for t in t_out if t >= t_start]
    ms = pd.DataFrame(ms).set_index('Time')

    msg = {'succ': True, 't_out': t_out, 'pars': pars}
    return ys, ms, msg


def simulate(model, pars, t_warmup, t_out, dfe=None):
    t_out = t_out[t_out > model.Year0]
    ys0, _, msg = warmup(model, pars, t_warmup, min(t_out), dfe=dfe)

    if not msg['succ']:
        return None, None, msg

    ys0 = ys0.y.T[-1]
    ys, ms, msg = update(model, ys0, pars, t_out=t_out, dfe=dfe)
    return ys, ms, msg


def bind_results(mss, keys=True, **kwargs):
    if bool(kwargs) and keys:
        mss = [ms.assign(Key=i, **kwargs) for i, ms in enumerate(mss) if ms is not None]
    elif bool(kwargs):
        mss = [ms.assign(**kwargs) for ms in mss if ms is not None]
    elif keys:
        mss = [ms.assign(Key=i) for i, ms in enumerate(mss) if ms is not None]
    else:
        mss = [ms for ms in mss if ms is not None]
    return pd.concat(mss)

