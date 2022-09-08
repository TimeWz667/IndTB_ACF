import json
import numpy as np
import pandas as pd
from sim.dy import Model
from mc.pymc_util import DataModel

__all__ = ['load_objectives']


def load_objectives(path_target, exo=None):
    m = Model(year0=1970)
    targets = json.load(open(path_target, 'r'))
    targets = pd.DataFrame(targets[:-3])

    def simulator(p):
        _, ms, msg = m.simulate_to_fit(p)
        if not msg['succ']:
            return np.zeros(6)
        return np.array([ms[k] for k in ['Prev', 'ARTI', 'PrDR_CNR',
                                         'PrAsym', 'PrPreCS', 'PrExCS',
                                         # 'PrSp_Asym', 'PrSp_PreCS', 'PrSp_ExCS'
                                         ]])

    obs, eps = targets.M, (targets.U - targets.L) / 1.96
    dm = DataModel(np.array(obs), np.array(eps), simulator, exo=exo)
    return dm
