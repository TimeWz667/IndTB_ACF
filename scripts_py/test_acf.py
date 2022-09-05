import json
import numpy as np
import pandas as pd
from sim.dy import Model
from sim import load_inputs


path_data = '../data/pars.json'

inputs = load_inputs(path_data)
m = Model(inputs, year0=1970)


with open('../out/dy_free/Post.json', 'r') as f:
    pars_post = json.load(f)


if __name__ == '__main__':
    p0 = pars_post[0]

    y0, ms, _ = m.simulate_to_baseline(p0, t_start=2023)
    print(y0)
    print(ms)
