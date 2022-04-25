import json

__author__ = 'Chu-Chang Ku'
__all__ = ['load_inputs']


def load_inputs(src):
    with open(src, 'r') as f:
        src = json.load(f)
    return src


if __name__ == '__main__':
    inp0 = load_inputs('../data/pars.json')
    print(inp0)
