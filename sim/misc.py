# from sims_pars import bayes_net_from_script

__author__ = 'Chu-Chang Ku'
__all__ = ['read_parameters']


def read_parameters(name, root, *filepaths):
    with open(root, mode='r') as f:
        script = f.read()
        bn = bayes_net_from_script(script)

    filepaths = list(filepaths)

    if len(filepaths) > 0:
        for i, filepath in enumerate(filepaths):
            with open(filepath, mode='r') as f:
                script = f.read()
                bn = bn.merge(f'{name}{i}', bayes_net_from_script(script))
    bn.Name = name
    return bn


if __name__ == '__main__':
    from sims_pars import sample

    bn = read_parameters('dy',
                         '../prior/demography.txt',
                         '../prior/transmission.txt',
                         '../prior/healthcare.txt',
                         '../prior/progression.txt')

    for k, v in sample(bn).items():
        print(k, v)
