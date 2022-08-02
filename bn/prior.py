from sims_pars import bayes_net_from_script

__author__ = 'Chu-Chang Ku'
__all__ = ['get_bn']


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


def get_bn(dir_prior='../../prior'):
    return read_parameters(
        'dy',
        f'{dir_prior}/demography.txt',
        f'{dir_prior}/transmission.txt',
        f'{dir_prior}/progression.txt',
        f'{dir_prior}/healthcare.txt'
    )


if __name__ == '__main__':
    from sims_pars import sample

    bn = get_bn()
    pc = sample(bn)

    for k, v in pc.items():
        try:
            print(k, ': ', v.reshape(-1))
        except AttributeError:
            print(k, ': ', v)
