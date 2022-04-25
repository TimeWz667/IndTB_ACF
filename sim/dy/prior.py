from sim.misc import read_parameters

__author__ = 'Chu-Chang Ku'
__all__ = ['get_bn']


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
