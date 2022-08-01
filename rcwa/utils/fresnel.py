# Fresnel equations definitions for TE and TM modes
from rcwa.matrices import *
import numpy as np
from rcwa.shorthand import sqrt

def k_vector(source, layer, normalize=False):
    kx = layer.n * sin(source.theta) * cos(source.phi)
    ky = layer.n * sin(source.theta) * sin(source.phi)
    kz = layer.n * cos(source.theta)
    k_vec = complexArray([kx, ky, kz])
    if not normalize:
        k_vec *= 2 * np.pi / source.wavelength
    return k_vec

def rTE(source, layer1, layer2):
    k_inc = k_vector(source, layer1)
    k_trans = k_vector(source, layer2)
    kz1 = k_inc[2]
    kz2 = k_trans[2]
    ur1 = layer1.ur
    ur2 = layer2.ur

    return (ur2 * kz1 - ur1 * kz2) / (ur2 * kz1 + ur1 * kz2)


def tTE(source, layer1, layer2):
    return 1 + rTE(source, layer1, layer2);


def rTM(source, layer1, layer2):
    k_inc = k_vector(source, layer1)
    k_trans = k_vector(source, layer2)
    kz1 = k_inc[2]
    kz2 = k_trans[2]
    er1 = layer1.er
    er2 = layer2.er

    return (er1 * kz2 - er2 * kz1) / (er1 * kz2 + er2 * kz1);

# TODO: TEST MORE VIGOROUSLY
def tTM(source, layer1, layer2):
    er1 = layer1.er
    er2 = layer2.er
    ur1 = layer1.ur
    ur2 = layer2.ur

    eta_1 = sqrt(ur1 / er1);
    eta_2 = sqrt(ur2 / er2);
    return eta_2 / eta_1 * (1 + rTM(source, layer1, layer2));