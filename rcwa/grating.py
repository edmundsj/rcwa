from rcwa import Layer, Slicer, Crystal
import numpy as np

class Grating(Layer):
    """
    Class used for simple generation of 1D gratings. By default oriented with periodicity along the x-direction.

    :param period: Spatial period of the grating
    :param t: Thickness of the grating along z
    :param er: Permittivity of the grating material
    :param ur: Permeability of the grating material
    :param n: Refractive index of the grating material. Overrides permittivity/permeability if used.
    :param er_void: Permittivity of the voids in the grating
    :param ur_void: Permeability of the voids in the grating
    :param n_void: Refractive index of the voids in the grating. Overrides permittivity / permeability iif used.
    :param groove_width: Width of the empty spaces (voids) in the grating. Must be smaller than period
    :param Nx: Number of points along x to divide the grating into
    :param shape: Shape of the grating. Current options are 'rectangular'
    """
    def __init__(self, period=1, er=2, ur=1, n=None, t=0.1,
                 er_void=1, ur_void=1, n_void=1,
                 groove_width=0.5, Nx=500, shape='rectangular'):
        self.period = period

        if groove_width > period:
            raise ValueError(f'Groove width {groove_width} must be larger than period {period}')

        groove_fraction = groove_width / period

        if shape == 'rectangular' or shape == 'rect':
            er_data, ur_data = self._groove_data(
                er=er, ur=ur, n=n, er_void=er_void, ur_void=ur_void, n_void=n_void,
                Nx=Nx, groove_fraction=groove_fraction)
        else:
            raise ValueError(f'shape {shape} is not available.')

        lattice_vector = np.array([period, 0])

        crystal = Crystal(er_data, ur_data, lattice_vector)
        super().__init__(L=t, crystal=crystal)

    def _groove_data(self, er=2, er_void=1, ur=1, ur_void=1, n=None, n_void=1, Nx=500,
                              groove_fraction=0.5):
        if n is not None:
            er_data = self._groove_data_single(np.square(n_void), np.square(n), Nx, groove_fraction)
            ur_data = np.ones(er_data.shape)
        else:
            er_data = self._groove_data_single(er_void, er, Nx, groove_fraction)
            ur_data = self._groove_data_single(ur_void, ur, Nx, groove_fraction)

        return er_data, ur_data

    def _groove_data_single(self, val1, val2, Nx, switch_fraction):
        positions = np.linspace(1/Nx, 1, Nx)
        void_positions = positions <= switch_fraction
        return (val1 - val2) * void_positions + val2

