from rcwa import Layer, Slicer, Crystal
import numpy as np


class Grating(Layer):
    """
    Base class that doesn't do much at all.
    """

    def _set_eun(self, n, n_void, er, er_void, ur, ur_void):
        if n is not None:
            self._er = np.square(n)
            self._er_void = np.square(n_void)
            self._ur = 1
            self._ur_void = 1
            self._n = n
            self._n_void = n_void
        else:
            self._er = er
            self._er_void = er_void
            self._ur = ur
            self._ur_void = ur_void
            self._n = np.sqrt(er*ur)
            self._n_void = np.sqrt(er_void*ur_void)

    def set_lv_period(self, period, lattice_vector):
        if lattice_vector is not None:
            self.period = np.linalg.norm(lattice_vector)
            self.lattice_vector = lattice_vector
        else:
            self.lattice_vector = np.array([period, 0])
            self.period = period


class TriangularGrating(Grating):
    """
    Class for one-dimensional triangular (ramp or blaze) gratings

    :param period: Spatial period of the grating
    :param thickness: Maximum thickness of grating along z (minimum thickness is zero)
    :param er: Permittivity of the grating material
    :param ur: Permeability of the grating material
    :param n: Refractive index of the grating material. Overrides permittivity/permeability if used.
    :param er_void: Permittivity of the voids in the grating
    :param ur_void: Permeability of the voids in the grating
    :param n_void: Refractive index of the voids in the grating. Overrides permittivity / permeability iif used.
    :param Nx: Number of points along x to divide the grating into
    :param Nx: Number of slices along z to divide the grating into
    :param lattice_vector: Explicit lattice vector for grating. Overrides period.
    """
    def __init__(self, period=1, er=2, ur=1, n=None, thickness=0.1,
                 er_void=1, ur_void=1, n_void=1, Nx=500, Nz=10,
                 lattice_vector=None):
        self.Nx = Nx
        self.Nz = Nz
        self.thickness = thickness

        self._set_eun(n=n, n_void=n_void, er=er, er_void=er_void, ur=ur, ur_void=ur_void)
        self.set_lv_period(period=period, lattice_vector=lattice_vector)
        super().__init__(thickness=thickness)

    def slice(self):
        er_slices, ur_slices = self._er_data()
        crystals = [Crystal(self.lattice_vector, er=er, ur=ur) \
                    for er, ur in zip(er_slices, ur_slices)]
        self.layers = [Layer(crystal=crystal, thickness=self.thickness / self.Nz) for crystal in crystals]

        return self.layers

    def _er_data(self):
        def triangle_func_er(x, y, z):
            in_void = z >= self.thickness * x / self.period
            er = in_void * (self._er_void - self._er) + self._er
            return er

        def triangle_func_ur(x, y, z):
            in_void = z >= self.thickness * x / self.period
            ur = in_void * (self._ur_void - self._ur) + self._ur
            return ur

        slicer_er = Slicer(
            func=triangle_func_er, xmin=0, xmax=self.period, ymin=0, ymax=1, zmin=0, zmax=self.thickness,
            Nx=self.Nx, Ny=1, Nz=self.Nz)
        slicer_ur = Slicer(
            func=triangle_func_ur, xmin=0, xmax=self.period, ymin=0, ymax=1, zmin=0, zmax=self.thickness,
            Nx=self.Nx, Ny=1, Nz=self.Nz)

        er_data = slicer_er.slice()
        ur_data = slicer_ur.slice()

        er_slices = [er_data[:,0, i] for i in range(er_data.shape[2])]
        ur_slices = [ur_data[:,0,i] for i in range(ur_data.shape[2])]
        er_slices.reverse()
        ur_slices.reverse()
        return er_slices, ur_slices


class RectangularGrating(Grating):
    """
    Class used for simple generation of 1D gratings. By default oriented with periodicity along the x-direction.

    :param period: Spatial period of the grating
    :param thickness: Thickness of the grating along z
    :param er: Permittivity of the grating material
    :param ur: Permeability of the grating material
    :param n: Refractive index of the grating material. Overrides permittivity/permeability if used.
    :param er_void: Permittivity of the voids in the grating
    :param ur_void: Permeability of the voids in the grating
    :param n_void: Refractive index of the voids in the grating. Overrides permittivity / permeability iif used.
    :param groove_width: Width of the empty spaces (voids) in the grating. Must be smaller than period
    :param nx: Number of points along x to divide the grating into
    :param lattice_vector: Explicit lattice vector for grating. Overrides period.
    """
    def __init__(self, period=1, er=2, ur=1, n=None, thickness=0.1,
                 er_void=1, ur_void=1, n_void=1,
                 groove_width=0.5, nx=500, lattice_vector=None):

        if groove_width > period:
            raise ValueError(f'Groove width {groove_width} must be larger than period {period}')

        self.thickness = thickness
        self.nx = nx

        self._set_eun(n=n, n_void=n_void, er=er, er_void=er_void, ur=ur, ur_void=ur_void)
        self.set_lv_period(period=period, lattice_vector=lattice_vector)

        groove_fraction = groove_width / period

        er_data, ur_data = self._er_data(
            er=er, ur=ur, n=n, er_void=er_void, ur_void=ur_void, n_void=n_void,
            Nx=nx, groove_fraction=groove_fraction)

        crystal = Crystal(self.lattice_vector, er=er_data, ur=ur_data)
        super().__init__(thickness=thickness, crystal=crystal)

    def _er_data(self, er=2, er_void=1, ur=1, ur_void=1, n=None, n_void=1, Nx=500,
                              groove_fraction=0.5):
        if n is not None:
            er_data = self._er_data_single(np.square(n_void), np.square(n), Nx, groove_fraction)
            ur_data = np.ones(er_data.shape)
        else:
            er_data = self._er_data_single(er_void, er, Nx, groove_fraction)
            ur_data = self._er_data_single(ur_void, ur, Nx, groove_fraction)

        return er_data, ur_data

    def _er_data_single(self, val1, val2, Nx, switch_fraction):
        positions = np.linspace(1/Nx, 1, Nx)
        void_positions = positions <= switch_fraction
        return (val1 - val2) * void_positions + val2
