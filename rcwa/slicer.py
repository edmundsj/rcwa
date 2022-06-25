import numpy as np

class Slicer:
    """
    Helper class that slices arbitrary functions of x, y, and z into several layers.

    :param func: Function f(x, y, z) to slice
    :param Nx: Number of desired evaluation points along x.
    :param Ny: Number of desired evaluation points along y.
    :param Nz: Number of desired slices along z.
    :param data: Raw 3-dimensional material property data

    :param xmin: Minimum x-value
    :param xmax: Maximum x value to be evaluated
    :param ymin: Minimum y-value to be evaluated
    :param ymax: Maximum y-value to be evaluated
    :param zmxn: Minimum z-value to be evaluated
    :param zmxn: Maximum z-value to be evaluated
    """
    def __init__(self, func=None, data=None, Nx=500, Ny=500, Nz=10,
                 xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, zmax=1):
        if func is None and data is None:
            raise ValueError('Must pass in a value for either func or data')

        self.func = func
        self.data = data

        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz

        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax

    def coordinates(self):
        x_coors = np.linspace(self.xmin, self.xmax, self.Nx)
        y_coors = np.linspace(self.ymin, self.ymax, self.Ny)
        z_coors = np.linspace(self.zmin, self.zmax, self.Nz)
        x, y, z = np.meshgrid(x_coors, y_coors, z_coors, indexing='ij')
        return x, y, z

    def slice(self):
        if self.func is not None:
            return self._slice_func()
        elif self.data is not None:
            return self._slice_data()

    def _slice_func(self):
        x, y, z = self.coordinates()
        return self.func(x, y, z)

    def _slice_data(self):
        raise NotImplementedError

