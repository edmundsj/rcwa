from rcwa.shorthand import complexArray
import numpy as np
from numpy.typing import ArrayLike
from typing import Union


class Crystal:
    """
    Class used for defining periodic structures in x and y

    :param er: 2D numpy array of permittivity values
    :param ur: 2D numpy array of permeability values
    :param lattice_vectors: Real-space lattice vectors
    """
    def __init__(self, *lattice_vectors,
                 er: Union[ArrayLike, float] = 1,
                 ur: Union[ArrayLike, float] = 1):
        self.permeabilityCellData = ur
        self.permittivityCellData = er

        self.dimensions = len(lattice_vectors)
        raw_lattice_vectors = np.array(lattice_vectors)
        self.lattice_vectors = []

        if len(raw_lattice_vectors[0]) < self.dimensions:
            raise ValueError('Lattice vector does not have enough dimensions. Needs at least {self.dimensions}')
        if self.dimensions > 3:
            raise ValueError('Crystal number of dimensions too high ({self.dimensions}). Can only implement up to 3D.')

        if self.dimensions == 1:
            self.lattice_vectors.append(raw_lattice_vectors[0, 0:2])

        if self.dimensions == 2:
            self.lattice_vectors.append(raw_lattice_vectors[0, 0:2])
            self.lattice_vectors.append(raw_lattice_vectors[1, 0:2])

        if self.dimensions == 3:
            self.lattice_vectors.append(raw_lattice_vectors[0])
            self.lattice_vectors.append(raw_lattice_vectors[1])
            self.lattice_vectors.append(raw_lattice_vectors[2])

        if self.dimensions > 0:
            self.reciprocal_lattice_vectors = self.calculateReciprocalLatticeVectors()
            self.crystalType = self._crystal_type()
            self.latticeConstant = np.linalg.norm(self.lattice_vectors[0]) # TODO: Make this more general

    def calculateReciprocalLatticeVectors(self) -> ArrayLike:
        if self.dimensions == 1:
            return self._reciprocal_lattice_vectors_1d()
        elif self.dimensions == 2:
            return self._reciprocal_lattice_vectors_2d()
        elif self.dimensions == 3:
            return self._reciprocal_lattice_vectors_3d()

    def _reciprocal_lattice_vectors_1d(self) -> ArrayLike:
        t1 = self.lattice_vectors[0]
        t1_direction = t1 / np.linalg.norm(t1)

        T1 = 2 * np.pi/ np.linalg.norm(t1) * t1_direction
        return (T1,)

    def _reciprocal_lattice_vectors_2d(self) -> ArrayLike:
        rotationMatirx90Degrees = complexArray([
            [0, -1],
            [1, 0]])
        t1 = self.lattice_vectors[0]
        t2 = self.lattice_vectors[1]

        T1 = 2 * np.pi * rotationMatirx90Degrees @ t2 / np.dot(t1, rotationMatirx90Degrees @ t2)
        T2 = 2 * np.pi * rotationMatirx90Degrees @ t1 / np.dot(t2, rotationMatirx90Degrees @ t1)
        return (T1, T2)

    def _reciprocal_lattice_vectors_3d(self) -> ArrayLike:
        t1 = self.lattice_vectors[0]
        t2 = self.lattice_vectors[1]
        t3 = self.lattice_vectors[2]
        T1 = 2 * np.pi * np.cross(t2, t3) / np.dot(t1, np.cross(t2, t3))
        T2 = 2 * np.pi * np.cross(t3, t1) / np.dot(t2, np.cross(t3, t1))
        T3 = 2 * np.pi * np.cross(t1, t2) / np.dot(t3, np.cross(t1, t2))

        return (T1, T2, T3)

    def _crystal_type(self) -> str:
        if self.dimensions == 1:
            crystalType = 'SQUARE'
        elif self.dimensions == 2:
            crystalType = self._crystal_type_2d()
        elif self.dimensions == 3:
            crystalType = self._crystal_type_3d()

        return crystalType

    def _crystal_type_2d(self) -> str:
        epsilon = 0.00001
        sideLengthDifference = abs(np.linalg.norm(self.reciprocal_lattice_vectors[0]) -
                                   np.linalg.norm(self.reciprocal_lattice_vectors[1]))
        latticeVectorProjection = abs(np.dot(self.reciprocal_lattice_vectors[0], self.reciprocal_lattice_vectors[1]))

        if sideLengthDifference < epsilon and latticeVectorProjection < epsilon:
            return "SQUARE"
        elif sideLengthDifference > epsilon and latticeVectorProjection < epsilon:
            return "RECTANGULAR"
        elif latticeVectorProjection > epsilon:
            return "OBLIQUE"
        else:
            raise NotImplementedError

    def _crystal_type_3d(self) -> str:
        epsilon = 0.00001
        difference_1 =  abs(np.linalg.norm(self.reciprocal_lattice_vectors[0]) - np.linalg.norm(self.reciprocal_lattice_vectors[1]))
        difference_2 = abs(np.linalg.norm(self.reciprocal_lattice_vectors[0]) - np.linalg.norm(self.reciprocal_lattice_vectors[2]))
        max_difference = max(difference_1, difference_2)

        proj_01 = abs(np.dot(self.reciprocal_lattice_vectors[0], self.reciprocal_lattice_vectors[1]))
        proj_02 = abs(np.dot(self.reciprocal_lattice_vectors[0], self.reciprocal_lattice_vectors[2]))
        proj_03 = abs(np.dot(self.reciprocal_lattice_vectors[1], self.reciprocal_lattice_vectors[2]))
        max_proj = max(proj_01, proj_02, proj_03)

        if max_difference < epsilon and max_proj < epsilon:
            return "SQUARE"
        elif max_difference > epsilon and max_proj < epsilon:
            return "RECTANGULAR"
        elif max_proj > epsilon:
            return "OBLIQUE"
        else:
            raise NotImplementedError

"""
    def _key_symmetry_points_2d(self):
        keySymmetryPoints = []
        keySymmetryNames = []
        T1 = self.reciprocal_lattice_vectors[0]
        T2 = self.reciprocal_lattice_vectors[1]

        if self.crystalType == "SQUARE":
            keySymmetryNames = ["X", "G", "M"]
            keySymmetryPoints = [0.5 * T1, 0 * T1, 0.5 * (T1 + T2)]
        elif self.crystalType == "RECTANGULAR":
            keySymmetryNames = ["X", "G", "Y", "S"]
            keySymmetryPoints = [0.5 * T1, 0 * T1, 0.5 * T2, 0.5 * (T1 + T2)]
        else:
            raise NotImplementedError

        return (keySymmetryPoints, keySymmetryNames)
"""