from rcwa.shorthand import *


class Crystal:
    """
    Class used for defining periodic structures in x and y

    :param er: 2D numpy array of permittivity values
    :param ur: 2D numpy array of permeability values
    :param lattice_vectors: Real-space lattice vectors
    """
    def __init__(self, *lattice_vectors, er=1, ur=1):
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
            self.reciprocal_lattice_vectors = self.calculateReciprocalLatticeVectors();
            self.crystalType = self._crystal_type();
            self.latticeConstant = norm(self.lattice_vectors[0]) # TODO: Make this more general

    def calculateReciprocalLatticeVectors(self):
        if self.dimensions == 1:
            return self._reciprocal_lattice_vectors_1d()
        elif self.dimensions == 2:
            return self._reciprocal_lattice_vectors_2d()
        elif self.dimensions == 3:
            return self._reciprocal_lattice_vectors_3d()

    def _reciprocal_lattice_vectors_1d(self):
        t1 = self.lattice_vectors[0]
        t1_direction = t1 / np.linalg.norm(t1)

        T1 = 2 * pi / np.linalg.norm(t1) * t1_direction
        return (T1,)

    def _reciprocal_lattice_vectors_2d(self):
        rotationMatirx90Degrees = complexArray([
            [0, -1],
            [1, 0]])
        t1 = self.lattice_vectors[0]
        t2 = self.lattice_vectors[1]

        T1 = 2 * pi * rotationMatirx90Degrees @ t2 / dot(t1, rotationMatirx90Degrees @ t2);
        T2 = 2 * pi * rotationMatirx90Degrees @ t1 / dot(t2, rotationMatirx90Degrees @ t1);
        return (T1, T2)

    def _reciprocal_lattice_vectors_3d(self):
        t1 = self.lattice_vectors[0]
        t2 = self.lattice_vectors[1]
        t3 = self.lattice_vectors[2]
        T1 = 2 * pi * cross(t2, t3) / dot(t1, cross(t2, t3))
        T2 = 2 * pi * cross(t3, t1) / dot(t2, cross(t3, t1))
        T3 = 2 * pi * cross(t1, t2) / dot(t3, cross(t1, t2))

        return (T1, T2, T3)

    def _crystal_type(self):
        if self.dimensions == 1:
            crystalType = 'SQUARE'
        elif self.dimensions == 2:
            crystalType = self._crystal_type_2d()
        elif self.dimensions == 3:
            crystalType = self._crystal_type_3d()

        return crystalType

    def _crystal_type_2d(self):
        epsilon = 0.00001
        sideLengthDifference = abs(norm(self.reciprocal_lattice_vectors[0]) -
                                   norm(self.reciprocal_lattice_vectors[1]))
        latticeVectorProjection = abs(dot(self.reciprocal_lattice_vectors[0], self.reciprocal_lattice_vectors[1]))

        if sideLengthDifference < epsilon and latticeVectorProjection < epsilon:
            return "SQUARE"
        elif sideLengthDifference > epsilon and latticeVectorProjection < epsilon:
            return "RECTANGULAR"
        elif latticeVectorProjection > epsilon:
            return "OBLIQUE"
        else:
            raise NotImplementedError;

    def _crystal_type_3d(self):
        epsilon = 0.00001
        difference_1 =  abs(norm(self.reciprocal_lattice_vectors[0]) - norm(self.reciprocal_lattice_vectors[1]))
        difference_2 = abs(norm(self.reciprocal_lattice_vectors[0]) - norm(self.reciprocal_lattice_vectors[2]))
        max_difference = max(difference_1, difference_2)

        proj_01 = abs(dot(self.reciprocal_lattice_vectors[0], self.reciprocal_lattice_vectors[1]))
        proj_02 = abs(dot(self.reciprocal_lattice_vectors[0], self.reciprocal_lattice_vectors[2]))
        proj_03 = abs(dot(self.reciprocal_lattice_vectors[1], self.reciprocal_lattice_vectors[2]))
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
            keySymmetryNames = ["X", "G", "M"];
            keySymmetryPoints = [0.5 * T1, 0 * T1, 0.5 * (T1 + T2)];
        elif self.crystalType == "RECTANGULAR":
            keySymmetryNames = ["X", "G", "Y", "S"];
            keySymmetryPoints = [0.5 * T1, 0 * T1, 0.5 * T2, 0.5 * (T1 + T2)];
        else:
            raise NotImplementedError;

        return (keySymmetryPoints, keySymmetryNames);
"""