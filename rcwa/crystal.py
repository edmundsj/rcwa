from rcwa.shorthand import *

class Crystal:
    """
    Class used for defining periodic structures in x and y

    :param permittivityCellData: 2D numpy array of permittivity values
    :param permeabilityCellData: 2D numpy array of permeability values
    :param latticeVectors: Real-space lattice vectors
    """
    def __init__(self, permittivityCellData=1, permeabilityCellData=1, *latticeVectors):
        self.permeabilityCellData = permeabilityCellData
        self.permittivityCellData = permittivityCellData

        self.dimensions = len(latticeVectors)
        rawLatticeVectors = np.array(latticeVectors)
        self.latticeVectors = []
        if(self.dimensions == 2 and len(latticeVectors[0]) >= 2):
            self.latticeVectors.append(rawLatticeVectors[0, 0:2])
            self.latticeVectors.append(rawLatticeVectors[1, 0:2])

        if(self.dimensions > 0):
            self.reciprocalLatticeVectors = self.calculateReciprocalLatticeVectors();
            self.crystalType = self.determineCrystalType();
            (self.keySymmetryPoints, self.keySymmetryNames) = self.generateKeySymmetryPoints()
            self.latticeConstant = norm(self.latticeVectors[0]) # TODO: Make this more general

    def calculateReciprocalLatticeVectors(self):
        if self.dimensions == 2:
            return self.calculateReciprocalLatticeVectors2D();
        elif self.dimensions == 3:
            return self.calculateReciprocalLatticeVectors3D();
        else:
            raise NotImplementedError(f"Cannot calculate reciprocal lattice for {self.dimensions}D." +
                    " Not currently implemented.")

    def calculateReciprocalLatticeVectors2D(self):
        rotationMatirx90Degrees = complexArray([
            [0,-1],
            [1,0]]);
        t1 = self.latticeVectors[0];
        t2 = self.latticeVectors[1];

        T1 = 2 * pi * rotationMatirx90Degrees @ t2 / dot(t1, rotationMatirx90Degrees @ t2);
        T2 = 2 * pi * rotationMatirx90Degrees @ t1 / dot(t2, rotationMatirx90Degrees @ t1);
        return (T1, T2);

    def calculateReciprocalLatticeVectors3D(self):
        t1 = self.latticeVectors[0];
        t2 = self.latticeVectors[1];
        t3 = self.latticeVectors[2];
        T1 = 2 * pi * cross(t2, t3) / dot(t1, cross(t2, t3));
        T2 = 2 * pi * cross(t3, t1) / dot(t2, cross(t3, t1));
        T3 = 2 * pi * cross(t1, t2) / dot(t3, cross(t1, t2));

        return (T1, T2, T3);

    def determineCrystalType(self):
        if self.dimensions == 2:
            crystalType = self.determineCrystalType2D()
            return crystalType
        else:
            raise NotImplementedError

    def determineCrystalType2D(self):
        epsilon = 0.00001
        sideLengthDifference = abs(norm(self.reciprocalLatticeVectors[0]) -
                norm(self.reciprocalLatticeVectors[1]))
        latticeVectorProjection = abs(dot(self.reciprocalLatticeVectors[0], self.reciprocalLatticeVectors[1]))

        if sideLengthDifference < epsilon and latticeVectorProjection < epsilon:
            return "SQUARE"
        elif sideLengthDifference > epsilon and latticeVectorProjection < epsilon:
            return "RECTANGULAR"
        elif latticeVectorProjection > epsilon:
            return "OBLIQUE"
        else:
            raise NotImplementedError;

    def generateKeySymmetryPoints(self):
        keySymmetryPoints = []
        keySymmetryNames = []
        T1 = self.reciprocalLatticeVectors[0]
        T2 = self.reciprocalLatticeVectors[1]

        if self.crystalType == "SQUARE":
            keySymmetryNames = ["X", "G", "M"];
            keySymmetryPoints = [0.5 * T1, 0 * T1, 0.5 * (T1 + T2)];
        elif self.crystalType == "RECTANGULAR":
            keySymmetryNames = ["X", "G", "Y", "S"];
            keySymmetryPoints = [0.5 * T1, 0 * T1, 0.5 * T2, 0.5 * (T1 + T2)];
        else:
            raise NotImplementedError;

        return (keySymmetryPoints, keySymmetryNames);
