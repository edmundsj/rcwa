from rcwa.shorthand import *
from rcwa import Material

# TODO: Convolution matrix generation must be refactored. It's a hot mess and hard to understand.


class Layer:
    """
    Class for defining a single layer of a layer stack used in a simulation

    :param er: Permittivity of the layer. Overridden by crystal permittivity if specified.
    :param ur: Permeability of the layer. Overridden by crystal permeability if specified.
    :param thickness: Thickness of the layer
    :param n: Refractive index of the layer. Overridden by cristal er/ur if specified.
    :param material: Material object containing the material's permittivity and permeability as a function of wavelength/angle.
    :param crystal: Crystal object if the layer is periodic in x and/or y. Overrides er, ur, n, and material
    """
    def __init__(self, er=1.0, ur=1.0, thickness=0.0, n=None, material=None, crystal=None):
        if material is None:
            self.material = Material(er=er, ur=ur, n=n)
        else:
            self.material = material

        self.thickness = thickness
        self.crystal = crystal

        if crystal is not None:
            self.homogenous = False
            #if numberHarmonics is not None:

                #self.setConvolutionMatrix(numberHarmonics)
        else:
            self.homogenous = True

    # Note: these are all just transparent wrappers for underlying material

    @property
    def er(self):
        return self.material.er

    @er.setter
    def er(self, er):
        self.material.er = er

    @property
    def ur(self):
        return self.material.ur

    @ur.setter
    def ur(self, ur):
        self.material.ur = ur

    @property
    def n(self):
        return self.material.n

    @n.setter
    def n(self, n):
        self.material.n = n

    @property
    def source(self):
        return self.material.source

    @source.setter
    def source(self, source):
        self.material.source = source

    def set_convolution_matrices(self, numberHarmonics):
        if self.crystal is not None:
            self.er = self._convolution_matrix(self.crystal.permittivityCellData, numberHarmonics)
            self.ur = self._convolution_matrix(self.crystal.permeabilityCellData, numberHarmonics)
        else:
            self.er = self.er * complexIdentity(prod(numberHarmonics))
            self.ur = self.ur * complexIdentity(prod(numberHarmonics))

    def _convolution_matrix(self, cellData, n_harmonics):
        dimension = self.crystal.dimensions;

        if isinstance(n_harmonics, int):
            n_harmonics = (n_harmonics,)

        if dimension == 1:
            n_harmonics = (n_harmonics + (1, 1))
        elif dimension == 2:
            n_harmonics = (n_harmonics + (1,))

        (P, Q, R) = n_harmonics

        convolutionMatrixSize = P*Q*R;
        convolutionMatrixShape = (convolutionMatrixSize, convolutionMatrixSize);
        convolutionMatrix = complexZeros(convolutionMatrixShape)

        cellData = reshapeLowDimensionalData(cellData);
        (Nx, Ny, Nz) = cellData.shape;
        zeroHarmonicsLocation = np.array([math.floor(Nx/2), math.floor(Ny/2), math.floor(Nz/2)])

        cellFourierRepresentation = fftn(cellData);
        for rrow in range(R):
            for qrow in range(Q):
                for prow in range(P):
                    row = rrow*Q*P + qrow*P + prow;
                    for rcol in range(R):
                        for qcol in range(Q):
                            for pcol in range(P):
                                col = rcol*Q*P + qcol*P + pcol;
                                # Get the desired harmonics relative to the 0th-order harmonic.
                                desiredHarmonics = np.array([prow - pcol, qrow - qcol, rrow - rcol])

                                # Get those harmonic locations from the zero harmonic location.
                                desiredHarmonicsLocation = zeroHarmonicsLocation + desiredHarmonics

                                convolutionMatrix[row][col] = \
                                    cellFourierRepresentation[desiredHarmonicsLocation[0],
                                    desiredHarmonicsLocation[1], desiredHarmonicsLocation[2]];
        if convolutionMatrix.shape == (1, 1):
            convolutionMatrix = convolutionMatrix[0][0]
        return convolutionMatrix;

    def __eq__(self, other):
        if not isinstance(other, Layer):
            return NotImplemented

        return self.er == other.er and self.ur == other.ur and self.thickness == other.thickness \
               and self.n == other.n and self.crystal == other.crystal

    def __str__(self):
        return f'Layer with\n\ter: {self.er}\n\tur: {self.ur}\n\tL: {self.thickness}\n\tn: {self.n}\n\tcrystal: {self.crystal}'

    def __repr__(self):
        return f'Layer with\n\ter: {self.er}\n\tur: {self.ur}\n\tL: {self.thickness}\n\tn: {self.n}\n\tcrystal: {self.crystal}'


freeSpaceLayer = Layer(1,1)


class LayerStack:
    """
    Class that defines overall geometry in terms of a stack of layers

    :param internal_layers: Layer objects, starting with the top-most layer (reflection region) and ending with the top-most region (substrate)
    :param incident_layer: Semi-infinite layer of incident region. Defaults to free space
    :param transmission_layer: Semi-infinite layer of transmission region. Defaults to free space
    """
    def __init__(self, *internal_layers, incident_layer=freeSpaceLayer, transmission_layer=freeSpaceLayer):
        if len(internal_layers) == 1:
            if isinstance(internal_layers[0], list):
                internal_layers = internal_layers[0]
        self.gapLayer = Layer(er=1,ur=1)
        self.incident_layer = incident_layer
        self.transmission_layer = transmission_layer

        self.internal_layers = list(internal_layers)

    def __eq__(self, other):
        if not isinstance(other, LayerStack):
            return NotImplemented

        reflection_layers_same = self.incident_layer == other.incident_layer
        transmission_layers_same = self.transmission_layer == other.transmission_layer
        internal_layers_same = False
        if len(self.internal_layers) == len(other.internal_layers):
            for i in range(len(self.internal_layers)):
                if self.internal_layers[i] != other.internal_layers[i]:
                    break
            internal_layers_same = True

        return internal_layers_same and reflection_layers_same and transmission_layers_same

    def __str__(self):
        top_string = f'\nReflection Layer:\n\t' + str(self.incident_layer) + \
                f'\nTransmissionLayer:\n\t' + str(self.transmission_layer) + \
                f'\nInternal Layer Count: {len(self.internal_layers)}\n'
        internal_string = ''
        for layer in self.internal_layers:
            internal_string += str(layer) + '\n'
        return top_string + internal_string

    def __repr__(self):
        return self.__str__()

    @property
    def source(self):
        return self._source

    @source.setter
    def source(self, source):
        self._source = source
        for layer in self.internal_layers:
            layer.source = self.source
        self.incident_layer.source = source
        self.transmission_layer.source = source

    def _set_gap_layer(self, kx, ky):
        self.gapLayer.er = 1 + sq(kx) + sq(ky)
        self.gapLayer.ur = 1

    # set all convolution matrices for all interior layers
    def set_convolution_matrices(self, numberHarmonics):
        for layer in self.internal_layers:
            layer.set_convolution_matrices(numberHarmonics)

    @property
    def crystal(self):
        for i in range(len(self.internal_layers)):
            if self.internal_layers[i].crystal is not None:
                return self.internal_layers[i].crystal
        return None


emptyStack = LayerStack()
