import unittest

import numpy as np

from rcwa import test_dir
from rcwa.testing import *
from rcwa.matrices import *
from rcwa.shorthand import *
from rcwa.testing import assert_almost_equal
from rcwa import Source, Layer, LayerStack, Crystal

class Test1DGrating(unittest.TestCase):

    def testSetConvolutionMatrix(self):
        t1 = complexArray([1.75,0])
        erData = self.grating_er
        urData = self.grating_ur
        crystal = Crystal(t1, er=erData, ur=urData)
        dummyLayer = Layer(crystal=crystal)
        dummyLayer.set_convolution_matrices(self.numberHarmonics)

        convolutionMatrixActual = self.erConvolutionMatrixLayer
        convolutionMatrixCalculated = dummyLayer.er
        assert_almost_equal(convolutionMatrixActual, convolutionMatrixCalculated, self.absoluteTolerance,
                            self.relativeTolerance, "UR convolution matrices for layer not equal")

    def testPMatrix(self):
        PActual = self.PLayer
        PCalculated = calculatePMatrix(self.Kx, self.Ky, self.layerStack.internal_layers[0])
        assert_almost_equal(PActual, PCalculated, self.absoluteTolerance, self.relativeTolerance,
                "P matrix layer");

    def testQMatrix(self):
        QActual = self.QLayer
        QCalculated = calculateQMatrix(self.Kx, self.Ky, self.layerStack.internal_layers[0])
        assert_almost_equal(QActual, QCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Q matrix Layer");

        QActual = self.QReflectionRegion
        QCalculated = calculateQMatrix(self.Kx, self.Ky, self.layerStack.incident_layer)
        assert_almost_equal(QActual, QCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Q Reflection Region");

        QActual = self.QTransmissionRegion
        QCalculated = calculateQMatrix(self.Kx, self.Ky, self.layerStack.transmission_layer)
        assert_almost_equal(QActual, QCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Q Transmission Region");

    def testOmegaSquaredMatrix(self):
        OmegaSquaredActual = self.OmegaSquaredLayer
        OmegaSquaredCalculated = calculateOmegaSquaredMatrix(self.PLayer, self.QLayer)
        assert_almost_equal(OmegaSquaredActual, OmegaSquaredCalculated,
                            self.absoluteTolerance, self.relativeTolerance);

    @classmethod
    def setUpClass(self):
        deg = pi / 180
        self.absoluteTolerance = 1e-4
        self.relativeTolerance = 1e-3
        self.theta = 60 * deg
        self.phi = 20 * deg
        self.pTE = 1
        self.pTM = 0j
        self.wavelength = 0.5
        erReflectionRegion = 1.0
        urReflectionRegion = 1.0
        erTransmissionRegion = 1.0
        urTransmissionRegion = 1.0
        self.grating_er = np.asarray([1,1,1,3,3,3])
        self.grating_ur = np.asarray([1,1,1,1,1,1])
        grating_thickness = 1
        N_harmonics_x = 2 
        N_harmonics_y = 1
        self.numberHarmonics = N_harmonics_x


        reflectionLayer = Layer(erReflectionRegion, urReflectionRegion)
        self.source = Source(wavelength=self.wavelength, theta=self.theta, phi=self.phi,
                pTEM=[self.pTE, self.pTM], layer=reflectionLayer)
        transmissionLayer = Layer(erTransmissionRegion, urTransmissionRegion)
        c = Crystal([1,0],er=self.grating_er,ur=self.grating_ur)
        l = Layer(crystal=c, thickness=grating_thickness)
        self.layerStack = LayerStack(l, incident_layer=reflectionLayer, transmission_layer=transmissionLayer)
        self.erConvolutionMatrixLayer = numpyArrayFromFile(
            test_dir + "/1dGrating/layer/erConvolutionData.txt")
        urConvolutionMatrixLayer = complexIdentity(2)

        self.Kx = np.diag(complexArray([7.096982989,0.813797681]))
        self.Ky = np.diag(complexArray([0.296198133,0.296198133]))
        # As with other tests this is a bit of a hack, but that's good for test purposes.
        self.layerStack.internal_layers[0].er = self.erConvolutionMatrixLayer
        self.layerStack.internal_layers[0].ur = urConvolutionMatrixLayer

        self.KzReflectionRegion = \
             self.KzTransmissionRegion = \
             self.KzFreeSpace = np.diag(complexArray([0-7.03241785j,0.5-0j]))

        self.QFreeSpace = numpyArrayFromFile(test_dir +
        "/1dGrating/freeSpace/QFreeSpace.txt")
        self.WFreeSpace = complexIdentity(4)
        self.LambdaFreeSpace = numpyArrayFromFile(test_dir +
                "/1dGrating/freeSpace/LambdaFreeSpace.txt")

        self.VFreeSpace = numpyArrayFromFile(test_dir +
        "/1dGrating/freeSpace/VFreeSpace.txt")

        self.S11Transparent = complexZeros((4, 4))
        self.S22Transparent = complexZeros((4, 4))
        self.S21Transparent = complexIdentity(4)
        self.S12Transparent = complexIdentity(4)

        self.PLayer = numpyArrayFromFile(test_dir +
                "/1dGrating/layer/PLayer.txt")
        self.QLayer = numpyArrayFromFile(test_dir +
                "/1dGrating/layer/QLayer.txt")
        self.OmegaSquaredLayer = numpyArrayFromFile(test_dir +
                "/1dGrating/layer/OmegaSquaredLayer.txt")
        self.LambdaLayer= numpyArrayFromFile(test_dir +
                "/1dGrating/layer/LambdaLayer.txt")

        self.VLayer = numpyArrayFromFile(test_dir +
                "/1dGrating/layer/VLayer.txt")
        self.ALayer = numpyArrayFromFile(test_dir +
                "/1dGrating/layer/ALayer.txt")
        self.BLayer = numpyArrayFromFile(test_dir +
                "/1dGrating/layer/BLayer.txt")
        self.XLayer = numpyArrayFromFile(test_dir +
                "/1dGrating/layer/XLayer.txt")
        self.WLayer = numpyArrayFromFile(test_dir +
                "/1dGrating/layer/WLayer.txt")
        self.S11Layer = numpyArrayFromFile(test_dir +
                "/1dGrating/layer/S11Layer.txt")
        self.S12Layer = numpyArrayFromFile(test_dir +
                "/1dGrating/layer/S12Layer.txt")
        self.S21Layer = self.S12Layer
        self.S22Layer = self.S11Layer
        self.SLayer = complexArray([
            [self.S11Layer, self.S12Layer],
            [self.S21Layer, self.S22Layer]])

        # For transmission and reflection their u and e match free space so derive from those
        self.QReflectionRegion = self.QFreeSpace
        self.LambdaReflectionRegion = self.LambdaFreeSpace
        self.WReflectionRegion = self.WFreeSpace
        self.VReflectionRegion = self.QReflectionRegion@np.diag(pinv(sqrt(self.LambdaFreeSpace)))
        self.LambdaReflectionRegion = self.LambdaFreeSpace
        self.AReflectionRegion= pinv(self.WFreeSpace) @ self.WReflectionRegion + inv(self.VFreeSpace) @ self.VReflectionRegion
        self.BReflectionRegion = pinv(self.WFreeSpace) @ self.WReflectionRegion - inv(self.VFreeSpace) @ self.VReflectionRegion
        self.S11ReflectionRegion = -pinv(self.AReflectionRegion)@self.BReflectionRegion
        self.S12ReflectionRegion = 2*pinv(self.AReflectionRegion)
        self.S21ReflectionRegion = 0.5*(self.AReflectionRegion-self.BReflectionRegion@pinv(self.AReflectionRegion)@self.BReflectionRegion)
        self.S22ReflectionRegion = self.BReflectionRegion@pinv(self.AReflectionRegion)
        self.SReflectionRegion = complexArray([
            [self.S11ReflectionRegion, self.S12ReflectionRegion],
            [self.S21ReflectionRegion, self.S22ReflectionRegion]])

        self.QTransmissionRegion = self.QFreeSpace
        self.LambdaTransmissionRegion = self.LambdaFreeSpace
        self.WTransmissionRegion = self.WFreeSpace
        self.VTransmissionRegion = self.QTransmissionRegion@np.diag(pinv(sqrt(self.LambdaFreeSpace)))
        self.LambdaTransmissionRegion = self.LambdaFreeSpace
        self.ATransmissionRegion= pinv(self.WFreeSpace) @ self.WTransmissionRegion + inv(self.VFreeSpace) @ self.VTransmissionRegion
        self.BTransmissionRegion = pinv(self.WFreeSpace) @ self.WTransmissionRegion - inv(self.VFreeSpace) @ self.VTransmissionRegion
        self.S11TransmissionRegion = self.BTransmissionRegion@pinv(self.ATransmissionRegion)
        self.S12TransmissionRegion = 0.5*(self.ATransmissionRegion-self.BTransmissionRegion@pinv(self.ATransmissionRegion)@self.BTransmissionRegion)
        self.S21TransmissionRegion = 2*pinv(self.ATransmissionRegion)
        self.S22TransmissionRegion = -pinv(self.ATransmissionRegion)@self.BTransmissionRegion

        self.STransmissionRegion = complexArray([
            [self.S11TransmissionRegion, self.S12TransmissionRegion],
            [self.S21TransmissionRegion, self.S22TransmissionRegion]])
        


if __name__ == '__main__':
    unittest.main()