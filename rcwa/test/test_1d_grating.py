import unittest
from rcwa.testing import *
from rcwa.shorthand import *
from rcwa import Source, Layer, LayerStack, Crystal
import numpy as np
from rcwa.matrices import *

class Test1DGrating(unittest.TestCase):
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
        grating_er = [1,1,1,3,3,3]
        grating_ur = [1,1,1,1,1,1]
        grating_thickness = 1
        N_harmonics_x = 2 # only x is used but include y for calrity
        N_harmocnis_y = 1


        reflectionLayer = Layer(erReflectionRegion, urReflectionRegion)
        self.source = Source(wavelength=self.wavelength, theta=self.theta, phi=self.phi,
                pTEM=[self.pTE, self.pTM], layer=reflectionLayer)
        transmissionLayer = Layer(erTransmissionRegion, urTransmissionRegion)
        c = Crystal([1,0],er=grating_er,ur=grating_ur)
        l = Layer(crystal=c, thickness=grating_thickness)
        self.layerStack = LayerStack(l, incident_layer=reflectionLayer, transmission_layer=transmissionLayer)

        self.Kx = np.diag(complexArray([7.096982989,0.813797681]))
        self.Ky = np.diag(complexArray([0.296198133,0.296198133]))

        self.KzReflectionRegion = \
             self.KzTransmissionRegion = \
             self.KzFreeSpace = np.diag(complexArray([0-7.03241785j,0.5-0j]))
        



if __name__ == '__main__':
    unittest.main()