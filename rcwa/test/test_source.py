"""
Module for testing the plane wave source class
"""
from shorthand import *
from rcwa import Layer, Source
import unittest
from shorthandTest import *


class testSource(unittest.TestCase):
    def testKIncident(self):
        kIncidentActual = complexArray([1.0607, 0.61237, 0.70711])
        kIncidentCalculated = self.source.kIncident
        assertAlmostEqual(kIncidentActual, kIncidentCalculated, self.absoluteTolerance, self.relativeTolerance,
                "kIncident in testSource")

    def setUp(self):
        self.absoluteTolerance = 1e-5
        self.relativeTolerance = 1e-4
        wavelength = 0.02
        theta = 60*deg
        phi = 30*deg
        pTEM = 1 / sqrt(2) * complexArray([1,1j])
        reflectionLayer = Layer(er=2,ur=1)
        self.source = Source(wavelength, theta, phi, pTEM, reflectionLayer)
