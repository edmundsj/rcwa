import context

import unittest
from material import Material
from shorthand import *
from shorthandTest import *

class TestMaterial(unittest.TestCase):

    def testParseCSV(self):
        wavelengthsToTest = np.array([0.25, 0.26, 0.27])
        nkDesired = complexArray([1.637 + 3.59j, 1.737 + 3.99j, 2.03 + 4.60j])
        nkCalculated = self.silicon.nkData[0:3]
        assertEqual(nkDesired, nkCalculated, errorMessage = "material: nk")

        erDesired = complexArray([-10.208331+11.75366j, -12.902931 + 13.86126j, -17.0391 + 18.676j])
        urDesired = complexArray([1, 1, 1])
        erCalculated = self.silicon.erData[0:3]
        urCalculated = self.silicon.urData[0:3]
        assertAlmostEqual(erDesired, erCalculated, errorMessage = "material: er")
        assertAlmostEqual(urDesired, urCalculated, errorMessage = "material: ur")

    def testnk(self):
        # Case 1: wavelength is exactly identical to one we have in database
        nDesired = 1.737 + 3.99j
        wavelength = 0.26
        nCalculated = self.silicon.n(wavelength)
        assertAlmostEqual(nCalculated, nDesired, errorMessage="material: testnk: n1")

        # Case 2: wavelength is nearly identical to one we have in database
        wavelength = 0.264
        nCalculated = self.silicon.n(wavelength)
        assertAlmostEqual(nCalculated, nDesired, errorMessage="material: testnk: n2")

        # Case 3: wavelength is larger than one we have in database - extrapolate linearly
        wavelength = 1.47
        nDesired = 3.485 + 1.09e-13j
        nCalculated = self.silicon.n(wavelength)
        assertAlmostEqual(nCalculated, nDesired, errorMessage="material: testnk: n3")

        # Case 4: wavelength is smaller than one we have in database
        wavelength = 0.23
        nDesired = 1.437 + 2.79j
        nCalculated = self.silicon.n(wavelength)
        assertEqual(nCalculated, nDesired, errorMessage="material: testnk: n4")

    def testEr(self):
        # Case 1: wavelength is exactly identical to one we have in database
        erDesired = sq(1.737 + 3.99j)
        wavelength = 0.26
        erCalculated = self.silicon.er(wavelength)
        assertAlmostEqual(erCalculated, erDesired, errorMessage="material: testnk: er1")

        # Case 2: wavelength is nearly identical to one we have in database
        wavelength = 0.264
        erCalculated = self.silicon.er(wavelength)
        assertAlmostEqual(erCalculated, erDesired, errorMessage="material: testnk: er2")

        # Case 3: wavelength is larger than one we have in database - extrapolate linearly
        wavelength = 1.47
        erDesired = sq(3.487 + 1.09e-13j) - 0.01395
        erCalculated = self.silicon.er(wavelength)
        assertAlmostEqual(erCalculated, erDesired, absoluteTolerance=1e-5, errorMessage="material: testnk: er3")

        # Case 4: wavelength is smaller than one we have in database
        wavelength = 0.23
        erDesired = sq(1.637 + 3.59j) + 5.3892 - 4.2152j
        erCalculated = self.silicon.er(wavelength)
        assertAlmostEqual(erCalculated, erDesired, absoluteTolerance = 1e-5,errorMessage="material: testnk: er4")

    def testUr(self):
        # Case 1: wavelength is exactly identical to one we have in database
        urDesired = 1
        wavelength = 0.26
        urDesired = 1
        urCalculated = self.silicon.ur(wavelength)
        assertAlmostEqual(urCalculated, urDesired, errorMessage="material: testnk: ur1")

        # Case 2: wavelength is nearly identical to one we have in database
        wavelength = 0.264
        urDesired = 1
        urCalculated = self.silicon.ur(wavelength)
        assertAlmostEqual(urCalculated, urDesired, errorMessage="material: testnk: ur2")

        # Case 3: wavelength is larger than one we have in database - extrapolate linearly
        wavelength = 1.47
        urDesired = 1
        urCalculated = self.silicon.ur(wavelength)
        assertAlmostEqual(urCalculated, urDesired, absoluteTolerance=1e-5, errorMessage="material: testnk: ur3")

        # Case 4: wavelength is smaller than one we have in database
        wavelength = 0.23
        urDesired = 1
        urCalculated = self.silicon.ur(wavelength)
        assertAlmostEqual(urCalculated, urDesired, absoluteTolerance = 1e-5,errorMessage="material: testnk: ur4")

    @classmethod
    def setUpClass(cls):
        cls.silicon = Material(str(context.nkLocation) + '/Si_Schinke.csv')
