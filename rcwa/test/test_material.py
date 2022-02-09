import unittest
from rcwa import Material, Source
from rcwa.testing import *
from rcwa.shorthand import *

class TestMaterial(unittest.TestCase):
    def testExtractMaterialDatabase(self):
        fake_material = Material()
        fake_material.extractMaterialDatabase()
        materials_to_check = ['Pt', 'Si', 'Ag', 'Ti', 'Au', 'SiO2']
        assert all(i in fake_material.materials.keys() for i in materials_to_check) == True

    @unittest.skip("deprecated and removed function")
    def testParseCSV(self):
        wavelengthsToTest = np.array([0.25, 0.26, 0.27])
        nkDesired = complexArray([1.637 + 3.59j, 1.737 + 3.99j, 2.03 + 4.60j])
        nkCalculated = self.silicon._n[0:3]
        assertAlmostEqual(nkDesired, nkCalculated, errorMessage = "material: nk", absoluteTolerance=1e-5)

        erDesired = complexArray([-10.208331+11.75366j, -12.902931 + 13.86126j, -17.0391 + 18.676j])
        urDesired = complexArray([1, 1, 1])
        erCalculated = self.silicon._er[0:3]
        urCalculated = self.silicon._ur[0:3]
        assertAlmostEqual(erDesired, erCalculated, errorMessage = "material: er", absoluteTolerance=1e-3)
        assertAlmostEqual(urDesired, urCalculated, errorMessage = "material: ur", absoluteTolerance=1e-3)

    def testLoadFromDatabase(self):
        wavelength = 0.1879 # Should load from Johnson nk table
        source = Source(wavelength=wavelength)
        Ag = Material('Ag', source=source)
        n_desired = 1.07+1.212j
        n_observed = Ag.n
        assertAlmostEqual(n_desired, n_observed, absoluteTolerance=1e-3)

        source.wavelength=0.2262
        n_desired = 1.26+1.344j
        n_observed = Ag.n
        assertAlmostEqual(n_desired, n_observed, absoluteTolerance=1e-3)


    def testnk(self):
        # Case 1: wavelength is exactly identical to one we have in database
        nDesired = 1.737 + 3.9932j
        self.silicon.source.wavelength = 0.26
        nCalculated = self.silicon.n
        assertAlmostEqual(nCalculated, nDesired, errorMessage="material: testnk: n1", absoluteTolerance=1e-3)


        # Case 3: wavelength is larger than one we have in database - extrapolate linearly
        self.silicon.source.wavelength = 1.47
        nDesired = 3.485 + 1.09e-13j
        nCalculated = self.silicon.n
        assertAlmostEqual(nCalculated, nDesired, errorMessage="material: testnk: n3", absoluteTolerance=1e-3)

        # Case 4: wavelength is smaller than one we have in database - extrapolate in the opposite direction
        self.silicon.source.wavelength = 0.23
        nDesired = 1.437 + 2.7803j
        nCalculated = self.silicon.n
        assertAlmostEqual(nCalculated, nDesired, errorMessage="material: testnk: n4", absoluteTolerance=1e-3)

    def testnkInterpolate(self):
        # Case 2: wavelength located between two points in the database - interpolate 
        nDesired = 4.2620000000000005 + 0.0461865j
        self.silicon.source.wavelength = 0.505
        nCalculated = self.silicon.n
        assertAlmostEqual(nCalculated, nDesired, errorMessage="material: testnk: n21", absoluteTolerance=1e-3)

        # Case 2 repeated: wavelength located between two points in the database - interpolate 
        nDesired = 4.2485 + 0.0450088J
        self.silicon.source.wavelength = 0.5075
        nCalculated = self.silicon.n
        assertAlmostEqual(nCalculated, nDesired, errorMessage="material: testnk: n22", absoluteTolerance=1e-3)

        # Case 2 software is having issues with: wavelength of 0.495 to 0.496. First 0.495
        nDesired = 4.319 + 0.051176j
        self.silicon.source.wavelength = 0.495
        nCalculated = self.silicon.n
        assertAlmostEqual(nCalculated, nDesired, errorMessage="material: testnk: n23", absoluteTolerance=1e-3)

        # Now 0.496
        nDesired = 4.313 + 0.0506492j
        self.silicon.source.wavelength = 0.496
        nCalculated = self.silicon.n
        assertAlmostEqual(nCalculated, nDesired, errorMessage="material: testnk: n24", absoluteTolerance=1e-3)

    """ Tests for discontinuities in the Aspnes data, which I haven't found in Schinke for Si"""
    def testAvoidDiscontinuities(self):
        source = Source(wavelength = 0.495)
        silicon = Material(filename='main/Si/Aspnes.yml', source=source)

        n_desired = 4.325778947368422 + 0.07380526315789473j
        n_observed = silicon.n
        assertAlmostEqual(n_observed, n_desired, errorMessage="material: testnk: n25", absoluteTolerance=1e-6)

        source.wavelength = 0.496
        n_desired = 4.319492753623189 + 0.07293719806763285j
        n_observed = silicon.n
        assertAlmostEqual(n_observed, n_desired, errorMessage="material: testnk: n26", absoluteTolerance=1e-6)

    def testEr(self):
        # Case 1: wavelength is exactly identical to one we have in database
        erDesired = sq(1.737 + 3.9932j)
        self.silicon.source.wavelength = 0.26
        erCalculated = self.silicon.er
        assertAlmostEqual(erCalculated, erDesired, errorMessage="material: testnk: er1")

        # Case 2: wavelength interpolated between two points in the database
        self.silicon.source.wavelength = 0.264
        erDesired = -14.557277+15.787005j
        erCalculated = self.silicon.er
        assertAlmostEqual(erCalculated, erDesired, errorMessage="material: testnk: er2", absoluteTolerance=1e-5)

        # Case 3: wavelength is larger than one we have in database - extrapolate linearly
        self.silicon.source.wavelength = 1.47
        erDesired = sq(3.487 + 1.09e-13j) - 0.01395
        erCalculated = self.silicon.er
        assertAlmostEqual(erCalculated, erDesired, absoluteTolerance=1e-5, errorMessage="material: testnk: er3")

        # Case 4: wavelength is smaller than one we have in database
        self.silicon.source.wavelength = 0.23
        erDesired = -4.744348+7.505422j
        erCalculated = self.silicon.er
        assertAlmostEqual(erCalculated, erDesired, absoluteTolerance = 1e-5,errorMessage="material: testnk: er4")

    def testUr(self):
        # Case 1: wavelength is exactly identical to one we have in database
        urDesired = 1
        self.silicon.wavelength = 0.27
        urCalculated = self.silicon.ur
        assertAlmostEqual(urCalculated, urDesired, errorMessage="material: testnk: ur1")

        # Case 2: wavelength is nearly identical to one we have in database
        wavelength = 0.264
        urDesired = 1
        urCalculated = self.silicon.ur
        assertAlmostEqual(urCalculated, urDesired, errorMessage="material: testnk: ur2")

        # Case 3: wavelength is larger than one we have in database - extrapolate linearly
        wavelength = 1.47
        urDesired = 1
        urCalculated = self.silicon.ur
        assertAlmostEqual(urCalculated, urDesired, absoluteTolerance=1e-5, errorMessage="material: testnk: ur3")

        # Case 4: wavelength is smaller than one we have in database
        wavelength = 0.23
        urDesired = 1
        urCalculated = self.silicon.ur
        assertAlmostEqual(urCalculated, urDesired, absoluteTolerance = 1e-5,errorMessage="material: testnk: ur4")

    def test_extract_dispersion_formula_2(self):
        src = Source(wavelength=0.5)
        SiO2 = Material(filename='main/SiO2/Ghosh-e.yml', source=src)
        n_desired = 1.5580
        n_actual = SiO2.n
        assertAlmostEqual(n_actual, n_desired, absoluteTolerance=1e-5)

    def test_extract_dispersion_formula_1(self):
        src = Source(wavelength=0.5)
        SiO2 = Material(filename='main/SiO2/Radhakrishnan-o.yml', source=src)
        n_desired = 1.548755
        n_actual = SiO2.n
        assertAlmostEqual(n_actual, n_desired, absoluteTolerance=1e-5)

    @classmethod
    def setUpClass(cls):
        source = Source(wavelength=1)
        cls.silicon = Material(filename='main/Si/Schinke.yml',source=source)
