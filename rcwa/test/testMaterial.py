import context

import unittest
from material import Material
from source import Source
from shorthand import *
from shorthandTest import *

class TestMaterial(unittest.TestCase):

	def testParseCSV(self):
		wavelengthsToTest = np.array([0.25, 0.26, 0.27])
		nkDesired = complexArray([1.637 + 3.59j, 1.737 + 3.99j, 2.03 + 4.60j])
		nkCalculated = self.silicon._n[0:3]
		assertEqual(nkDesired, nkCalculated, errorMessage = "material: nk")

		erDesired = complexArray([-10.208331+11.75366j, -12.902931 + 13.86126j, -17.0391 + 18.676j])
		urDesired = complexArray([1, 1, 1])
		erCalculated = self.silicon._er[0:3]
		urCalculated = self.silicon._ur[0:3]
		assertAlmostEqual(erDesired, erCalculated, errorMessage = "material: er", absoluteTolerance=1e-3)
		assertAlmostEqual(urDesired, urCalculated, errorMessage = "material: ur", absoluteTolerance=1e-3)

	def testnk(self):
		# Case 1: wavelength is exactly identical to one we have in database
		nDesired = 1.737 + 3.99j
		self.silicon.source.wavelength = 0.26
		nCalculated = self.silicon.n
		assertAlmostEqual(nCalculated, nDesired, errorMessage="material: testnk: n1", absoluteTolerance=1e-3)

		# Case 2: wavelength interpolated between two points in the database
		nDesired = 1.8542+4.234j
		self.silicon.source.wavelength = 0.264
		nCalculated = self.silicon.n
		assertAlmostEqual(nCalculated, nDesired, errorMessage="material: testnk: n2", absoluteTolerance=1e-3)

		# Case 3: wavelength is larger than one we have in database - extrapolate linearly
		self.silicon.source.wavelength = 1.47
		nDesired = 3.485 + 1.09e-13j
		nCalculated = self.silicon.n
		assertAlmostEqual(nCalculated, nDesired, errorMessage="material: testnk: n3", absoluteTolerance=1e-3)

		# Case 4: wavelength is smaller than one we have in database
		self.silicon.source.wavelength = 0.23
		nDesired = 1.437 + 2.79j
		nCalculated = self.silicon.n
		assertEqual(nCalculated, nDesired, errorMessage="material: testnk: n4", absoluteTolerance=1e-3)

	def testEr(self):
		# Case 1: wavelength is exactly identical to one we have in database
		erDesired = sq(1.737 + 3.99j)
		self.silicon.source.wavelength = 0.26
		erCalculated = self.silicon.er
		assertAlmostEqual(erCalculated, erDesired, errorMessage="material: testnk: er1")

		# Case 2: wavelength interpolated between two points in the database
		self.silicon.source.wavelength = 0.264
		erDesired = -14.557399+15.787156j
		erCalculated = self.silicon.er
		assertAlmostEqual(erCalculated, erDesired, errorMessage="material: testnk: er2", absoluteTolerance=1e-5,relativeTolerance=1e-4)

		# Case 3: wavelength is larger than one we have in database - extrapolate linearly
		self.silicon.source.wavelength = 1.47
		erDesired = sq(3.487 + 1.09e-13j) - 0.01395
		erCalculated = self.silicon.er
		assertAlmostEqual(erCalculated, erDesired, absoluteTolerance=1e-5, errorMessage="material: testnk: er3", absoluteTolerance=1e-5, relativeTolerance=1e-4)

		# Case 4: wavelength is smaller than one we have in database
		self.silicon.source.wavelength = 0.23
		erDesired = sq(1.637 + 3.59j) + 5.3892 - 4.2152j
		erCalculated = self.silicon.er
		assertAlmostEqual(erCalculated, erDesired, absoluteTolerance = 1e-5,errorMessage="material: testnk: er4", absoluteTolerance=1e-5, relativeTolerance=1e-4)

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

	@classmethod
	def setUpClass(cls):
		source = Source(wavelength=1)
		cls.silicon = Material(material_name='Si',source=source)
