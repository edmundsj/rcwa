"""
Run all the core unit tests, not the lengthy and major integration tests
"""
import unittest
import sys
sys.path.append('source')

import test.testCrystal as testCrystal
import test.testHarmonics as testHarmonics
import test.testLayer as testLayer
import test.testMaterial as testMaterial
import test.testMatrices1x1Harmonics as testMatrices1x1Harmonics
import test.testMatrices3x3Harmonics as testMatrices3x3Harmonics
import test.testNetlistParser as testNetlistParser
import test.testSolver1x1Harmonics as testSolver1x1Harmonics
import test.testSource as testSource

loader = unittest.TestLoader()
suite = unittest.TestSuite()
suite.addTests(loader.loadTestsFromModule(testCrystal))
suite.addTests(loader.loadTestsFromModule(testHarmonics))
suite.addTests(loader.loadTestsFromModule(testLayer))
suite.addTests(loader.loadTestsFromModule(testMaterial))
suite.addTests(loader.loadTestsFromModule(testSource))
suite.addTests(loader.loadTestsFromModule(testMatrices1x1Harmonics))
suite.addTests(loader.loadTestsFromModule(testMatrices3x3Harmonics))
suite.addTests(loader.loadTestsFromModule(testNetlistParser))
suite.addTests(loader.loadTestsFromModule(testSolver1x1Harmonics))

runner = unittest.TextTestRunner(verbosity=3)
result = runner.run(suite)
