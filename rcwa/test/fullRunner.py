"""
Run all the core unit tests, not the lengthy and major integration tests
"""
import sys
import unittest
import rcwa

import rcwa.test.testCrystal as testCrystal
import rcwa.test.testHarmonics as testHarmonics
import rcwa.test.testLayer as testLayer
import rcwa.test.testMaterial as testMaterial
import rcwa.test.testMatrices1x1Harmonics as testMatrices1x1Harmonics
import rcwa.test.testMatrices3x3Harmonics as testMatrices3x3Harmonics
import rcwa.test.testNetlistParser as testNetlistParser
import rcwa.test.testSolver1x1Harmonics as testSolver1x1Harmonics
import rcwa.test.testSolver3x3Harmonics as testSolver3x3Harmonics
import rcwa.test.testSource as testSource

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
suite.addTests(loader.loadTestsFromModule(testSolver3x3Harmonics))

runner = unittest.TextTestRunner(verbosity=3)
result = runner.run(suite)
numberFailures = len(result.errors)
numberErrors= len(result.failures)
numberIssues = numberFailures + numberErrors

sys.exit(numberIssues)
