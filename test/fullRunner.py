"""
Run all the core unit tests, not the lengthy and major integration tests
"""
import unittest
import sys
sys.path.append('test')

import testCrystal as testCrystal
import testHarmonics as testHarmonics
import testLayer as testLayer
import testMaterial as testMaterial
import testMatrices1x1Harmonics as testMatrices1x1Harmonics
import testMatrices3x3Harmonics as testMatrices3x3Harmonics
import testNetlistParser as testNetlistParser
import testSolver1x1Harmonics as testSolver1x1Harmonics
import testSolver3x3Harmonics as testSolver3x3Harmonics
import testSource as testSource

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
