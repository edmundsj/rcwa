"""
Run all the core unit tests, not the lengthy and major integration tests
"""
import context
import unittest
import sys

import testCrystal
import testHarmonics
import testLayer
import testMaterial
import testMatrices1x1Harmonics
import testMatrices3x3Harmonics
import testNetlistParser
import testSolver1x1Harmonics
import testSolver3x3Harmonics
import testSource

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
numberFailures = len(result.errors)
numberErrors= len(result.failures)
numberIssues = numberFailures + numberErrors

sys.exit(numberIssues)
