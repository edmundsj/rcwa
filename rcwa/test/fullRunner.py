"""
Run all the core unit tests, not the lengthy and major integration tests
"""
import sys
import unittest
import rcwa

import rcwa.test.test_crystal as testCrystal
import rcwa.test.test_harmonics as testHarmonics
import rcwa.test.test_layer as testLayer
import rcwa.test.test_material as testMaterial
import rcwa.test.test_matrices_1x1_harmonics as testMatrices1x1Harmonics
import rcwa.test.test_matrices_3x3_harmonics as testMatrices3x3Harmonics
import rcwa.test.test_netlist_parser as testNetlistParser
import rcwa.test.test_solver_1x1_harmonics as testSolver1x1Harmonics
import rcwa.test.test_solver_3x3_harmonics as testSolver3x3Harmonics
import rcwa.test.test_source as testSource

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
