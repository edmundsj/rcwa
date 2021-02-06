"""
Run all the core unit tests, not the lengthy and major integration tests
"""
import context
import unittest
import sys

import testMaterial

loader = unittest.TestLoader()
suite = unittest.TestSuite()
suite.addTests(loader.loadTestsFromModule(testMaterial))

runner = unittest.TextTestRunner(verbosity=0)
result = runner.run(suite)
numberFailures = len(result.errors)
numberErrors= len(result.failures)
numberIssues = numberFailures + numberErrors

sys.exit(numberIssues)
