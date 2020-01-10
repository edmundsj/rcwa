# It looks like the unit testing framework unit-test is going to cause me more headaches than
# it's worth, so I'm just going to write my own tests.
import unittest
import sys


sys.path.append('../core')

# I believe this assumes our function will be called from the top directory.
# Apparently we need to prepend all of our test files with the text 'test' and then
# python's unit testing functionality can automatically discover and run all of our tests.

from matrices import *

class TestCoreMatrices(unittest.TestCase):

    def setUp(self):
        self.NUM_PARAMS = 7;
        self.theta_test = 0;
        self.phi_test = 0;
        self.er_test = 2;
        self.ur_test = 2;
        self.n_test = sqrt(self.ur_test * self.er_test);
        self.kx_test = self.n_test * sin(self.theta_test) * cos(self.phi_test);
        self.ky_test = self.n_test * sin(self.theta_test) * sin(self.phi_test);


    def test_PMatrix(self):
        actual_val = 1 / self.er_test * np.array([
            [self.kx_test * self.ky_test, self.ur_test * self.er_test - sq(self.kx_test)],
            [sq(self.ky_test) - self.ur_test * self.er_test, -self.kx_test * self.ky_test]]);

        P_val = Pi_gen(self.kx_test, self.ky_test, self.er_test, self.ur_test);
        np.testing.assert_array_equal(P_val, actual_val);

    def test_Qmatrix(self):
        actual_val = 1 / self.ur_test * np.array([
            [self.kx_test * self.ky_test, self.ur_test * self.er_test - sq(self.kx_test)],
            [sq(self.ky_test) - self.er_test * self.ur_test, - self.kx_test * self.ky_test]]);

        Q_val = Qi_gen(self.kx_test, self.ky_test, self.er_test, self.ur_test);
        np.testing.assert_array_equal(Q_val, actual_val);

if __name__ == '__main__':
    # In the future, we may want to try the same test functions for a bunch of different values, and so
    # we may not want to call the main function, but I'm not sure. 
    unittest.main();
