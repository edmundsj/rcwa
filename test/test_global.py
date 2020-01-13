# I decided to write my own unit tests rather than use python's unit testing framework because
# it was causing me more trouble than it was worth, since it doesn't have a built-in capability
# to use multiple datasets to run the same test. It's extremely annoying.

import sys
sys.path.append('core');

from matrices import *
from fresnel import *

statuses = [];
messages = [];

class TestClass:
    def __init__(self):
        self.messages = []; # Messages sent back from our tests (strings)
        self.statuses = []; # statuses sent back from our tests (boolean)
        self.unit_tests_enabled = True;
        self.integration_tests_enabled = False;

    def printResults(self):
        for s, i in zip(self.statuses, range(len(self.statuses))):
            if(s == False):
                print(self.messages[i]);
        print(f"{self.statuses.count(True)} PASSED, {self.statuses.count(False)} FAILED");

    def testCaller(self, testFunction, *args):
        """
        Handles the actual calling of test functions, manages them with try/catch blocks. Maybe
        not the most elegant way to do things, but the best idea I currently have without wasting
        an inordinate amount of time.
        """
        test_status = False; # By default assume we failed the test.
        test_message = f"{testFunction.__name__}({args}): ";

        try:
            print(f"Calling function {testFunction.__name__} ... ", end=" ");
            testFunction(*args);
            print("OK");
            test_status = True;
            self.statuses.append(test_status);
            self.messages.append(test_message);
        except AssertionError as ae:
            print("FAIL");
            test_message += "FAILED";
            test_message += str(ae);
            self.statuses.append(test_status);
            self.messages.append(test_message);

    def runUnitTests(self):
        print("--------- RUNNING UNIT TESTS... ----------");
        self.testCaller(self.testKz);
        self.testCaller(self.testQMatrix);
        self.testCaller(self.testOMatrix);
        self.testCaller(self.testVMatrix);
        self.testCaller(self.testXMatrix);
        self.testCaller(self.testAMatrix);
        #self.testCaller(testAMatrix);
        #self.testCaller(testBMatrix);
        #self.testCaller(testDMatrix);

    def runIntegrationTests(self):
        """
        Runs integration tests (to test s parameters and reflectance and their composition and stuff)
        """
        theta = 0.3;
        phi = 0.0;
        er = 3;
        ur = 2;

        n1 = sqrt(ur * er);
        kx_n = n1 * sin(theta) * cos(phi);
        ky_n = n1 * sin(theta) * sin(phi);
        kz1_n = sqrt(sq(n1) - sq(kx_n) - sq(ky_n));

        aTE, aTM = aTEM_gen(kx_n, ky_n, kz1_n);
        ATEM = np.transpose(np.array([aTM, aTE]));

        testCaller(testS11Ref, kx_n, ky_n, er, ur, ATEM);
        testCaller(testS22Ref, kx_n, ky_n, er, ur, ATEM);
        testCaller(testS11Trn, kx_n, ky_n, er, ur, ATEM);
        testCaller(testS22Trn, kx_n, ky_n, er, ur, ATEM);
        testCaller(testS21Ref, kx_n, ky_n, er, ur, ATEM);
        testCaller(testS21Trn, kx_n, ky_n, er, ur, ATEM);


    # BEGIN UNIT TESTS SECTION OF THE CLASS
    def testKz(self):
        """
        Tests that we are correctly computing the z-component of our wavevector given some
        material values and a kx and ky.
        """
        abstol = 0.0001;
        reltol = 0.001;
        kx = 1.0006;
        ky = 0.4247;

        # First, we have some data for layer 1
        er = 2.0;
        ur = 1.0;
        kz_n_actual = 0.9046;
        kz_n_calc = kzGen(kx, ky, er, ur);
        np.testing.assert_allclose(kz_n_actual, kz_n_calc, rtol=reltol, atol=abstol);

        # Now, we have some data for layer 2.
        er = 1.0;
        ur = 3.0;

        kz_n_actual = 1.3485;
        kz_n_calc = kzGen(kx, ky, er, ur);
        np.testing.assert_allclose(kz_n_actual, kz_n_calc, rtol=reltol, atol=abstol);

    def testVMatrix(self):
        # The data we have available is only accurate to the 4th decimal place. This should
        # be sufficient. kx and ky are given in the setup, fixed by our angles theta and phi.
        abstol = 0.0001;

    def testQMatrix(self):
        """
        Tests our Q matrix (the matrix we get from the differential equation in terms of Hx and Hy).
        Uses data available on Raymond Rumpf's website on computational electromagnetics.
        """
        # The data we have available is only accurate to the 4th decimal place. This should
        # be sufficient. kx and ky are given in the setup, fixed by our angles theta and phi.
        abstol = 0.0001;
        reltol = 0.001;
        kx = 1.0006;
        ky = 0.4247;

        # Zeroth, we actually have data for our gap layer
        er = 1.0 + sq(kx) + sq(ky);
        ur = 1.0;
        Q_actual = np.array([[0.4250, 1.1804],[-2.0013, -0.4250]], dtype=np.cdouble);
        Q_calc = Qi_gen(kx, ky, er, ur);
        np.testing.assert_allclose(Q_actual, Q_calc, rtol=reltol, atol=abstol);

        # First, we have some data for layer 1
        er = 2.0;
        ur = 1.0;
        Q_actual = np.array([[0.4250, 0.9987],[-1.8196, -0.4250]], dtype=np.cdouble);
        Q_calc = Qi_gen(kx, ky, er, ur);
        np.testing.assert_allclose(Q_actual, Q_calc, rtol=reltol, atol=abstol);

        # Now, we have some data for layer 2.
        er = 1.0;
        ur = 3.0;

        Q_actual = np.array([[0.1417, 0.6662],[-0.9399, -0.1417]], dtype=np.cdouble);
        Q_calc = Qi_gen(kx, ky, er, ur);
        np.testing.assert_allclose(Q_actual, Q_calc, rtol=reltol, atol=abstol);

    def testOMatrix(self):
        """
        Tests the omega matrix (aka the lambda matrix).
        """
        abstol = 0.0001;
        reltol = 0.001;
        kx = 1.0006;
        ky = 0.4247;

        # LAYER 1 DATA
        er = 2.0;
        ur = 1.0;
        kz = 0.9046;

        O_actual = np.array([[0 + 0.9046j, 0+0j],[0+0j,0+0.9046j]], dtype=np.cdouble);
        O_calc = Omega_gen(kz);
        np.testing.assert_allclose(O_actual, O_calc, rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        er = 1.0;
        ur = 3.0;
        kz = 1.3485;

        O_actual = np.array([[0 + 1.3485j, 0+0j],[0+0j,0+1.3485j]], dtype=np.cdouble);
        O_calc = Omega_gen(kz);
        np.testing.assert_allclose(O_actual, O_calc, rtol=reltol, atol=abstol);

    def testVMatrix(self):
        """
        Tests the V matrix (the matrix that relates the magnetic field and the electric field modes)
        """
        abstol = 0.0001;
        reltol = 0.001;
        kx = 1.0006;
        ky = 0.4247;

        # GAP DATA
        er = 1.0 + sq(kx) + sq(ky); # Apparently this is a convenient gap choice.
        ur = 1.0;
        kz = 1.0; # er being the above value makes this true

        (V_calc, W) = VWX_gen(kx, ky, kz, er, ur);
        V_actual = np.array([[0-0.4250j, 0 - 1.1804j], [0 + 2.0013j, 0 + 0.4250j]],dtype=np.cdouble);
        np.testing.assert_allclose(V_actual, V_calc, rtol=reltol, atol=abstol);

        # LAYER 1 DATA
        er = 2.0;
        ur = 1.0;
        kz = 0.9046;

        (V_calc, W) = VWX_gen(kx, ky, kz, er, ur);
        V_actual = np.array([[0-0.4698j,0-1.1040j],[0+2.0114j,0+0.4698j]], dtype=np.cdouble);
        np.testing.assert_allclose(V_actual, V_calc, rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        er = 1.0;
        ur = 3.0;
        kz = 1.3485;

        (V_calc, W) = VWX_gen(kx, ky, kz, er, ur);
        V_actual = np.array([[0-0.1051j,0-0.4941j],[0+0.6970j,0+0.1051j]], dtype=np.cdouble);
        np.testing.assert_allclose(V_actual, V_calc, rtol=reltol, atol=abstol);

    def testXMatrix(self):
        """
        Tests the X matrix (the matrix exponential of the omega matrix)
        """
        abstol = 0.0001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.001; # Relative error tolerance (probably not necessary)
        kx = 1.0006;    # x component of k vector
        ky = 0.4247;    # y component of k vector
        l0 = 2.7;       # Free-space wavelength
        k0 = 2.3271;    # Free-space wavenumber

        # LAYER 1 DATA
        er = 2.0;
        ur = 1.0;
        kz = 0.9046;
        L = 0.25*l0;

        (V, W, X_calc) = VWX_gen(kx, ky, kz, er, ur, k0, L);
        X_actual = np.array([[0.1493+0.9888j, 0+0j],[0+0j,0.1493+0.9888j]],dtype=np.cdouble);
        np.testing.assert_allclose(X_actual, X_calc, rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        er = 1.0;
        ur = 3.0;
        kz = 1.3485;
        L = 0.5*l0;

        (V, W, X_calc) = VWX_gen(kx, ky, kz, er, ur, k0, L);
        X_actual = np.array([[-0.4583 - 0.8888j, 0+0j],[0+0j, -0.4583 - 0.8888j]],dtype=np.cdouble);
        np.testing.assert_allclose(X_actual, X_calc, rtol=reltol, atol=abstol);

    def testAMatrix(self):
        """
        Tests the A matrix (an intermediate matrix in calculating the scattering parameters of a given layer)
        """
        abstol = 0.0001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.001; # Relative error tolerance (probably not necessary)
        kx = 1.0006;    # x component of k vector
        ky = 0.4247;    # y component of k vector
        l0 = 2.7;       # Free-space wavelength
        k0 = 2.3271;    # Free-space wavenumber

        # LAYER 1 DATA
        er = 2.0;
        ur = 1.0;
        kz = 0.9046;
        L = 0.25*l0;
        W1 = np.identity(2, dtype=np.cdouble);
        Wg = np.identity(2, dtype=np.cdouble);
        V1 = np.array([[0 - 0.4698j, 0 - 1.1040j],[0 + 2.0114j, 0 + 0.4698j]], dtype=np.cdouble);
        Vg = np.array([[0 - 0.4250j, 0 - 1.1804j], [0 + 2.0013j, 0 + 0.4250j]], dtype=np.cdouble);

        A_calc = Aij_gen(W1, Wg, V1, Vg);
        A_actual = np.array([[2.0049, -0.0427], [-0.0427, 2.0873]], dtype=np.cdouble);
        np.testing.assert_allclose(A_actual, A_calc, rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        er = 1.0;
        ur = 3.0;
        kz = 1.3485;
        L = 0.5*l0;

        W2 = np.identity(2, dtype=np.cdouble);
        Wg = np.identity(2, dtype=np.cdouble);
        V2 = np.array([[0 - 0.1051j, 0 - 0.4941j],[0 + 0.6970j, 0 + 0.1051j]], dtype=np.cdouble);
        Vg = np.array([[0 - 0.4250j, 0 - 1.1804j],[0 + 2.0013j, 0 + 0.4250j]], dtype=np.cdouble);

        A_calc = Aij_gen(W2, Wg, V2, Vg);
        A_actual = np.array([[3.8324, 0.2579],[0.2579, 3.3342]], dtype=np.cdouble);
        np.testing.assert_allclose(A_actual, A_calc, rtol=reltol, atol=abstol);

    def testAMatrix(self):
        """
        Tests the B matrix (an intermediate matrix in calculating the scattering parameters for
        a given layer).
        """
        abstol = 0.0001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.001; # Relative error tolerance (probably not necessary)
        kx = 1.0006;    # x component of k vector
        ky = 0.4247;    # y component of k vector
        l0 = 2.7;       # Free-space wavelength
        k0 = 2.3271;    # Free-space wavenumber

        # LAYER 1 DATA
        er = 2.0;
        ur = 1.0;
        kz = 0.9046;
        L = 0.25*l0;
        W1 = np.identity(2, dtype=np.cdouble);
        Wg = np.identity(2, dtype=np.cdouble);
        V1 = np.array([[0 - 0.4698j, 0 - 1.1040j],[0 + 2.0114j, 0 + 0.4698j]], dtype=np.cdouble);
        Vg = np.array([[0 - 0.4250j, 0 - 1.1804j], [0 + 2.0013j, 0 + 0.4250j]], dtype=np.cdouble);

        A_calc = Aij_gen(W1, Wg, V1, Vg);
        A_actual = np.array([[2.0049, -0.0427], [-0.0427, 2.0873]], dtype=np.cdouble);
        np.testing.assert_allclose(A_actual, A_calc, rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        er = 1.0;
        ur = 3.0;
        kz = 1.3485;
        L = 0.5*l0;

        W2 = np.identity(2, dtype=np.cdouble);
        Wg = np.identity(2, dtype=np.cdouble);
        V2 = np.array([[0 - 0.1051j, 0 - 0.4941j],[0 + 0.6970j, 0 + 0.1051j]], dtype=np.cdouble);
        Vg = np.array([[0 - 0.4250j, 0 - 1.1804j],[0 + 2.0013j, 0 + 0.4250j]], dtype=np.cdouble);

        A_calc = Aij_gen(W2, Wg, V2, Vg);
        A_actual = np.array([[3.8324, 0.2579],[0.2579, 3.3342]], dtype=np.cdouble);
        np.testing.assert_allclose(A_actual, A_calc, rtol=reltol, atol=abstol);

def main():
    test_class = TestClass(); # Create a new test class
    if(test_class.unit_tests_enabled == True):
        test_class.runUnitTests();
    if(test_class.integration_tests_enabled == True):
        test_class.runIntegrationTests();
    test_class.printResults();

def testaTEM():
    # First, we want to test the case where theta = 0, phi = 0;
    kx_n = 0;
    ky_n = 0;
    kz_n = 4;

    aTE_calc, aTM_calc = aTEM_gen(kx_n, ky_n, kz_n);
    aTE_actual = np.array([0, 1]);
    aTM_actual = np.array([1, 0]);

    np.testing.assert_array_equal(aTE_actual, aTE_calc);
    np.testing.assert_array_equal(aTM_actual, aTM_calc);

# This appears to be working
def testS11Ref(kx_n, ky_n, er1, ur1, ATEM):
    """
    Tests only the reference matrix. This is easy to do because we know it should just
    be the same thing as the composite interface matrix but with the second media set to that of
    free space.
    kx_n: Normalized x component of k vector
    ky_n: Normalized y component of k vector
    ATEM: Matrix that converts TE/TM polarized light to x/y polarized light. Columns are the aTE aTM unit vectors
    """
    er2 = 1;
    ur2 = 1;
    S11_actual, S12_actual, S21_actual, S22_actual = \
            fresnelSMatrixInterface(kx_n, ky_n, er1, er2, ur1, ur2)

    reltol = 1e-10;
    abstol = 1e-10;
    erg = 1;

def main():
    test_class = TestClass(); # Create a new test class
    if(test_class.unit_tests_enabled == True):
        test_class.runUnitTests();
    if(test_class.integration_tests_enabled == True):
        test_class.runIntegrationTests();
    test_class.printResults();

def testaTEM():
    # First, we want to test the case where theta = 0, phi = 0;
    kx_n = 0;
    ky_n = 0;
    kz_n = 4;

    aTE_calc, aTM_calc = aTEM_gen(kx_n, ky_n, kz_n);
    aTE_actual = np.array([0, 1]);
    aTM_actual = np.array([1, 0]);

    np.testing.assert_array_equal(aTE_actual, aTE_calc);
    np.testing.assert_array_equal(aTM_actual, aTM_calc);

# This appears to be working
def testS11Ref(kx_n, ky_n, er1, ur1, ATEM):
    """
    Tests only the reference matrix. This is easy to do because we know it should just
    be the same thing as the composite interface matrix but with the second media set to that of
    free space.
    kx_n: Normalized x component of k vector
    ky_n: Normalized y component of k vector
    ATEM: Matrix that converts TE/TM polarized light to x/y polarized light. Columns are the aTE aTM unit vectors
    """
    er2 = 1;
    ur2 = 1;
    S11_actual, S12_actual, S21_actual, S22_actual = \
            fresnelSMatrixInterface(kx_n, ky_n, er1, er2, ur1, ur2)

    reltol = 1e-10;
    abstol = 1e-10;
    erg = 1;
    urg = 1;
    Vg, Wg = VWX_gen(kx_n, ky_n, erg, urg);

    #print(f"\nATEM: {ATEM}");
    (Sref, Wref) = SWref_gen(kx_n, ky_n, er1, ur1, Wg, Vg);

    S11_calc = Sref[0,0];
    S11_calc = inv(ATEM) @ S11_calc @ ATEM

    np.testing.assert_allclose(S11_actual, S11_calc, rtol=reltol, atol=abstol);

def testS22Ref(kx_n, ky_n, er1, ur1, ATEM):
    """
    Tests only the reference matrix. This is easy to do because we know it should just
    be the same thing as the composite interface matrix but with the second media set to that of
    free space.
    kx_n: Normalized x component of k vector
    ky_n: Normalized y component of k vector
    ATEM: Matrix that converts TE/TM polarized light to x/y polarized light. Columns are the aTE aTM unit vectors
    """
    er2 = 1;
    ur2 = 1;
    S11_actual, S12_actual, S21_actual, S22_actual = \
            fresnelSMatrixInterface(kx_n, ky_n, er1, er2, ur1, ur2)

    reltol = 1e-10;
    abstol = 1e-10;
    erg = 1;
    urg = 1;
    Vg, Wg = VWX_gen(kx_n, ky_n, erg, urg);

    #print(f"\nATEM: {ATEM}");
    (Sref, Wref) = SWref_gen(kx_n, ky_n, er1, ur1, Wg, Vg);

    S22_calc = Sref[1,1];
    S22_calc = inv(ATEM) @ S22_calc @ ATEM

    np.testing.assert_allclose(S22_actual, S22_calc, rtol=reltol, atol=abstol);

def testS21Ref(kx_n, ky_n, er1, ur1, ATEM):
    """
    Tests only the reference matrix. This is easy to do because we know it should just
    be the same thing as the composite interface matrix but with the second media set to that of
    free space.
    kx_n: Normalized x component of k vector
    ky_n: Normalized y component of k vector
    ATEM: Matrix that converts TE/TM polarized light to x/y polarized light. Columns are the aTE aTM unit vectors
    """
    er2 = 1;
    ur2 = 1;
    S11_actual, S12_actual, S21_actual, S22_actual = \
            fresnelSMatrixInterface(kx_n, ky_n, er1, er2, ur1, ur2)

    reltol = 1e-10;
    abstol = 1e-10;
    erg = 1;
    urg = 1;
    Vg, Wg = VWX_gen(kx_n, ky_n, erg, urg);

    #print(f"\nATEM: {ATEM}");
    (Sref, Wref) = SWref_gen(kx_n, ky_n, er1, ur1, Wg, Vg);

    S21_calc = Sref[1,0];
    S21_calc = inv(ATEM) @ S21_calc @ ATEM

    np.testing.assert_allclose(S21_actual, S21_calc, rtol=reltol, atol=abstol);

def testS11Trn(kx_n, ky_n, er2, ur2, ATEM):
    """
    Tests only the reference matrix. This is easy to do because we know it should just
    be the same thing as the composite interface matrix but with the second media set to that of
    free space.
    kx_n: Normalized x component of k vector
    ky_n: Normalized y component of k vector
    ATEM: Matrix that converts TE/TM polarized light to x/y polarized light. Columns are the aTE aTM unit vectors
    """
    er1 = 1;
    ur1 = 1;
    S11_actual, S12_actual, S21_actual, S22_actual = \
            fresnelSMatrixInterface(kx_n, ky_n, er1, er2, ur1, ur2)

    reltol = 1e-10;
    abstol = 1e-10;
    erg = 1;
    urg = 1;
    Vg, Wg = VWX_gen(kx_n, ky_n, erg, urg);

    #print(f"\nATEM: {ATEM}");
    (Strn, Wtrn) = SWtrn_gen(kx_n, ky_n, er2, ur2, Wg, Vg);

    S11_calc = Strn[0,0];
    S11_calc = inv(ATEM) @ S11_calc @ ATEM

    np.testing.assert_allclose(S11_actual, S11_calc, rtol=reltol, atol=abstol);

def testS22Trn(kx_n, ky_n, er2, ur2, ATEM):
    """
    Tests only the reference matrix. This is easy to do because we know it should just
    be the same thing as the composite interface matrix but with the second media set to that of
    free space.
    kx_n: Normalized x component of k vector
    ky_n: Normalized y component of k vector
    ATEM: Matrix that converts TE/TM polarized light to x/y polarized light. Columns are the aTE aTM unit vectors
    """
    er1 = 1;
    ur1 = 1;
    S11_actual, S12_actual, S21_actual, S22_actual = \
            fresnelSMatrixInterface(kx_n, ky_n, er1, er2, ur1, ur2)

    reltol = 1e-10;
    abstol = 1e-10;
    erg = 1;
    urg = 1;
    Vg, Wg = VWX_gen(kx_n, ky_n, erg, urg);

    #print(f"\nATEM: {ATEM}");
    (Strn, Wtrn) = SWtrn_gen(kx_n, ky_n, er2, ur2, Wg, Vg);

    S22_calc = Strn[1,1];
    S22_calc = inv(ATEM) @ S22_calc @ ATEM

    np.testing.assert_allclose(S22_actual, S22_calc, rtol=reltol, atol=abstol);

def testS21Trn(kx_n, ky_n, er2, ur2, ATEM):
    """
    Tests only the reference matrix. This is easy to do because we know it should just
    be the same thing as the composite interface matrix but with the second media set to that of
    free space.
    kx_n: Normalized x component of k vector
    ky_n: Normalized y component of k vector
    ATEM: Matrix that converts TE/TM polarized light to x/y polarized light. Columns are the aTE aTM unit vectors
    """
    er1 = 1;
    ur1 = 1;
    S11_actual, S12_actual, S21_actual, S22_actual = \
            fresnelSMatrixInterface(kx_n, ky_n, er1, er2, ur1, ur2)

    reltol = 1e-10;
    abstol = 1e-10;
    erg = 1;
    urg = 1;
    Vg, Wg = VWX_gen(kx_n, ky_n, erg, urg);

    #print(f"\nATEM: {ATEM}");
    (Strn, Wref) = SWtrn_gen(kx_n, ky_n, er1, ur1, Wg, Vg);

    S21_calc = Strn[1,0];
    S21_calc = inv(ATEM) @ S21_calc @ ATEM

    np.testing.assert_allclose(S21_actual, S21_calc, rtol=reltol, atol=abstol);

main();
