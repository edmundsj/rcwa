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
        self.testCaller(self.testBMatrix);
        self.testCaller(self.testDiMatrix);
        self.testCaller(self.testS11i);
        self.testCaller(self.testS12i);
        self.testCaller(self.testS21i);
        self.testCaller(self.testS22i);
        self.testCaller(self.testDRed);
        self.testCaller(self.testFRed);
        self.testCaller(self.testRedhefferProductS11);
        self.testCaller(self.testRedhefferProductS12);
        self.testCaller(self.testRedhefferProductS21);
        self.testCaller(self.testRedhefferProductS22);
        self.testCaller(self.testSrefFull);
        self.testCaller(self.testStrnFull);

    def runIntegrationTests(self):
        """
        Runs integration tests to verify s-parameters for composite code, to verify the output field
        for a given input field, and to verify the reflectance/transmittance and enforce power 
        conservation.
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

        # First, test the s-parameters

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
        V_actual = np.array([[0 - 0.4250j, 0 - 1.1804j], [0 + 2.0013j, 0 + 0.4250j]],dtype=np.cdouble);
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

        # REFERENCE REGION DATA
        er = 1.4;
        ur = 1.2;
        kz = 0.705995; # Calculated manually using er and ur above.
        (V_calc, W_ref) = VWX_gen(kx, ky, kz, er, ur);
        V_actual = np.array([
            [0 - 0.5017j, 0 - 0.8012j],
            [0 + 1.7702j, 0 + 0.5017j]], dtype=np.cdouble);
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

    def testBMatrix(self):
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

        B_calc = Bij_gen(W1, Wg, V1, Vg);
        B_actual = np.array([[-0.0049, 0.0427],[0.0427, -0.0873]], dtype=np.cdouble);
        np.testing.assert_allclose(B_actual, B_calc, rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        er = 1.0;
        ur = 3.0;
        kz = 1.3485;
        L = 0.5*l0;

        W2 = np.identity(2, dtype=np.cdouble);
        Wg = np.identity(2, dtype=np.cdouble);
        V2 = np.array([[0 - 0.1051j, 0 - 0.4941j],[0 + 0.6970j, 0 + 0.1051j]], dtype=np.cdouble);
        Vg = np.array([[0 - 0.4250j, 0 - 1.1804j],[0 + 2.0013j, 0 + 0.4250j]], dtype=np.cdouble);

        B_calc = Bij_gen(W2, Wg, V2, Vg);
        B_actual = np.array([[-1.8324, -0.2579],[-0.2579, -1.3342]], dtype=np.cdouble);
        np.testing.assert_allclose(B_actual, B_calc, rtol=reltol, atol=abstol);

    def testDiMatrix(self):
        """
        Tests the composite D matrix (one of the matrices we use directly in the calculation
        of scattering matrices. At this point, we have to make a decision. Unfortunately since we
        only have the intermediate matrices to 4 decimal places, and we only have this matrix to
        4 decimal places (and it contains some nearly-zero terms), we are going to incur appreciable
        error. For now, I will tolerate that error, because we have one test case that we can test
        to 4 decimal places.
        """
        abstol = 0.003;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.1; # Relative error tolerance (probably not necessary)
        kx = 1.0006;    # x component of k vector
        ky = 0.4247;    # y component of k vector
        l0 = 2.7;       # Free-space wavelength
        k0 = 2.3271;    # Free-space wavenumber

        # LAYER 1 DATA
        er = 2.0;
        ur = 1.0;
        kz = 0.9046;
        A = np.array([[2.0049, -0.0427], [-0.0427, 2.0873]], dtype=np.cdouble);
        B = np.array([[-0.0049, 0.0427], [0.0427, -0.0873]], dtype=np.cdouble);
        X = np.array([[0.1493 + 0.9888j, 0+0j],[0+0j, 0.4193 + 0.9888j]], dtype=np.cdouble);

        D_calc = DiGen(A, B, X);
        D_actual = np.array([[2.0057 - 0.0003j, -0.0445 + 0.0006j],[-0.0445 + 0.0006j, 2.0916 - 0.0013j]],
                dtype=np.cdouble);
        np.testing.assert_allclose(D_actual, D_calc, rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        # Since now we have the d-matrix to higher precision we can test it more strongly.
        abstol = 0.0001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.001; # Relative error tolerance (probably not necessary)
        er = 1.0;
        ur = 3.0;
        kz = 1.3485;
        L = 0.5*l0;

        A = np.array([[3.8324, 0.2579],[0.2579, 3.3342]], dtype=np.cdouble);
        B = np.array([[-1.8324, -0.2579], [-0.2579, -1.3342]], dtype=np.cdouble);
        X = np.array([[-0.4583 - 0.8888j, 0+0j],[0+0j, -0.4583 - 0.8888j]], dtype=np.cdouble);

        D_calc = DiGen(A, B, X);
        D_actual = np.array([[4.3436 - 0.7182j, 0.3604 - 0.1440j], [0.3604 - 0.1440j, 3.6475 - 0.4401j]], dtype=np.cdouble);
        np.testing.assert_allclose(D_actual, D_calc, rtol=reltol, atol=abstol);

    def testS11i(self):
        """
        Tests the S11 element of an inner layer (the ith layer)
        """
        abstol = 0.03;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 1; # Relative error tolerance (probably not necessary)

        # LAYER 1 DATA
        A = np.array([[2.0049, -0.0427], [-0.0427, 2.0873]], dtype=np.cdouble);
        B = np.array([[-0.0049, 0.0427], [0.0427, -0.0873]], dtype=np.cdouble);
        X = np.array([[0.1493 + 0.9888j, 0+0j],[0+0j, 0.4193 + 0.9888j]], dtype=np.cdouble);
        D = np.array([[2.0057 - 0.0003j, -0.0445 + 0.0006j],[-0.0445 + 0.0006j, 2.0916 - 0.0013j]],
                dtype=np.cdouble);

        S_calc = Si_gen(A, B, X, D)
        S11_calc = S_calc[0,0];
        S11_actual = np.array([[0.0039 - 0.0006j, -0.0398 + 0.0060j],[-0.0398 + 0.0060j, 0.0808 - 0.0121j]],
                dtype=np.cdouble);
        np.testing.assert_allclose(S11_actual, S11_calc, rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        # Since now we have the S-matrix to higher precision we can test it more strongly.
        abstol = 0.0001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.001; # Relative error tolerance (probably not necessary)

        A = np.array([[3.8324, 0.2579],[0.2579, 3.3342]], dtype=np.cdouble);
        B = np.array([[-1.8324, -0.2579], [-0.2579, -1.3342]], dtype=np.cdouble);
        X = np.array([[-0.4583 - 0.8888j, 0+0j],[0+0j, -0.4583 - 0.8888j]], dtype=np.cdouble);
        D = np.array([[4.3436 - 0.7182j, 0.3604 - 0.1440j],[0.3604 - 0.1440j, 3.6475 - 0.4401j]],
            dtype=np.cdouble);

        S_calc = Si_gen(A, B, X, D);
        S11_calc = S_calc[0,0];
        S11_actual = np.array([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]], dtype = np.cdouble);
        np.testing.assert_allclose(S11_actual, S11_calc, rtol=reltol, atol=abstol);

    def testS12i(self):
        """
        Tests the S11 element of an inner layer (the ith layer)
        """
        abstol = 0.001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 1; # Relative error tolerance (probably not necessary)

        # LAYER 1 DATA
        A = np.array([[2.0049, -0.0427], [-0.0427, 2.0873]], dtype=np.cdouble);
        B = np.array([[-0.0049, 0.0427], [0.0427, -0.0873]], dtype=np.cdouble);
        X = np.array([[0.1493 + 0.9888j, 0+0j],[0+0j, 0.4193 + 0.9888j]], dtype=np.cdouble);
        D = np.array([[2.0057 - 0.0003j, -0.0445 + 0.0006j],[-0.0445 + 0.0006j, 2.0916 - 0.0013j]],
                dtype=np.cdouble);

        S_calc = Si_gen(A, B, X, D)
        S12_calc = S_calc[0,1];
        S12_actual = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);
        np.testing.assert_allclose(S12_actual, S12_calc, rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        # Since now we have the S-matrix to higher precision we can test it more strongly.
        abstol = 0.0001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.001; # Relative error tolerance (probably not necessary)

        A = np.array([[3.8324, 0.2579],[0.2579, 3.3342]], dtype=np.cdouble);
        B = np.array([[-1.8324, -0.2579], [-0.2579, -1.3342]], dtype=np.cdouble);
        X = np.array([[-0.4583 - 0.8888j, 0+0j],[0+0j, -0.4583 - 0.8888j]], dtype=np.cdouble);
        D = np.array([[4.3436 - 0.7182j, 0.3604 - 0.1440j],[0.3604 - 0.1440j, 3.6475 - 0.4401j]],
            dtype=np.cdouble);

        S_calc = Si_gen(A, B, X, D);
        S12_calc = S_calc[0,1];
        S12_actual = np.array([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]], dtype=np.cdouble);
        np.testing.assert_allclose(S12_actual, S12_calc, rtol=reltol, atol=abstol);

    def testS21i(self):
        """
        Tests the S11 element of an inner layer (the ith layer)
        """
        abstol = 0.001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 1; # Relative error tolerance (probably not necessary)

        # LAYER 1 DATA
        A = np.array([[2.0049, -0.0427], [-0.0427, 2.0873]], dtype=np.cdouble);
        B = np.array([[-0.0049, 0.0427], [0.0427, -0.0873]], dtype=np.cdouble);
        X = np.array([[0.1493 + 0.9888j, 0+0j],[0+0j, 0.4193 + 0.9888j]], dtype=np.cdouble);
        D = np.array([[2.0057 - 0.0003j, -0.0445 + 0.0006j],[-0.0445 + 0.0006j, 2.0916 - 0.0013j]],
                dtype=np.cdouble);

        S_calc = Si_gen(A, B, X, D)
        S21_calc = S_calc[1,0];
        S21_actual = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);
        np.testing.assert_allclose(S21_actual, S21_calc, rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        # Since now we have the S-matrix to higher precision we can test it more strongly.
        abstol = 0.0001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.001; # Relative error tolerance (probably not necessary)

        A = np.array([[3.8324, 0.2579],[0.2579, 3.3342]], dtype=np.cdouble);
        B = np.array([[-1.8324, -0.2579], [-0.2579, -1.3342]], dtype=np.cdouble);
        X = np.array([[-0.4583 - 0.8888j, 0+0j],[0+0j, -0.4583 - 0.8888j]], dtype=np.cdouble);
        D = np.array([[4.3436 - 0.7182j, 0.3604 - 0.1440j],[0.3604 - 0.1440j, 3.6475 - 0.4401j]],
            dtype=np.cdouble);

        S_calc = Si_gen(A, B, X, D);
        S21_calc = S_calc[1,0];
        S21_actual = np.array([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]], dtype=np.cdouble);
        np.testing.assert_allclose(S21_actual, S21_calc, rtol=reltol, atol=abstol);

    def testS22i(self):
        """
        Tests the S11 element of an inner layer (the ith layer)
        """
        abstol = 0.03;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 1; # Relative error tolerance (probably not necessary)

        # LAYER 1 DATA
        A = np.array([[2.0049, -0.0427], [-0.0427, 2.0873]], dtype=np.cdouble);
        B = np.array([[-0.0049, 0.0427], [0.0427, -0.0873]], dtype=np.cdouble);
        X = np.array([[0.1493 + 0.9888j, 0+0j],[0+0j, 0.4193 + 0.9888j]], dtype=np.cdouble);
        D = np.array([[2.0057 - 0.0003j, -0.0445 + 0.0006j],[-0.0445 + 0.0006j, 2.0916 - 0.0013j]],
                dtype=np.cdouble);

        S_calc = Si_gen(A, B, X, D)
        S22_calc = S_calc[1,1];
        S22_actual = np.array([[0.0039 - 0.0006j, -0.0398 + 0.0060j],[-0.0398 + 0.0060j, 0.0808 - 0.0121j]],
                dtype=np.cdouble);
        np.testing.assert_allclose(S22_actual, S22_calc, rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        # Since now we have the S-matrix to higher precision we can test it more strongly.
        abstol = 0.0001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.001; # Relative error tolerance (probably not necessary)

        A = np.array([[3.8324, 0.2579],[0.2579, 3.3342]], dtype=np.cdouble);
        B = np.array([[-1.8324, -0.2579], [-0.2579, -1.3342]], dtype=np.cdouble);
        X = np.array([[-0.4583 - 0.8888j, 0+0j],[0+0j, -0.4583 - 0.8888j]], dtype=np.cdouble);
        D = np.array([[4.3436 - 0.7182j, 0.3604 - 0.1440j],[0.3604 - 0.1440j, 3.6475 - 0.4401j]],
            dtype=np.cdouble);

        S_calc = Si_gen(A, B, X, D);
        S22_calc = S_calc[1,1];
        S22_actual = np.array([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]], dtype = np.cdouble);
        np.testing.assert_allclose(S22_actual, S22_calc, rtol=reltol, atol=abstol);

    def testDRed(self):
        abstol = 0.0001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.001; # Relative error tolerance (probably not necessary)

        # LAYER 1 DATA - This is the data 
        # Current global data. Before applying the Redheffer star product to
        # any values, S11/S22 should be zero and S12/S21 should be the identity.
        S11A = np.zeros(2);
        S22A = np.zeros(2);
        S12A = np.identity(2);
        S21A = np.identity(2);

        S11B = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.0060j, 0.0808 - 0.0121j]], dtype=np.cdouble);
        S12B = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);
        S21B = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);
        S22B = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.0060j, 0.0808 - 0.0121j]], dtype=np.cdouble);

        Dred_calc = DredGen(S12A, S22A, S11B)
        Dred_actual = np.array([[1,0],[0,1]], dtype=np.cdouble);
        np.testing.assert_allclose(Dred_actual, Dred_calc, rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        # Current global data. After applying the Redheffer star product once, we will end up
        # with a new set of values for the scattering matrix, which are given in the text.
        # They are the SG11/SG12/SG21/SG22 from the first layer (after the update)
        S11A = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)
        S12A = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        S21A = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        S22A = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)

        # The S11/etc. for the second layer.
        S11B = np.array([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]], dtype = np.cdouble);
        S12B = np.array([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]], dtype=np.cdouble);
        S21B = np.array([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]], dtype=np.cdouble);
        S22B = np.array([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]], dtype = np.cdouble);

        Dred_calc = DredGen(S12A, S22A, S11B)
        Dred_actual = np.array([
            [0.1506 + 0.9886j, -0.0163 - 0.0190j],
            [-0.0163 - 0.0190j, 0.1822 + 1.0253j]], dtype=np.cdouble);
        np.testing.assert_allclose(Dred_actual, Dred_calc, rtol=reltol, atol=abstol);

    def testFRed(self):
        abstol = 0.001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.01; # Relative error tolerance (probably not necessary)

        # LAYER 1 DATA - This is the data 
        # Current global data. Before applying the Redheffer star product to
        # any values, S11/S22 should be zero and S12/S21 should be the identity.
        S11A = np.zeros(2);
        S22A = np.zeros(2);
        S12A = np.identity(2);
        S21A = np.identity(2);

        S11B = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.0060j, 0.0808 - 0.0121j]], dtype=np.cdouble);
        S12B = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);
        S21B = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);
        S22B = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.0060j, 0.0808 - 0.0121j]], dtype=np.cdouble);

        Fred_calc = FredGen(S22A, S11B, S21B)
        Fred_actual = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.148 + 0.9848j]], dtype=np.cdouble);
        np.testing.assert_allclose(Fred_actual, Fred_calc, rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        # Current global data. After applying the Redheffer star product once, we will end up
        # with a new set of values for the scattering matrix, which are given in the text.
        # They are the SG11/SG12/SG21/SG22 from the first layer (after the update)
        abstol = 0.0001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.001; # Relative error tolerance (probably not necessary)

        S11A = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)
        S12A = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        S21A = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        S22A = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)

        # The S11/etc. for the second layer.
        S11B = np.array([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]], dtype = np.cdouble);
        S12B = np.array([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]], dtype=np.cdouble);
        S21B = np.array([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]], dtype=np.cdouble);
        S22B = np.array([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]], dtype = np.cdouble);

        Fred_calc = FredGen(S22A, S11B, S21B)
        Fred_actual = np.array([
            [-0.2117 - 0.6413j, 0.0471 + 0.0518j],
            [0.0471 + 0.0518j, -0.3027 - 0.7414j]], dtype=np.cdouble);
        np.testing.assert_allclose(Fred_actual, Fred_calc, rtol=reltol, atol=abstol);

    def testRedhefferProductS11(self):
        abstol = 0.0001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.001; # Relative error tolerance (probably not necessary)

        # LAYER 1 DATA - This is the data 
        # Current global data. Before applying the Redheffer star product to
        # any values, S11/S22 should be zero and S12/S21 should be the identity.
        SA = np.zeros((2,2,2,2));
        SA[0,1] = np.identity(2);
        SA[1,0] = np.identity(2);

        SB = np.zeros((2,2,2,2), dtype=np.cdouble);
        SB[0,0] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.0060j, 0.0808 - 0.0121j]], dtype=np.cdouble);
        SB[0,1] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);
        SB[1,0] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);
        SB[1,1] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.0060j, 0.0808 - 0.0121j]], dtype=np.cdouble);
        D = np.array([
            [1,0],
            [0,1]], dtype=np.cdouble);
        F = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);

        SAB_actual = np.zeros((2,2,2,2), dtype=np.cdouble);
        SAB_actual[0,0] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)
        SAB_actual[0,1] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SAB_actual[1,0] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SAB_actual[1,1] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)
        SAB_calc = redhefferProduct(SA, SB);
        np.testing.assert_allclose(SAB_actual[0,0], SAB_calc[0,0], rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        # Current global data. After applying the Redheffer star product once, we will end up
        # with a new set of values for the scattering matrix, which are given in the text.
        # They are the SG11/SG12/SG21/SG22 from the first layer (after the update)
        SA = np.zeros((2,2,2,2), dtype=np.cdouble);
        SA[0,0] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)
        SA[0,1] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SA[1,0] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SA[1,1] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)

        # The S11/etc. for the second layer.
        SB = np.zeros((2,2,2,2), dtype=np.cdouble);
        SB[0,0] = np.array([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]], dtype = np.cdouble);
        SB[0,1] = np.array([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]], dtype=np.cdouble);
        SB[1,0] = np.array([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]], dtype=np.cdouble);
        SB[1,1] = np.array([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]], dtype = np.cdouble);
        D = np.array([
            [0.1506 + 0.9886j, -0.0163 - 0.0190j],
            [-0.0163 - 0.0190j, 0.1822 + 1.0253j]], dtype=np.cdouble);
        F = np.array([
            [-0.2117 - 0.6413j, 0.0471 + 0.0518j],
            [0.0471 + 0.0518j, -0.3027 - 0.7414j]], dtype=np.cdouble);

        SAB_calc = redhefferProduct(SA, SB);
        # THIS STILL NEEDS TO BE FILLED IN. IT'S THE LAST GLOBAL S-MATRIX.
        SAB_actual = np.zeros((2,2,2,2), dtype=np.cdouble);
        SAB_actual[0,0] = np.array([
            [-0.5961 + 0.4214j, -0.0840 + 0.0085j],
            [-0.0840 + 0.0085j, -0.4339 + 0.4051j]], dtype=np.cdouble);
        SAB_actual[0,1] = np.array([
            [0.6020 - 0.3046j, -0.0431 + 0.0534j],
            [-0.0431 + 0.0534j, 0.6852 - 0.4078j]], dtype=np.cdouble);
        SAB_actual[1,0] = np.array([
            [0.6020 - 0.3046j, -0.0431 + 0.0534j],
            [-0.0431 + 0.0534j, 0.6852 - 0.4078j]], dtype=np.cdouble);
        SAB_actual[1,1] = np.array([
            [0.6971 - 0.2216j, 0.0672 - 0.0211j],
            [0.0672 - 0.0211j, 0.5673 - 0.1808j]], dtype=np.cdouble);

        np.testing.assert_allclose(SAB_actual[0,0], SAB_calc[0,0], rtol=reltol, atol=abstol);

    def testRedhefferProductS12(self):
        abstol = 0.0001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.001; # Relative error tolerance (probably not necessary)

        # LAYER 1 DATA - This is the data 
        # Current global data. Before applying the Redheffer star product to
        # any values, S11/S22 should be zero and S12/S21 should be the identity.
        SA = np.zeros((2,2,2,2));
        SA[0,1] = np.identity(2);
        SA[1,0] = np.identity(2);

        SB = np.zeros((2,2,2,2), dtype=np.cdouble);
        SB[0,0] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.0060j, 0.0808 - 0.0121j]], dtype=np.cdouble);
        SB[0,1] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);
        SB[1,0] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);
        SB[1,1] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.0060j, 0.0808 - 0.0121j]], dtype=np.cdouble);
        D = np.array([
            [1,0],
            [0,1]], dtype=np.cdouble);
        F = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);

        SAB_actual = np.zeros((2,2,2,2), dtype=np.cdouble);
        SAB_actual[0,0] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)
        SAB_actual[0,1] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SAB_actual[1,0] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SAB_actual[1,1] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)
        SAB_calc = redhefferProduct(SA, SB);
        np.testing.assert_allclose(SAB_actual[0,1], SAB_calc[0,1], rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        # Current global data. After applying the Redheffer star product once, we will end up
        # with a new set of values for the scattering matrix, which are given in the text.
        # They are the SG11/SG12/SG21/SG22 from the first layer (after the update)
        SA = np.zeros((2,2,2,2), dtype=np.cdouble);
        SA[0,0] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)
        SA[0,1] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SA[1,0] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SA[1,1] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)

        # The S11/etc. for the second layer.
        SB = np.zeros((2,2,2,2), dtype=np.cdouble);
        SB[0,0] = np.array([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]], dtype = np.cdouble);
        SB[0,1] = np.array([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]], dtype=np.cdouble);
        SB[1,0] = np.array([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]], dtype=np.cdouble);
        SB[1,1] = np.array([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]], dtype = np.cdouble);
        D = np.array([
            [0.1506 + 0.9886j, -0.0163 - 0.0190j],
            [-0.0163 - 0.0190j, 0.1822 + 1.0253j]], dtype=np.cdouble);
        F = np.array([
            [-0.2117 - 0.6413j, 0.0471 + 0.0518j],
            [0.0471 + 0.0518j, -0.3027 - 0.7414j]], dtype=np.cdouble);

        SAB_calc = redhefferProduct(SA, SB);
        # THIS STILL NEEDS TO BE FILLED IN. IT'S THE LAST GLOBAL S-MATRIX.
        SAB_actual = np.zeros((2,2,2,2), dtype=np.cdouble);
        SAB_actual[0,0] = np.array([
            [-0.5961 + 0.4214j, -0.0840 + 0.0085j],
            [-0.0840 + 0.0085j, -0.4339 + 0.4051j]], dtype=np.cdouble);
        SAB_actual[0,1] = np.array([
            [0.6020 - 0.3046j, -0.0431 + 0.0534j],
            [-0.0431 + 0.0534j, 0.6852 - 0.4078j]], dtype=np.cdouble);
        SAB_actual[1,0] = np.array([
            [0.6020 - 0.3046j, -0.0431 + 0.0534j],
            [-0.0431 + 0.0534j, 0.6852 - 0.4078j]], dtype=np.cdouble);
        SAB_actual[1,1] = np.array([
            [0.6971 - 0.2216j, 0.0672 - 0.0211j],
            [0.0672 - 0.0211j, 0.5673 - 0.1808j]], dtype=np.cdouble);

        np.testing.assert_allclose(SAB_actual[0,1], SAB_calc[0,1], rtol=reltol, atol=abstol);

    def testRedhefferProductS21(self):
        abstol = 0.0001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.001; # Relative error tolerance (probably not necessary)

        # LAYER 1 DATA - This is the data 
        # Current global data. Before applying the Redheffer star product to
        # any values, S11/S22 should be zero and S12/S21 should be the identity.
        SA = np.zeros((2,2,2,2));
        SA[0,1] = np.identity(2);
        SA[1,0] = np.identity(2);

        SB = np.zeros((2,2,2,2), dtype=np.cdouble);
        SB[0,0] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.0060j, 0.0808 - 0.0121j]], dtype=np.cdouble);
        SB[0,1] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);
        SB[1,0] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);
        SB[1,1] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.0060j, 0.0808 - 0.0121j]], dtype=np.cdouble);
        D = np.array([
            [1,0],
            [0,1]], dtype=np.cdouble);
        F = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);

        SAB_actual = np.zeros((2,2,2,2), dtype=np.cdouble);
        SAB_actual[0,0] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)
        SAB_actual[0,1] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SAB_actual[1,0] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SAB_actual[1,1] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)
        SAB_calc = redhefferProduct(SA, SB);
        np.testing.assert_allclose(SAB_actual[1,0], SAB_calc[1,0], rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        # Current global data. After applying the Redheffer star product once, we will end up
        # with a new set of values for the scattering matrix, which are given in the text.
        # They are the SG11/SG12/SG21/SG22 from the first layer (after the update)
        SA = np.zeros((2,2,2,2), dtype=np.cdouble);
        SA[0,0] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)
        SA[0,1] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SA[1,0] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SA[1,1] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)

        # The S11/etc. for the second layer.
        SB = np.zeros((2,2,2,2), dtype=np.cdouble);
        SB[0,0] = np.array([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]], dtype = np.cdouble);
        SB[0,1] = np.array([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]], dtype=np.cdouble);
        SB[1,0] = np.array([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]], dtype=np.cdouble);
        SB[1,1] = np.array([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]], dtype = np.cdouble);
        D = np.array([
            [0.1506 + 0.9886j, -0.0163 - 0.0190j],
            [-0.0163 - 0.0190j, 0.1822 + 1.0253j]], dtype=np.cdouble);
        F = np.array([
            [-0.2117 - 0.6413j, 0.0471 + 0.0518j],
            [0.0471 + 0.0518j, -0.3027 - 0.7414j]], dtype=np.cdouble);

        SAB_calc = redhefferProduct(SA, SB);
        # THIS STILL NEEDS TO BE FILLED IN. IT'S THE LAST GLOBAL S-MATRIX.
        SAB_actual = np.zeros((2,2,2,2), dtype=np.cdouble);
        SAB_actual[0,0] = np.array([
            [-0.5961 + 0.4214j, -0.0840 + 0.0085j],
            [-0.0840 + 0.0085j, -0.4339 + 0.4051j]], dtype=np.cdouble);
        SAB_actual[0,1] = np.array([
            [0.6020 - 0.3046j, -0.0431 + 0.0534j],
            [-0.0431 + 0.0534j, 0.6852 - 0.4078j]], dtype=np.cdouble);
        SAB_actual[1,0] = np.array([
            [0.6020 - 0.3046j, -0.0431 + 0.0534j],
            [-0.0431 + 0.0534j, 0.6852 - 0.4078j]], dtype=np.cdouble);
        SAB_actual[1,1] = np.array([
            [0.6971 - 0.2216j, 0.0672 - 0.0211j],
            [0.0672 - 0.0211j, 0.5673 - 0.1808j]], dtype=np.cdouble);

        np.testing.assert_allclose(SAB_actual[1,0], SAB_calc[1,0], rtol=reltol, atol=abstol);

    def testRedhefferProductS22(self):
        abstol = 0.0001;# Absolute error tolerance for test data (we only have it to 4 digits)
        reltol = 0.001; # Relative error tolerance (probably not necessary)

        # LAYER 1 DATA - This is the data 
        # Current global data. Before applying the Redheffer star product to
        # any values, S11/S22 should be zero and S12/S21 should be the identity.
        SA = np.zeros((2,2,2,2));
        SA[0,1] = np.identity(2);
        SA[1,0] = np.identity(2);

        SB = np.zeros((2,2,2,2), dtype=np.cdouble);
        SB[0,0] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.0060j, 0.0808 - 0.0121j]], dtype=np.cdouble);
        SB[0,1] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);
        SB[1,0] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);
        SB[1,1] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.0060j, 0.0808 - 0.0121j]], dtype=np.cdouble);
        D = np.array([
            [1,0],
            [0,1]], dtype=np.cdouble);
        F = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype = np.cdouble);

        SAB_actual = np.zeros((2,2,2,2), dtype=np.cdouble);
        SAB_actual[0,0] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)
        SAB_actual[0,1] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SAB_actual[1,0] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SAB_actual[1,1] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)
        SAB_calc = redhefferProduct(SA, SB);
        np.testing.assert_allclose(SAB_actual[1,1], SAB_calc[1,1], rtol=reltol, atol=abstol);

        # LAYER 2 DATA
        # Current global data. After applying the Redheffer star product once, we will end up
        # with a new set of values for the scattering matrix, which are given in the text.
        # They are the SG11/SG12/SG21/SG22 from the first layer (after the update)
        SA = np.zeros((2,2,2,2), dtype=np.cdouble);
        SA[0,0] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)
        SA[0,1] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SA[1,0] = np.array([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]], dtype=np.cdouble)
        SA[1,1] = np.array([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.006j, 0.0808 - 0.0121j]], dtype=np.cdouble)

        # The S11/etc. for the second layer.
        SB = np.zeros((2,2,2,2), dtype=np.cdouble);
        SB[0,0] = np.array([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]], dtype = np.cdouble);
        SB[0,1] = np.array([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]], dtype=np.cdouble);
        SB[1,0] = np.array([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]], dtype=np.cdouble);
        SB[1,1] = np.array([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]], dtype = np.cdouble);
        D = np.array([
            [0.1506 + 0.9886j, -0.0163 - 0.0190j],
            [-0.0163 - 0.0190j, 0.1822 + 1.0253j]], dtype=np.cdouble);
        F = np.array([
            [-0.2117 - 0.6413j, 0.0471 + 0.0518j],
            [0.0471 + 0.0518j, -0.3027 - 0.7414j]], dtype=np.cdouble);

        SAB_calc = redhefferProduct(SA, SB);
        # THIS STILL NEEDS TO BE FILLED IN. IT'S THE LAST GLOBAL S-MATRIX.
        SAB_actual = np.zeros((2,2,2,2), dtype=np.cdouble);
        SAB_actual[0,0] = np.array([
            [-0.5961 + 0.4214j, -0.0840 + 0.0085j],
            [-0.0840 + 0.0085j, -0.4339 + 0.4051j]], dtype=np.cdouble);
        SAB_actual[0,1] = np.array([
            [0.6020 - 0.3046j, -0.0431 + 0.0534j],
            [-0.0431 + 0.0534j, 0.6852 - 0.4078j]], dtype=np.cdouble);
        SAB_actual[1,0] = np.array([
            [0.6020 - 0.3046j, -0.0431 + 0.0534j],
            [-0.0431 + 0.0534j, 0.6852 - 0.4078j]], dtype=np.cdouble);
        SAB_actual[1,1] = np.array([
            [0.6971 - 0.2216j, 0.0672 - 0.0211j],
            [0.0672 - 0.0211j, 0.5673 - 0.1808j]], dtype=np.cdouble);

        np.testing.assert_allclose(SAB_actual[1,1], SAB_calc[1,1], rtol=reltol, atol=abstol);

    def testSrefFull(self):
        # I generated this Aref myself.
        # Tolerances relaxed due to small elements in the matrix
        abstol = 0.007;
        reltol = 0.03;
        Aref = np.array([
            [1.86002, 0.113614],
            [0.115376, 1.64547]], dtype=np.cdouble);
        Bref = np.array([
            [0.139976, -0.113614],
            [-0.115376, 0.354529]], dtype=np.cdouble);
        Sref_calc = Sref_gen(Aref, Bref);

        Sref_actual = np.zeros((2,2,2,2), dtype=np.cdouble);
        Sref_actual[0,0] = np.array([
            [-0.0800, 0.0761],
            [0.0761, -0.2269]], dtype=np.cdouble);
        Sref_actual[0,1] = np.array([
            [1.0800, -0.0761],
            [-0.0761, 1.2269]], dtype=np.cdouble);
        Sref_actual[1,0] = np.array([
            [0.9200, 0.0761],
            [0.0761, 0.7731]], dtype=np.cdouble);
        Sref_actual[1,1] = np.array([
            [0.0800, -0.0761],
            [-0.0761, 0.2269]], dtype=np.cdouble);
        np.testing.assert_allclose(Sref_calc, Sref_actual, atol=abstol, rtol=reltol);

    def testStrnFull(self):
        """
        WARNING: THIS HAS TO BE MODIFIED AND SHOULD CURRENTLY BE AN INTEGRATION TEST. I DO NOT HAVE
        RAW INPUT DATA FOR THIS FUNCTION. I AM RELYING ON AIJ/BIJ, and VWX to generate my input.
        I STILL HAVE TO WRITE THIS DAMN TEST.
        """
        abstol = 0.0001;
        reltol = 0.001;

        # I generated these myself from known-working Aij, Bij, Wtrn. They are now hard-coded.
        Atrn = np.array([
            [1.660774, -0.0652574],
            [-0.06525902, 1.786816]], dtype=np.cdouble);
        Btrn = np.array([
            [0.339226, 0.0652574],
            [0.06525902, 0.21318382]], dtype=np.cdouble);

        Strn_calc = Strn_gen(Atrn, Btrn);

        Strn_actual = np.zeros((2,2,2,2), dtype=np.cdouble);
        Strn_actual[0,0] = np.array([
            [0.2060, 0.0440],
            [0.0440, 0.1209]], dtype=np.cdouble);
        Strn_actual[0,1] = np.array([
            [0.7940, -0.0440],
            [-0.0440, 0.8791]], dtype=np.cdouble);
        Strn_actual[1,0] = np.array([
            [1.2060, 0.0440],
            [0.0440, 1.1209]], dtype=np.cdouble);
        Strn_actual[1,1] = np.array([
            [-0.2060, -0.0440],
            [-0.0440, -0.1209]], dtype=np.cdouble);
        np.testing.assert_allclose(Strn_actual, Strn_calc, atol=abstol, rtol=reltol);

    def SglobalIntegration(self):
        """
        Tests that the global s-matrices obtained for a given input angle and material properties
        are correct. Using data from Rumpf, in particular, the final global scattering matrices.
        This is the big daddy integration test. It requires calculation of all intermediate matrices,
        s parameters, and computation of redheffer star products. If this is correct, all our code is
        correct.
        """
        abstol=0.0001;
        reltol=0.001;
        l0 = 2.7;
        k0 = 2*np.pi / l0;
        theta_deg = 57.0;
        phi_deg = 23.0;
        theta = theta_deg * np.pi / 180.0;
        phi = phi_deg * np.pi / 180.0;
        pTE = 1/sqrt(2); # Our input polarization is circular
        pTM = 1j/sqrt(2);

        urref = 1.2;
        erref = 1.4;
        urtrn = 1.6;
        ertrn = 1.8;

        uri = [1.0, 3.0];
        eri = [2.0, 1.0]; # The inner layers
        Li = [0.25*l0, 0.5*l0];

        Sglobal_actual = np.zeros((2,2,2,2), dtype=np.cdouble);
        Sglobal_actual[0,0] = np.array([
            [-0.6018 + 0.3062j, -0.0043 + 0.0199j],
            [-0.0043 + 0.0199j, -0.5935 + 0.2678j]], dtype=np.cdouble);
        Sglobal_actual[0,1] = np.array([
            [0.5766 - 0.3110j, -0.0919 + 0.0469j],
            [-0.0919 + 0.0469j, 0.7542 - 0.4016j]], dtype=np.cdouble);
        Sglobal_actual[1,0] = np.array([
            [0.7415 - 0.4007j, 0.0716 - 0.0409j],
            [0.0716 - 0.0409j, 0.6033 - 0.3218j]], dtype=np.cdouble);
        Sglobal_actual[1,1] = np.array([
            [0.5861 - 0.3354j, 0.0170 + 0.0042j],
            [0.0170 + 0.0042j, 0.5533 - 0.3434j]], dtype=np.cdouble);

        np.testing.assert_all_close(Sglobal_actual, Sglobal_calc, atol=abstol, rtol=reltol);

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
