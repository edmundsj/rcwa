# I decided to write my own unit tests rather than use python's unit testing framework because
# it was causing me more trouble than it was worth, since it doesn't have a built-in capability
# to use multiple datasets to run the same test. It's extremely annoying.

import sys
sys.path.append('core');

from matrices import *
from fresnel import *

def main():
    print("--------- RUNNING UNIT TESTS... ----------");
    NUM_PARAMS = 5;
    theta_test = np.linspace(0, np.pi/2, NUM_PARAMS);
    phi_test = np.linspace(0, np.pi*2, NUM_PARAMS);
    er_test = np.linspace(1, 10, NUM_PARAMS);
    ur_test = np.linspace(1, 10, NUM_PARAMS);
    statuses = [];
    messages = [];
    unit_tests_enabled = False;
    integration_tests_enabled = True;

    if(unit_tests_enabled == True):
        for theta, i in zip(theta_test, range(len(theta_test))):
            for phi in phi_test:
                for er in er_test:
                    for ur in ur_test:
                        # FIRST, RUN ALL OF OUR UNIT TESTS. ESP. SCATTERING MATRIX..
                        print(f"\nUsing parameters theta={theta}")
                        # First, setup our local dependent variables.
                        n = sqrt(ur * er);
                        kx = n * sin(theta) * cos(phi);
                        ky = n * sin(theta) * sin(phi);

                        (test_status, test_message) = testCaller(testPMatrix, kx, ky, er, ur);
                        statuses.append(test_status);
                        messages.append(test_message);

                        (test_status, test_message) = testCaller(testQMatrix, kx, ky, er, ur);
                        statuses.append(test_status);
                        messages.append(test_message);


    # NOW, RUN OUR INTEGRATION TESTS.
    # First integration test: making sure reflections at an interface between
    # two materials are consistent with Fresnel's equation, and check for energy
    # conservation 

    if(integration_tests_enabled == True):
        # Run our integration tests. Start with the S-parameters.
        # There is a known problem with these tests when we are beyond the critical angle.
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

        # NOW, RUN OUR INTEGRATION TESTS.
        # First integration test: making sure reflections at an interface between
        # two materials are consistent with Fresnel's equation, and check for energy
        # conservation 
        (test_status, test_message) = testCaller(testS11Ref, kx_n, ky_n, er, ur, ATEM);
        statuses.append(test_status);
        messages.append(test_message);

        (test_status, test_message) = testCaller(testS22Ref, kx_n, ky_n, er, ur, ATEM);
        statuses.append(test_status);
        messages.append(test_message);

        (test_status, test_message) = testCaller(testS11Trn, kx_n, ky_n, er, ur, ATEM);
        statuses.append(test_status);
        messages.append(test_message);

        (test_status, test_message) = testCaller(testS22Trn, kx_n, ky_n, er, ur, ATEM);
        statuses.append(test_status);
        messages.append(test_message);

        (test_status, test_message) = testCaller(testS21Ref, kx_n, ky_n, er, ur, ATEM);
        statuses.append(test_status);
        messages.append(test_message);

        (test_status, test_message) = testCaller(testS21Trn, kx_n, ky_n, er, ur, ATEM);
        statuses.append(test_status);
        messages.append(test_message);

    print("----------- DONE ------------");
    print(f"{statuses.count(True)} passed, {statuses.count(False)} failed")

    for s, i in zip(statuses, range(len(statuses))):
        if(s == False):
            print(messages[i]);

def testCaller(testFunction, *args):
    """
    This function is the wrapper for all the functions we want to test. By passing in all the
    function arguments and the test function itself, we generate test reports and handle
    assertion error exceptions.
    """
    test_status = False; # By default assume we failed the test.
    test_msg = f"{testFunction.__name__}({args}): ";

    try:
        print(f"Calling function {testFunction.__name__} ... ", end=" ");
        testFunction(*args);
        print("OK");
        test_status = True;
        return (test_status, test_msg);
    except AssertionError as ae:
        print("FAIL");
        test_msg += "FAILED";
        test_msg += str(ae);
        return (test_status, test_msg);

def testPMatrix(kx, ky, er, ur):
    actual_val = 1/er *  np.array([
        [kx * ky, ur * er - sq(kx)],
        [sq(ky) - ur * er, - kx * ky]]);

    P_val = Pi_gen(kx, ky, er, ur);
    np.testing.assert_array_equal(P_val, actual_val);

def testQMatrix(kx, ky, er, ur):
    actual_val = 1 / ur * np.array([
        [kx * ky, ur * er - sq(kx)],
        [sq(ky) - er * ur, - kx * ky]]);


    Q_val = Qi_gen(kx, ky, er, ur);
    np.testing.assert_array_equal(Q_val, actual_val);

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

def testRTInterface():
    """

    """
    print("Function not implemented");


def testRTEtalon():
    # Test the reflectivity and transmission of an etalon using the equations we can derive directly
    # from Fresnel's equations. This will be a little tricky as I may have to derive them again for
    # permittivities and permeabilities not equal to that of free space.
    print("Not implemented yet");

main();
