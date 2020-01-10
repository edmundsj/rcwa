# I decided to write my own unit tests rather than use python's unit testing framework because
# it was causing me more trouble than it was worth, since it doesn't have a built-in capability
# to use multiple datasets to run the same test. It's extremely annoying.

import sys
sys.path.append('core');

from matrices import *

def testCaller(testFunction, *args):
    test_status = False; # By default assume we failed the test.
    test_msg = f"{testFunction.__name__}({args}): ";

    try:
        print(f"Calling function {testFunction.__name__} ... ", end=" ");
        testFunction(*args);
        print("OK");
        test_status = True;
    except AssertionError as ae:
        print("FAIL");
        test_msg += "FAILED";
        test_msg += str(ae);
    except Exception as e:
        print(e)
    finally:
        return (test_status, test_msg);



def test_PMatrix(kx, ky, er, ur):
    actual_val = 1/er *  np.array([
        [kx * ky, ur * er - sq(kx)],
        [sq(ky) - ur * er, - kx * ky]]);

    P_val = Pi_gen(kx, ky, er, ur);
    np.testing.assert_array_equal(P_val, actual_val);

def test_QMatrix(kx, ky, er, ur):
    actual_val = 1 / ur * np.array([
        [kx * ky, ur * er - sq(kx)],
        [sq(ky) - er * ur, - kx * ky]]);


    Q_val = Qi_gen(kx, ky, er, ur);
    np.testing.assert_array_equal(Q_val, actual_val);

def main():
    print("--------- RUNNING UNIT TESTS... ----------");
    NUM_PARAMS = 5;
    theta_test = np.linspace(0, np.pi/2, NUM_PARAMS);
    phi_test = np.linspace(0, np.pi*2, NUM_PARAMS);
    er_test = np.linspace(1, 10, NUM_PARAMS);
    ur_test = np.linspace(1, 10, NUM_PARAMS);
    statuses = [];
    messages = [];

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

                    (test_status, test_message) = testCaller(test_PMatrix, kx, ky, er, ur);
                    statuses.append(test_status);
                    messages.append(test_message);

                    (test_status, test_message) = testCaller(test_QMatrix, kx, ky, er, ur);
                    statuses.append(test_status);
                    messages.append(test_message);


    print("----------- DONE ------------");
    print(f"{statuses.count(True)} passed, {statuses.count(False)} failed")

    for s, i in zip(statuses, range(len(statuses))):
        if(s == False):
            print(messages[i]);

main();
