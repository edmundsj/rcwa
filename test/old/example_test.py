# This is our test case, where we test the sum() function for an arbitrary set of inputs.
import sys
sys.path.append("..");

def test_sum_list(x, y):
    assert sum([x, y]) == x + y, f"Sum error is {sum([x,y])-(x+y)}";

def test_sum_tuple(x, y):
    assert sum((x,y,2)) == x + y, f"Sum error is {sum((x,y,2))-(x+y)}";
# What's the name == main doing here?
if __name__ == "__main__":
    test_sum_list(1, 2); # This is the line where we actually test something.
    test_sum_tuple(1,3);
    print("Everything passed");
