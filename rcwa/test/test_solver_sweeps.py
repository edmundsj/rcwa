import pytest
from rcwa import Solver, LayerStack, Source
from numpy.testing import assert_equal

@pytest.fixture
def solver():
    stack = LayerStack()
    sol = Solver(LayerStack(), Source())
    yield sol

def test_sweeps_1D(solver):
    sweep_vars, sweep_vals = solver.sweeps(wavelength=[1, 2, 3])
    desired_vars = ['wavelength']
    desired_vals = [(1,), (2,), (3,)]
    assert_equal(desired_vars, list(sweep_vars))
    assert_equal(desired_vals, list(sweep_vals))

def test_val_not_found_sweep(solver):
    with pytest.raises(ValueError) as e:
        solver.assign_sweep_vars(['nonsense'], (1,))

def test_assign_sweep_vars(solver):
    sweep_vars, sweep_vals = ['wavelength'], [(1,),(2,),(3,)]
    solver.assign_sweep_vars(sweep_vars, sweep_vals[1])
    assert_equal(solver.source.wavelength, 2)
