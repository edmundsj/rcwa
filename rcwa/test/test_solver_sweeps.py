import pytest
from rcwa import Solver, Layer, LayerStack, Source
from numpy.testing import assert_equal

@pytest.fixture
def solver():
    sol = Solver(LayerStack(), Source())
    yield sol

@pytest.fixture
def layer():
    layer = Layer(thickness=0.1)
    yield layer

@pytest.fixture
def solver_2(layer):
    stack = LayerStack(layer)
    source = Source()
    sol = Solver(stack, source)
    yield layer, source, sol

def test_sweeps_kw_1D(solver):
    sweep_vars, sweep_vals = solver._sweeps(wavelength=[1, 2, 3])
    desired_vars = ['wavelength']
    desired_vals = [(1,), (2,), (3,)]
    assert_equal(desired_vars, list(sweep_vars))
    assert_equal(desired_vals, list(sweep_vals))

def test_sweeps_kw_2D(solver):
    sweep_vars, sweep_vals = solver._sweeps(wavelength=[1, 2], phi=[0, 1, 2])
    desired_vars = ['wavelength', 'phi']
    desired_vals = ((1, 0), (1, 1), (1, 2), (2,0), (2, 1), (2, 2))
    assert_equal(list(sweep_vars), desired_vars)
    assert_equal(list(sweep_vals), desired_vals)

def test_sweep_args(solver_2):
    layer, source, solver = solver_2

def test_val_not_found_sweep(solver):
    with pytest.raises(ValueError) as e:
        solver._assign_sweep_vars(['nonsense'], (1,))

def test_assign_sweep_vars(solver):
    sweep_vars, sweep_vals = ['wavelength'], [(1,),(2,),(3,)]
    solver._assign_sweep_vars(sweep_vars, sweep_vals[1])
    assert_equal(solver.source.wavelength, 2)

