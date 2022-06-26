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
    _, sweep_vars, sweep_vals = solver._sweeps(wavelength=[1, 2, 3])
    desired_vars = ['wavelength']
    desired_vals = [(1,), (2,), (3,)]
    assert_equal(desired_vars, list(sweep_vars))
    assert_equal(desired_vals, list(sweep_vals))

def test_sweeps_kw_2D(solver):
    _, sweep_vars, sweep_vals = solver._sweeps(wavelength=[1, 2], phi=[0, 1, 2])
    desired_vars = ['wavelength', 'phi']
    desired_vals = ((1, 0), (1, 1), (1, 2), (2,0), (2, 1), (2, 2))
    assert_equal(list(sweep_vars), desired_vars)
    assert_equal(list(sweep_vals), desired_vals)

def test_setup_sweep_arg_1d(solver_2):
    layer, source, solver = solver_2
    thickness = [0.1, 0.2, 0.3]
    desired_sweeps = [(0.1,), (0.2,), (0.3,)]
    sweep_objs, sweep_vars, sweep_vals = solver._sweeps((layer, {'thickness': thickness}))
    assert sweep_objs[0] is layer
    assert sweep_vars[0] == 'thickness'
    assert_equal(sweep_vals, desired_sweeps)

def test_setup_sweep_arg_kw_2d(solver_2):
    layer, source, solver = solver_2
    thickness = [0.1, 0.2, 0.3]
    desired_sweeps = [(0.1,1), (0.1, 2), (0.2,1), (0.2, 2), (0.3,1), (0.3, 2)]
    desired_vars = ['thickness', 'wavelength']
    desired_objs = [layer, None]
    sweep_objs, sweep_vars, sweep_vals = solver._sweeps((layer, {'thickness': thickness}), wavelength=[1, 2])
    assert sweep_objs[0] is desired_objs[0]
    assert sweep_objs[1] is desired_objs[1]
    assert_equal(sweep_vars, desired_vars)
    assert_equal(sweep_vals, desired_sweeps)

def test_assign_sweep_vars(solver_2):
    layer, source, solver = solver_2
    thickness = [0.1, 0.2, 0.3]
    solver.sweep_vars = ['thickness']
    solver.sweep_objects = [layer]
    solver.sweep_vals = [(0.1,), (0.2,), (0.3,)]
    for sweep in solver.sweep_vals:
        solver._assign_sweep_vars(sweep)
        assert layer.thickness == sweep[0]
