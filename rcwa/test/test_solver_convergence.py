import pytest
from rcwa import Solver, Layer, LayerStack, Source
from numpy.testing import assert_equal

@pytest.fixture
def solver():
    sol = Solver(LayerStack(), Source())
    yield sol

@pytest.mark.parametrize(
    'factor,input,desired',
    [(1, 1, 3), (1, 3, 5),
     (1, (1, 1), (3, 3)), (1, (3, 3), (5, 5))])
def test_increase_harmonics(solver, factor, input, desired):
    solver.n_harmonics = input
    desired_harmonic = desired
    solver._increase_harmonics(factor=factor)
    actual_harmonics = solver.n_harmonics
    assert_equal(actual_harmonics, desired_harmonic)