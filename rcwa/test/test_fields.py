import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal
from rcwa import Solver, LayerStack, Layer, Source, k_vector
from rcwa.utils import rTE

@pytest.fixture
def source():
    src = Source(wavelength=1, theta=np.pi/6, phi=np.pi/4, pTEM=[1, 0])
    yield src

@pytest.fixture
def stack_interface():
    """
    Layer stack for an interface
    """
    incident_layer = Layer(er=2, ur=1)
    transmission_layer = Layer(er=1, ur=2)
    stack = LayerStack(incident_layer=incident_layer, transmission_layer=transmission_layer)
    yield stack

@pytest.fixture
def solver_interface(stack_interface, source):
    sol = Solver(layer_stack=stack_interface, source=source)
    sol.solve()
    yield sol

@pytest.mark.skip
def test_Ex_point(solver_interface):
    stack = solver_interface.layer_stack
    source = solver_interface.source
    er1, er2 = stack.incident_layer.er, stack.transmission_layer.er
    ur1, ur2 = stack.incident_layer.ur, stack.transmission_layer.ur
    k_inc = k_vector(source, stack.incident_layer)
    k_tran = k_vector(source, stack.transmission_layer)
    z0 = -0.1
    field_actual = solver_interface.fields(component='Ex', z_min=z0, z_max=z0)
    field_desired = 1 + rTE(k_inc[2], k_tran[2], er1, er2, ur1, ur2)
    breakpoint()
    assert_almost_equal(field_actual, field_desired)

