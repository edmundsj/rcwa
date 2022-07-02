import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_almost_equal
from rcwa import Solver, LayerStack, Layer, Source
from rcwa.shorthand import complexArray
from rcwa.utils import rTE, k_vector

@pytest.fixture
def source():
    src = Source(wavelength=1, theta=0, phi=0, pTEM=[1, 0])
    yield src

@pytest.fixture
def stack_interface():
    """
    Layer stack for an interface
    """
    incident_layer = Layer(er=5, ur=1)
    transmission_layer = Layer(er=1, ur=7)
    stack = LayerStack(incident_layer=incident_layer, transmission_layer=transmission_layer)
    yield stack

@pytest.fixture
def solver_interface(stack_interface, source):
    sol = Solver(layer_stack=stack_interface, source=source)
    sol.solve()
    yield sol

def test_E_point(solver_interface):
    stack = solver_interface.layer_stack
    source = solver_interface.source
    z0 = -0.0

    rte = rTE(source, stack.incident_layer, stack.transmission_layer)

    # NEED TO CONVERT FROM ANGLE OF SOURCE TO POLARIZATION ANGLE - I DONT TRUST THIS CURRENTLY.
    # WHY DOES PHI=0 CORRESPOND TO Y-POLARIZATION? I WOULD EXPECT IT TO CORRESPOND TO X.

    field_desired_x = (1 + rte) * np.cos(source.theta) * np.sin(source.phi)
    field_desired_y = (1 + rte) * np.cos(source.theta) * np.cos(source.phi)
    field_desired = complexArray([field_desired_x, field_desired_y])

    field_actual = solver_interface.fields(layer=stack.incident_layer, component='E', z_min=z0, z_max=z0)

    assert_almost_equal(field_actual, field_desired)
