from rcwa.utils import rTE, rTM, k_vector
from rcwa.shorthand import complexArray
from rcwa import Source, Layer, LayerStack
import pytest
from numpy.testing import assert_equal, assert_allclose
import numpy as np

@pytest.fixture
def normal_source():
    source = Source(phi=0, theta=0)
    yield source

@pytest.fixture
def transparent_stack():
    incident_layer = Layer(er=2, ur=3)
    transmission_layer = Layer(er=2, ur=3)
    stack = LayerStack(incident_layer=incident_layer, transmission_layer=transmission_layer)
    yield stack

@pytest.fixture
def opaque_stack_n():
    incident_layer = Layer(er=1, ur=1)
    transmission_layer = Layer(er=1e8, ur=1)
    stack = LayerStack(incident_layer=incident_layer, transmission_layer=transmission_layer)
    yield stack

@pytest.fixture
def opaque_stack_p():
    incident_layer = Layer(er=1e8, ur=1)
    transmission_layer = Layer(er=1, ur=1)
    stack = LayerStack(incident_layer=incident_layer, transmission_layer=transmission_layer)
    yield stack

def test_rTE_transparent(normal_source, transparent_stack):
    desired_rTE = 0
    actual_rTE = rTE(normal_source, transparent_stack.incident_layer, transparent_stack.transmission_layer)
    assert_equal(actual_rTE, desired_rTE)

def test_rTE_opaque_neg(normal_source, opaque_stack_n):
    stack = opaque_stack_n
    desired_rTE = -1
    actual_rTE = rTE(normal_source, stack.incident_layer, stack.transmission_layer)
    assert_allclose(actual_rTE, desired_rTE, atol=2e-4)

def test_rTE_opaque_pos(normal_source, opaque_stack_p):
    stack = opaque_stack_p
    desired_rTE = 1
    actual_rTE = rTE(normal_source, stack.incident_layer, stack.transmission_layer)
    assert_allclose(actual_rTE, desired_rTE, atol=2e-4)

@pytest.mark.parametrize('theta,phi,desired_vec', [
    (0, 0, [0, 0, 1]),
    (0, np.pi/2, [0, 0, 1]),
    (np.pi/2, 0, [1, 0, 0]),
    (-np.pi/2, 0, [-1, 0, 0]),
    (np.pi/2, np.pi/2, [0, 1, 0])
])
def test_k_vector_varying_angles(theta, phi, desired_vec):
    layer = Layer()
    wavelength = 1
    source = Source(theta=theta, phi=phi, wavelength=wavelength)
    desired_k_vector = complexArray(desired_vec) * 2 * np.pi / wavelength
    actual_k_vector = k_vector(source, layer)
    assert_allclose(actual_k_vector, desired_k_vector, atol=1e-10)
