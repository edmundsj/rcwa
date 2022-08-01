import pytest
import numpy as np
from rcwa import Results
from numpy.testing import assert_allclose, assert_equal
from rcwa.examples.bragg_mirror import solve_system as solve_bragg
from rcwa.examples.diffraction_grating_1D import solve_system as solve_rectangular_1d
from rcwa.examples.diffraction_grating_triangular_1D import solve_system as solve_triangular_1d
from rcwa.examples.si_dispersive import solve_system as solve_si_dispersive
from rcwa.examples.SiO2_dispersive import solve_system as solve_siO2_dispersive
from rcwa.examples.si_ellipsometry import solve_system as solve_si_ellipsometry
from rcwa.examples.thin_film_dispersive import solve_system as solve_thin_film_dispersive
from rcwa.examples.triangular_photonic_crystal_2D import solve_system as solve_triangular_photonic_crystal
from rcwa.examples.grating_sweep import solve_system as solve_grating_sweep
from rcwa.examples.wavelength_angle_sweep import solve_system as wavelength_angle_sweep


@pytest.mark.example
def test_bragg_mirror():
    results = solve_bragg()
    assert isinstance(results, Results)
    assert 'wavelength' in results.keys()
    assert 'RTot' in results.keys()
    assert 'TTot' in results.keys()
    assert len(results['wavelength']) == 321
    RTot_desired_subset = np.array([
        0.3196300529185759, 0.2914974105926992, 0.22330554235257308,
        0.13178599256399812, 0.054328653032396854, 0.03580068889990351])
    RTot_actual_subset = results['RTot'][0:len(RTot_desired_subset)]
    assert_allclose(RTot_actual_subset, RTot_desired_subset, atol=1e-7)


@pytest.mark.example
def test_rectangular_1d():
    results = solve_rectangular_1d()
    assert 'rx' in results.keys()
    rx_actual =  results['rx']
    rx_desired = np.array([-0.01165311-0.00861393j, -0.03160784-0.01872845j,
       -0.06969903-0.00039625j, -0.07897799+0.06488948j,
        0.02983048+0.10016968j, -0.01012311-0.04538567j,
        0.01214772+0.09141299j,  0.08429551-0.00862236j,
        0.06702831-0.02112818j,  0.03356285-0.02818963j,
        0.01247477-0.01675098j])
    assert_allclose(rx_actual, rx_desired, atol=1e-7)


@pytest.mark.example
def test_triangular_1d():
    results = solve_triangular_1d()
    assert 'rx' in results.keys()
    rx_actual = results['rx']
    rx_desired = np.array([
        -1.12636392e-06+2.34933208e-07j, -6.91060869e-06+4.53669748e-06j,
        -5.59890347e-05+6.17713486e-05j, -5.80226897e-04+6.80539196e-04j,
        -6.53509916e-03+7.24777873e-03j,  3.52339625e-02+7.19222657e-02j,
        3.08931024e-02-3.50408518e-02j,  3.59406795e-03-1.56088999e-03j,
        3.52427618e-04-5.51611054e-05j,  2.71156851e-05-1.96513775e-06j,
        1.42429830e-06-1.31545273e-06j])
    assert_allclose(rx_actual, rx_desired, atol=1e-7)


@pytest.mark.example
def test_solve_si_dispersive():
    results = solve_si_dispersive()
    assert 'wavelength' in results.keys()

    wavelengths_subset_desired = np.array([0.25 , 0.251, 0.252, 0.253, 0.254, 0.255, 0.256, 0.257, 0.258])
    wavelengths_subset_actual = results['wavelength'][0:len(wavelengths_subset_desired)]
    assert_equal(wavelengths_subset_actual, wavelengths_subset_desired)


@pytest.mark.example
def test_solve_si_ellipsometry():
    with pytest.warns(Warning):
        fig, ax = solve_si_ellipsometry()
        assert fig is not None
        assert ax is not None


@pytest.mark.example
def test_solve_siO2_dispersive():
    results = solve_siO2_dispersive()
    assert 'wavelength' in results.keys()

    wavelengths_subset_desired = np.array([0.25 , 0.251, 0.252, 0.253, 0.254, 0.255, 0.256, 0.257, 0.258])
    wavelengths_subset_actual = results['wavelength'][0:len(wavelengths_subset_desired)]
    assert_equal(wavelengths_subset_actual, wavelengths_subset_desired)


@pytest.mark.example
def test_solve_thin_film_dispersive():
    with pytest.warns(Warning):
        results = solve_thin_film_dispersive()
        assert 'wavelength' in results.keys()

        wavelengths_subset_desired = np.array([0.25 , 0.251, 0.252, 0.253, 0.254, 0.255, 0.256, 0.257, 0.258])
        wavelengths_subset_actual = results['wavelength'][0:len(wavelengths_subset_desired)]
        assert_equal(wavelengths_subset_actual, wavelengths_subset_desired)


@pytest.mark.example
def test_solve_triangular_pc():
    results = solve_triangular_photonic_crystal()
    for x in ['conservation', 'rx', 'RTot']:
        assert x in results.keys()
    rx_desired = np.array([ 0.0194899 +0.01175673j, -0.05570714+0.03010019j,
       -0.03797483-0.03124727j, -0.03920141+0.01425774j,
        0.20357306-0.00134404j, -0.01196725-0.02682893j,
       -0.00221492-0.00226153j, -0.0298306 -0.00790429j,
        0.0283686 +0.01849218j])

    assert_allclose(results['conservation'], 1)
    assert_allclose(results['RTot'], 0.15593548090823117)
    assert_allclose(results['TTot'], 0.844064519091769)
    assert_allclose(results['rx'], rx_desired, atol=1e-8, rtol=1e-5)


@pytest.mark.example
def test_grating_sweep():
    results = solve_grating_sweep()
    assert 'thickness' in results.keys()
    assert len(np.unique(results['RTot'])) == len(results['thickness'])


@pytest.mark.example
def test_wavelength_angle():
    results = wavelength_angle_sweep()
    assert 'theta' in results.keys()
    assert 'wavelength' in results.keys()
    assert len(np.unique(results['RTot'])) == len(results['theta'])

