from rcwa.shorthand import *
from rcwa.matrices import *
from rcwa.harmonics import *
from rcwa import Layer, LayerStack
from copy import deepcopy
from progressbar import ProgressBar, Bar, Counter, ETA
from itertools import product


class Solver:
    """ Main class that invokes all methods necessary to solve an RCWA/TMM simulation

    :param layer_stack: layerStack: Stack of layers to simulate
    :param source: Source object which includes wavelength and direction information
    :param n_harmonics: The number of harmonics in x, y to simulate (number of Fourier components). For planar films this should be 1. For 1D diffraction gratings this should be a single integer. For 2D periodic films this should be a 2-tuple. Must be an odd number.
    """
    def __init__(self, layer_stack, source, n_harmonics=1):
        self.n_harmonics = n_harmonics
        self.layer_stack = layer_stack
        self.source = source
        self.source.layer = layer_stack.incident_layer
        self.layer_stack.source = source

        self.layer_stack._set_convolution_matrices(n_harmonics)
        self.baseCrystal = self.layer_stack.crystal
        self._k_matrices()
        self._gap_matrices()
        self._outer_matrices()
        self.results = []

    def solve(self, *sweep_args, **sweep_kw):
        """
        Solves the simulation or performs a simulation sweep of the desired parameters

        :param sweep_args: Objects along with their parameters to sweep. i.e. (layer, {'thickness': [1, 2, 2.5]})
        :param sweep_kw: Source variables to sweep (theta, phi, wavelength, etc.). Can either be a single keyword argument or several. If several are used, all combinations of the two parameters will be made
        """
        self.results = []
        self.sweep_vars, self.sweep_vals = self._sweeps(*sweep_args, **sweep_kw)
        n_sweeps = len(self.sweep_vals)

        bar = ProgressBar(widgets=[Counter(), f'/{n_sweeps} ', Bar(), ETA()], max_value=n_sweeps).start()

        for i, sweep in enumerate(self.sweep_vals):
            self._assign_sweep_vars(self.sweep_vars, sweep)
            self._couple_source()

            self._k_matrices()
            self._gap_matrices()
            self._outer_matrices()
            self._initialize()
            self._inner_s_matrix()
            self._global_s_matrix()
            self._rt_quantities()
            self._append_results()
            bar.update(i)

        bar.finish()
        self.results = self._package_results()

    @staticmethod
    def _sweeps(**sweep_kw):
        sweep_vars = sweep_kw.keys()
        sweep_vals = list(product(*sweep_kw.values()))
        return sweep_vars, sweep_vals

    def _couple_source(self):
        self.source.layer = self.layer_stack.incident_layer

    def _assign_sweep_vars(self, sweep_vars, sweep):
        for var, val in zip(sweep_vars, sweep):
            if not hasattr(self.source, var):
                raise ValueError(f'Source does not have attribute {var}. Invalid sweep variable. Available variables are "phi", "theta", "wavelength", "pTEM"')
            setattr(self.source, var, val)

    def _rt_quantities(self):
        self.rx, self.ry, self.rz = calculateReflectionCoefficient(self.SGlobal, self.Kx, self.Ky,
                                                                   self.KzReflectionRegion, self.WReflectionRegion, self.source, self.n_harmonics)
        self.tx, self.ty, self.tz = calculateTransmissionCoefficient(self.SGlobal, self.Kx, self.Ky,
                                                                     self.KzTransmissionRegion, self.WTransmissionRegion, self.source, self.n_harmonics)
        self.R = calculateDiffractionReflectionEfficiency(self.rx, self.ry, self.rz, self.source,
                                                          self.KzReflectionRegion, self.layer_stack, self.n_harmonics)
        self.T = calculateDiffractionTransmissionEfficiency(self.tx, self.ty, self.tz, self.source,
                                                            self.KzTransmissionRegion, self.layer_stack, self.n_harmonics)
        self.RTot = np.sum(self.R)
        self.TTot = np.sum(self.T)
        self.conservation = self.RTot + self.TTot

        if self.TMMSimulation is True:
            self.rTEM = calculateTEMReflectionCoefficientsFromXYZ(self.source, self.rx, self.ry, self.rz)

    def _package_results(self):
        """
        Turns the list of simulation results created during simulation into something more useful
        """

        n_results = len(self.results)
        result_keys = self.results[0].keys()
        new_results = {}

        if n_results > 1:
            for key in result_keys:
                new_results[key] = []
                for result in self.results:
                    new_results[key].append(result[key])

            for i, key in enumerate(self.sweep_vars):
                new_results[key] = []
                for sweep in self.sweep_vals:
                    new_results[key].append(sweep[i])
        else:
            new_results = self.results[0]

        return new_results

    def _append_results(self):
        """
        Packages the results from the simulation into a dictionary
        """
        tempResults = {}
        tempResults['rx'], tempResults['ry'], tempResults['rz'] = deepcopy((self.rx, self.ry, self.rz))
        tempResults['tx'], tempResults['ty'], tempResults['tz'] = deepcopy((self.tx, self.ty, self.tz))
        tempResults['R'], tempResults['T'] = deepcopy((self.R, self.T))
        tempResults['RTot'], tempResults['TTot'], tempResults['conservation'] = \
                deepcopy((self.RTot, self.TTot, self.conservation))
        tempResults['crystal'] = deepcopy(self.baseCrystal)
        tempResults['source'] = deepcopy(self.source)
        tempResults['S'] = deepcopy(self.SGlobal)

        if self.TMMSimulation is True:
            tempResults['rTE'] = self.rTEM[0]
            tempResults['rTM'] = self.rTEM[1]
            rho = tempResults['rTM'] / tempResults['rTE']
            tempResults['tanPsi'] = np.abs(rho)
            tempResults['cosDelta'] = cos(np.angle(rho))
            tempResults['delta'] = np.angle(rho)

        self.results.append(tempResults)

    def _k_matrices(self):
        """
        Sets up the Kx, Ky, and Kz matrices for solving the simulation once the source, crystal, and
        number harmonics are known.
        """
        self.Kx = kx_matrix(self.source, self.baseCrystal, self.n_harmonics)
        self.Ky = ky_matrix(self.source, self.baseCrystal, self.n_harmonics)
        if isinstance(self.Kx, np.ndarray):
            self.KDimension = self.Kx.shape[0]
            self.TMMSimulation = False
        else:
            self.KDimension = 1
            # Ensure that Kz for the gap layer is 1
            self.layer_stack.gapLayer = Layer(er=1 + sq(self.Kx) + sq(self.Ky), ur=1, thickness=0)
            self.TMMSimulation = True

        self.KzReflectionRegion = calculateKzBackward(self.Kx, self.Ky, self.layer_stack.incident_layer)
        self.KzTransmissionRegion = calculateKzForward(self.Kx, self.Ky, self.layer_stack.transmission_layer)
        self.KzGapRegion = calculateKzForward(self.Kx, self.Ky, self.layer_stack.gapLayer)

        self.scatteringElementDimension = self.KDimension * 2
        self.scatteringElementShape = (self.scatteringElementDimension, self.scatteringElementDimension)

    def _outer_matrices(self):
        self.WReflectionRegion = complexIdentity(self.scatteringElementDimension)
        self.WTransmissionRegion = complexIdentity(self.scatteringElementDimension)

    def _gap_matrices(self):
        self.WGap = complexIdentity(self.scatteringElementDimension)
        QGap = calculateQMatrix(self.Kx, self.Ky, self.layer_stack.gapLayer)
        LambdaGap = calculateLambdaMatrix(self.KzGapRegion)
        self.VGap = QGap @ inv(LambdaGap)

    def _inner_s_matrix(self):
        for i in range(len(self.layer_stack.internal_layers)):
            self.Si[i] = calculateInternalSMatrix(self.Kx, self.Ky, self.layer_stack.internal_layers[i],
                                                  self.source, self.WGap, self.VGap)
            self.SGlobal = calculateRedhefferProduct(self.SGlobal, self.Si[i])

    def _global_s_matrix(self):
        self.STransmission = calculateTransmissionRegionSMatrix(self.Kx, self.Ky, self.layer_stack,
                                                                self.WGap, self.VGap)
        self.SReflection = calculateReflectionRegionSMatrix(self.Kx, self.Ky, self.layer_stack,
                                                            self.WGap, self.VGap)
        self.SGlobal = calculateRedhefferProduct(self.SGlobal, self.STransmission)
        self.SGlobal = calculateRedhefferProduct(self.SReflection, self.SGlobal)

    def _initialize(self):
        self.SGlobal = generateTransparentSMatrix(self.scatteringElementShape)
        self.rx, self.ry, self.rz = None, None, None
        self.tx, self.ty, self.tz = None, None, None
        self.R, self.T, self.RTot, self.TTot, self.CTot = None, None, None, None, None
        self.Si = [None for _ in range(len(self.layer_stack.internal_layers))]
