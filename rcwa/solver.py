from rcwa.shorthand import *
from rcwa.matrices import *
from rcwa.harmonics import kx_matrix, ky_matrix
from rcwa import Layer, LayerStack, Results
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
        self.atol = None
        self.rtol = None
        self.max_iters = None
        self.check_convergence = None

        self.n_harmonics = n_harmonics
        self.layer_stack = layer_stack
        self.source = source
        self.source.layer = layer_stack.incident_layer
        self.layer_stack.source = source

        self._initialize()
        self._k_matrices()
        self._gap_matrices()
        self._outer_matrices()
        self.results = []

    def solve(self, *sweep_args, max_iters=50, atol=1e-3, rtol=1e-2, check_convergence=False, **sweep_kw):
        """
        Solves the simulation or performs a simulation sweep of the desired parameters

        :param sweep_args: Objects along with their parameters to sweep. i.e. (layer, {'thickness': [1, 2, 2.5]})
        :param sweep_kw: Source variables to sweep (theta, phi, wavelength, etc.). Can either be a single keyword argument or several. If several are used, all combinations of the two parameters will be made
        :param max_iters: Maximum number of iterations to complete before convergence
        :param atol: Absolute tolerance threshold for total reflectance at which simulation has converged
        :param rtol: Relative tolerance threshold for total reflectance at which simulation has converged
        :param check_convergence: If True, perform convergence testing for non-TMM simulations
        """
        self.atol = atol
        self.rtol = rtol
        self.max_iters = max_iters

        self.check_convergence = check_convergence
        self.converged = False
        self.iters = 0
        self.last_RTot = 1

        self.results = []
        self.sweep_objects, self.sweep_vars, self.sweep_vals = self._sweeps(*sweep_args, **sweep_kw)
        n_sweeps = len(self.sweep_vals)

        bar = ProgressBar(widgets=[Counter(), f'/{n_sweeps} ', Bar(), ETA()], max_value=n_sweeps).start()

        for i, sweep in enumerate(self.sweep_vals):

            self._assign_sweep_vars(sweep)

            while not self.converged:

                self._initialize()
                self._inner_s_matrix()
                self._global_s_matrix()
                self._rt_quantities()
                self.iters += 1
                self.converged = self._check_converged()

                if not self.converged:
                    self.last_RTot = self.RTot
                    self._increase_harmonics()

            self._append_results()
            self.iters = 0
            self.last_RTot = 1
            self.converged = False
            bar.update(i)

        bar.finish()
        self.results = self._package_results()
        return self.results

    def fields(self, component='Ex', layer=None, x_min=0, x_max=0, y_min=0, y_max=0, z_min=0, z_max=0, N_x=1, N_y=1, N_z=1):
        # First, we find the forward- and backward propagating waves in the incident region
        V_inc, W_inc, L_inc, _ = self.layer_stack.incident_layer.VWLX_matrices()
        c_incident = np.linalg.inv(W_inc) @ s_incident(self.source, self.n_harmonics)
        c_reflected = self.SGlobal[0, 0] @ c_incident

        if layer is self.layer_stack.incident_layer:
            c_forward_target = c_incident
            c_backward_target = c_reflected
        elif layer is self.layer_stack.transmission_layer:
            c_forward_target = self.SGlobal[1, 0] @ c_incident
            c_backward_target = 0 * c_incident
        else:
            raise NotImplementedError

        V_target, W_target, L_target, _ = self.layer_stack.incident_layer.VWLX_matrices()
        z = z_min

        if 'E' in component:
            field_target = W_target @ matrixExponentiate(-1 * L_target * self.source.k0 * z) @ c_forward_target + \
                           W_target @ matrixExponentiate(L_target * self.source.k0 * z) @ c_backward_target

        return field_target

    @property
    def base_crystal(self):
        return self.layer_stack.crystal

    def grad(self, loss_func, obj, attribute):
        """
        Computes the gradient of a user-specified loss function with respect to an attribute
        of an object in the simulation.

        :param loss_func: The loss function you are trying to optimize. Should take a single argument
        of the Solver Results object.
        :param obj: The object whose attribute you want to tweak in order to do the optimization (i.e. layer3)
        :param attribute: The attribute of the object you want to tweak (i.e. 'thickness' or 'er')
        """
        raise NotImplementedError


    def _increase_harmonics(self, factor=1):
        n_harmonics = np.array(self.n_harmonics)
        n_harmonics *= factor
        n_harmonics += 2
        even_elements = np.logical_not((n_harmonics % 2).astype(bool))
        n_harmonics[even_elements] -= 1
        if n_harmonics.size == 1:
            self.n_harmonics = int(n_harmonics)
        else:
            self.n_harmonics = tuple(n_harmonics)

    def _check_converged(self):
        converged = False
        if self.iters >= self.max_iters:
            raise RuntimeError('Exceeded maximum number of iterations {self.max_iters} without convergence. Aborting.')

        self.relative_error = np.abs((self.last_RTot - self.RTot)/self.last_RTot)
        self.absolute_error = np.abs(self.RTot - self.last_RTot)

        if self.TMMSimulation and self.iters > 0:
            converged = True

        if not self.check_convergence:
            converged = True

        if self.relative_error < self.rtol and self.absolute_error < self.atol:
            converged = True

        return converged

    @staticmethod
    def _sweeps(*sweep_args, **sweep_kw):
        sweep_objects = []
        sweep_vars = []
        sweep_vectors = []
        for pair in sweep_args:
            obj, param_dict = pair
            for key, val in param_dict.items():
                sweep_objects.append(obj)
                sweep_vars.append(key)
                sweep_vectors.append(val)
        for key, val in sweep_kw.items():
            sweep_objects.append(None)
            sweep_vars.append(key)
            sweep_vectors.append(val)

        sweep_vals = list(product(*sweep_vectors))

        return sweep_objects, sweep_vars, sweep_vals

    def _assign_sweep_vars(self, sweep):
        for obj, var, val in zip(self.sweep_objects, self.sweep_vars, sweep):
            if obj is None:
                obj = self.source

            if not hasattr(obj, var):
                raise ValueError(f"""Object {obj} does not have attribute {var}.
                                 Invalid sweep variable. Available default variables are 
                                 "phi", "theta", "wavelength", "pTEM"'
                                 """)
            setattr(obj, var, val)

    def _couple_source(self):
        self.source.layer = self.layer_stack.incident_layer
        self.layer_stack.source = self.source

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

        if self.TMMSimulation:
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

        new_results = Results(new_results)
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
        tempResults['crystal'] = deepcopy(self.base_crystal)
        tempResults['source'] = deepcopy(self.source)
        tempResults['S'] = deepcopy(self.SGlobal)
        tempResults['Si'] = deepcopy(self.Si)

        if self.TMMSimulation:
            tempResults['rTE'] = self.rTEM[0]
            tempResults['rTM'] = self.rTEM[1]
            rho = tempResults['rTM'] / tempResults['rTE']
            tempResults['tanPsi'] = np.abs(rho)
            tempResults['cosDelta'] = cos(np.angle(rho))
            tempResults['delta'] = np.angle(rho)

        self.results.append(tempResults)

    @property
    def _k_dimension(self):
        if self.TMMSimulation:
            k_dim = 1
        else:
            k_dim = np.prod(self.n_harmonics)

        return k_dim

    @property
    def _s_element_dimension(self):
        s_dim = self._k_dimension * 2
        return s_dim
        
    @property
    def _s_element_shape(self):
        s_dim = self._s_element_dimension
        s_shape = (s_dim, s_dim)
        return s_shape



    def _k_matrices(self):
        """
        Sets up the Kx, Ky, and Kz matrices for solving the simulation once the source, crystal, and
        number harmonics are known.
        """
        self.Kx = kx_matrix(self.source, self.base_crystal, self.n_harmonics)
        self.Ky = ky_matrix(self.source, self.base_crystal, self.n_harmonics)
        self.layer_stack.Kx = self.Kx
        self.layer_stack.Ky = self.Ky

        self.KzReflectionRegion = self.layer_stack.incident_layer.Kz_backward()
        self.KzTransmissionRegion = self.layer_stack.transmission_layer.Kz_forward()

    def _outer_matrices(self):
        self.WReflectionRegion = complexIdentity(self._s_element_dimension)
        self.WTransmissionRegion = complexIdentity(self._s_element_dimension)

    def _gap_matrices(self):
        self.layer_stack.set_gap_layer()
        self.KzGapRegion = self.layer_stack.gapLayer.Kz_gap()

    def _inner_s_matrix(self):
        for i, layer in enumerate(self.layer_stack.internal_layers):
            self.Si[i] = layer.S_matrix()
            self.SGlobal = redheffer_product(self.SGlobal, self.Si[i])

    def _global_s_matrix(self):
        self.STransmission = self.layer_stack.transmission_layer.S_matrix()
        self.SReflection = self.layer_stack.incident_layer.S_matrix()
        self.SGlobal = redheffer_product(self.SGlobal, self.STransmission)
        self.SGlobal = redheffer_product(self.SReflection, self.SGlobal)

    def _initialize(self):
        if self.base_crystal is None:
            self.TMMSimulation = True
        else:
            self.TMMSimulation = False

        self.SGlobal = S_matrix_transparent(self._s_element_shape)
        self.rx, self.ry, self.rz = None, None, None
        self.tx, self.ty, self.tz = None, None, None
        self.R, self.T, self.RTot, self.TTot, self.CTot = None, None, None, None, None
        self.Si = [None for _ in range(len(self.layer_stack.internal_layers))]

        self._couple_source()
        self.layer_stack.set_convolution_matrices(self.n_harmonics)
        self._k_matrices()
        self._gap_matrices()
        self._outer_matrices()
