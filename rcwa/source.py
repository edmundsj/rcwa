from rcwa.shorthand import *
from rcwa.utils import k_vector
from rcwa import Layer


class Source:
    """
    Class for defining monochromatic excitation source

    :param wavelength: The wavelength (in microns, or your preferred length unit due to the scale-invariance of Maxwell's equations.
    :param theta: Angle with respect to the vector normal to the layer stack, in radians
    :param phi: Rotation angle amount the vector normal to the layer stack
    :param pTEM: Polarization vector for TE/TM polarization fraction (can be complex)
    :param layer: Layer source is located in
    """
    def __init__(self, wavelength=2*pi, theta=0, phi=0, pTEM=[1, 1], layer=Layer(er=1, ur=1)):
        self._free_space_wavelength = wavelength
        self.layer = layer
        self.k0 = 2*pi / self._free_space_wavelength
        self._phi = phi
        self._theta = theta
        self._set_k_incident()

        self._pTEM = pTEM / norm(pTEM)

        self._set_tem_vectors()
        self.pX = (self.pTE*self.aTE + self.pTM*self.aTM)[0]
        self.pY = (self.pTE*self.aTE + self.pTM*self.aTM)[1]

    def __eq__(self, other):
        if not isinstance(other, Source):
            raise ValueError(f'Cannot compare Source() and non-source object {type(other)}')
        return self._free_space_wavelength == other._free_space_wavelength and \
               self.k0 == other.k0 and \
               self.theta == other.theta and \
               self.phi == other.phi and \
               np.all(self.pTE == other.pTE) and \
               np.all(self.pTM == other.pTM) and \
               self.pX == other.pX and \
               self.pY == other.pY and \
               np.all(self.k_incident == other.k_incident) and \
               np.all(self.aTE == other.aTE) and \
               np.all(self.aTM == other.aTM)

    def __str__(self):
        return f'wavelength: {self._free_space_wavelength:.3f}, (theta, phi) = ({self.theta:.4f}, {self.phi:.4f})\n' + \
                f'kIncident: {self.k_incident}, polarization: ({self.pTE:.3f}, {self.pTM:.3f})\n' + \
                f'TEM vector: ({self.aTE}, {self.aTM})\n'

    @property
    def wavelength(self):
        return self._free_space_wavelength

    @wavelength.setter
    def wavelength(self, wavelength):
        self._free_space_wavelength = wavelength
        self.k0 = 2 * pi / wavelength
        self._set_k_incident()

    def _set_tem_vectors(self):
        deviceNormalUnitVector = complexArray([0, 0, -1])
        epsilon = 1e-3

        k_norm = self.k_incident / np.linalg.norm(self.k_incident)

        if abs(k_norm[0]) < epsilon and abs(k_norm[0]) < epsilon:
            self.aTE = np.array([0, 1, 0])
        else:
            self.aTE = - np.cross(deviceNormalUnitVector, k_norm)
            self.aTE = self.aTE / norm(self.aTE)

        self.aTM = np.cross(self.aTE, k_norm)
        self.aTM /= norm(self.aTM)
        self.ATEM = np.vstack((self.aTE, self.aTM)) # matrix that goes from x/y basis to TE/TM basis

    @property
    def pTEM(self):
        return self._pTEM

    @property
    def pTE(self):
        return self.pTEM[0]

    @property
    def pTM(self):
        return self.pTEM[1]

    @pTEM.setter
    def pTEM(self, pTEM):
        self._pTEM = pTEM / np.linalg.norm(pTEM)

    @property
    def phi(self):
        return self._phi

    @phi.setter
    def phi(self, phi):
        self._phi = phi
        self._set_k_incident()

    @property
    def theta(self):
        return self._theta

    @theta.setter
    def theta(self, theta):
        self._theta = theta
        self._set_k_incident()

    @property
    def k_incident(self):
        return self._k_incident

    def _set_k_incident(self):
        self._k_incident = k_vector(self, self.layer, normalize=True)

zeroSource = Source(float("inf"))