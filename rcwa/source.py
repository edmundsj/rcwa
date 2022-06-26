from rcwa.shorthand import *
from rcwa import freeSpaceLayer


class Source:
    """
    Class for defining monochromatic excitation source

    :param wavelength: The wavelength (in microns, or your preferred length unit due to the scale-invariance of Maxwell's equations.
    :param theta: Angle with respect to the vector normal to the layer stack, in radians
    :param phi: Rotation angle amount the vector normal to the layer stack
    :param pTEM: Polarization vector for TE/TM polarization fraction (can be complex)
    :param layer: Layer source is located in
    """
    def __init__(self, wavelength=2*pi, theta=0, phi=0, pTEM=[1, 1], layer=freeSpaceLayer):
        self.freeSpaceWavelength = wavelength
        self.layer = layer
        self.k0 = 2*pi / self.freeSpaceWavelength
        self._phi = phi
        self._theta = theta
        self._set_k_incident()

        self._pTEM = pTEM / norm(pTEM)

        self._set_tem_vectors()
        self.pX = (self.pTE*self.aTE + self.pTM*self.aTM)[0]
        self.pY = (self.pTE*self.aTE + self.pTM*self.aTM)[1]

    def __eq__(self, other):
        if not isinstance(other, Source):
            return NotImplemented
        return self.freeSpaceWavelength == other.freeSpaceWavelength and \
               self.k0 == other.k0 and \
               self.theta == other.theta and \
               self.phi == other.phi and \
               self.pTE == other.pTE and \
               self.pTM == other.pTM and \
               self.pX == other.pX and \
               self.pY == other.pY and \
               self.k_incident == other.k_incident and \
               self.aTE == other.aTE and \
               self.aTM == other.aTM

    def __str__(self):
        return f'wavelength: {self.freeSpaceWavelength:.3f}, (theta, phi) = ({self.theta:.4f}, {self.phi:.4f})\n' + \
                f'kIncident: {self.k_incident}, polarization: ({self.pTE:.3f}, {self.pTM:.3f})\n' + \
                f'TEM vector: ({self.aTE}, {self.aTM})\n'

    def __repr__(self):
        return f'wavelength: {self.freeSpaceWavelength:.3f}, (theta, phi) = ({self.theta:.4f}, {self.phi:.4f})' + \
                f'kIncident: {self.k_incident}, polarization: ({self.pTE:.3f}, {self.pTM:.3f})' + \
                f'TEM vector: ({self.aTE}, {self.aTM})'

    def getRepresentationVector(self):
        return complexArray([self.wavelength, self.k0, self.theta, self.phi, self.pTE, self.pTM,
                             self.pX, self.pY, self.k_incident[0], self.k_incident[1], self.k_incident[2],
                             self.aTE[0], self.aTE[1], self.aTM[0], self.aTM[1]])

    @property
    def wavelength(self):
        if isinstance(self.freeSpaceWavelength, np.ndarray):
            return self.freeSpaceWavelength[0]
        else:
            return self.freeSpaceWavelength

    @wavelength.setter
    def wavelength(self, wavelength):
        self.freeSpaceWavelength = wavelength
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
        self._set_k_incident()
        self._set_tem_vectors()

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
        n = sqrt(self.layer.er*self.layer.ur);
        kx = n * sin(self.theta) * cos(self.phi);
        ky = n * sin(self.theta) * sin(self.phi);
        kz = n * cos(self.theta);
        self._k_incident = complexArray([kx, ky, kz])

zeroSource = Source(float("inf"))