from shorthand import *
from layer import *

# WHEN WE GET BACK - THIS IS WHERE WE SHOULD PUT OUR CONVOLUTION MATRIX.
class Source:
    def __init__(self, wavelength=2*pi,theta=0, phi=0, pTEM=[1,0], layer=freeSpaceLayer):
        self.wavelength = wavelength
        self.k0 = 2*pi / wavelength
        self.theta = theta
        self.phi = phi
        self.pTE = pTEM[0]
        self.pTM = pTEM[1]
        self.setKIncident(layer)
        self.setTEMVectors()
        self.pX = (self.pTE*self.aTE + self.pTM*self.aTM)[0]
        self.pY = (self.pTE*self.aTE + self.pTM*self.aTM)[1]

    def setTEMVectors(self):
        deviceNormalUnitVector = complexArray([0,0,-1]);
        epsilon = 1e-3;

        if(abs(self.kIncident[0]) < epsilon and abs(self.kIncident[1]) < epsilon):
            self.aTE = np.array([0,1,0]);
        else:
            self.aTE = - np.cross(deviceNormalUnitVector, self.kIncident);
            self.aTE = self.aTE / norm(self.aTE);

        self.aTM = (np.cross(self.aTE, self.kIncident));
        self.aTM = (self.aTM / norm(self.aTM));

    def setKIncident(self, layer):
        n = sqrt(layer.er*layer.ur);
        kx = n * sin(self.theta) * cos(self.phi);
        ky = n * sin(self.theta) * sin(self.phi);
        kz = n * cos(self.theta);
        self.kIncident = complexArray([kx, ky, kz])

zeroSource = Source(float("inf"))
