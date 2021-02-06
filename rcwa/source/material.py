import context
import csv
import numpy as np
import pandas as pd
from shorthand import *

class Material:
	"""
	Material class for defining materials permittivity / permeability / refractive index as a function of wavelength / angle.

	:param filename: File containing n/k data for the material in question
	"""

	def __init__(self, material_name='', er=1, ur=1, n=None, filename='', source=None):
		self.name = ''
		self.source=source
		filename = context.nkLocation + '/' + material_name + '.csv'

		if material_name != '':
			self.parseCSV(filename)

		elif material_name == '':
			self.dispersive = False
			if n == None: # If the refractive index is not defined, go with the permittivity
				self._er = er
				self._ur = ur
				self._n = np.sqrt(er*ur)
			else: # If the refractive index is defined, ignore the permittivity and permeability
				self._n = n
				self._er = np.square(n)
				self._ur = 1

	"""
	Parses data from a CSV file into a set of numpy arrays.

	:param filename: File containing n/k data for material in question
	"""
	def parseCSV(self, filename):

		data = pd.read_csv(filename, skiprows=2)
		self.wavelengths = data['Wavelength (um)'].to_numpy()

		self._n = (data['n'] + 1j*data['k']).to_numpy()
		data_size = len(self._n)
		self._er = sq(self._n)
		self._ur = np.ones(data_size)
		self.dispersive = True

	@property
	def n(self):
		if self.dispersive == False:
			return self._n
		else:
			return self.lookupParameter(self._n)

	@n.setter
	def n(self, n):
		self._n = n
		self._er = np.square(n)
		self._ur = 1
		self.dispersive = False

	@property
	def er(self):
		if self.dispersive == False:
			return self._er
		else:
			return self.lookupParameter(self._er)

	@er.setter
	def er(self, er):
		self._er = er
		self._n = np.sqrt(self._er * self._ur)
		self.dispersive = False

	@property
	def ur(self):
		if self.dispersive == False:
			return self._ur
		else:
			return self.lookupParameter(self._ur)

	@ur.setter
	def ur(self, ur):
		self._ur = ur
		self._n = np.sqrt(self._ur*self._er)
		self.dispersive = False

	def lookupParameter(self, parameter):
		wavelength = self.source.wavelength
		indexOfWavelength = np.searchsorted(self.wavelengths, wavelength)
		return_value = 0

		if wavelength > self.wavelengths[-1]: # Extrapolate if necessary
			slope = (parameter[-1] - parameter[-2]) / (self.wavelengths[-1] - self.wavelengths[-2])
			deltaWavelength = wavelength - self.wavelengths[-1]
			return_value = parameter[-1] + slope * deltaWavelength

		elif wavelength < self.wavelengths[0]: # Extrapolate the other direction if necessary
			slope = (parameter[1] - parameter[0]) / (self.wavelengths[1] - self.wavelengths[0])
			deltaWavelength = self.wavelengths[0] - wavelength
			return_value = parameter[0] - slope * deltaWavelength

		else: # Our wavelength is in the range over which we have data
			if wavelength == self.wavelengths[indexOfWavelength]: # We found the EXACT wavelength
				return_value = parameter[indexOfWavelength]
			else: # We need to interpolate the wavelength
				slope = (parameter[indexOfWavelength] - parameter[indexOfWavelength-1]) / (self.wavelengths[1] - self.wavelengths[0]) # wavelength spacing between two points
				deltaWavelength = wavelength - self.wavelengths[indexOfWavelength]
				return_value = parameter[indexOfWavelength] + slope * deltaWavelength

		return return_value
