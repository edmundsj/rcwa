"""
I think the way this currently works is too convoluted. It needs to be refactored to be understandable.

"""
import csv
import numpy as np
import pandas as pd
import yaml
import rcwa
import os
from rcwa.shorthand import *

class Material:
    """
    Material class for defining materials permittivity / permeability / refractive index as a function of wavelength / angle.

    :param filename: File containing n/k data for the material in question
    """

    def __init__(self, material_name=None, er=1, ur=1, n=None, filename=None, source=None):
        self.name = ''
        self.source=source
        self.dispersive = False
        self.dispersion_formula = None

        if material_name != None or filename !=None:
            self.dispersive = True
            self.loadFromDatabase(material_name, filename=filename)

        elif material_name == None:
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
    def loadFromDatabase(self, material_name, filename=None):

        self.extractMaterialDatabase()

        if filename != None:
            file_to_load = os.path.join(rcwa.nkLocation, 'data', filename)

        if material_name in self.materials.keys():
            file_to_load = os.path.join(rcwa.nkLocation, 'data', self.materials[material_name])

        self.loadMaterialFile(file_to_load)


    def extractMaterialDatabase(self):
        database_filename = rcwa.nkLocation + '/library.yml'
        self.materials = {}
        with open(database_filename) as database_file:
            database_list = yaml.load(database_file, Loader=yaml.FullLoader)

        main_content = database_list[0]['content']
        for i in range(len(main_content)):
            book_or_divider = database_list[0]['content'][i]

            if 'BOOK' in book_or_divider.keys():
                material = main_content[i]['BOOK']
                material_content = main_content[i]['content']

                for j in range(len(material_content)):
                    if 'PAGE' in material_content[j].keys():
                        file_location = material_content[j]['data']
                        self.materials[material] = file_location
                        break

    def loadMaterialFile(self, filename):
        with open(filename) as fn:
            material_file = yaml.load(fn, Loader=yaml.FullLoader)
            material_data = material_file['DATA'][0]
            if material_data['type'] == 'tabulated nk':
                self.dispersion_type = 'tabulated'
                self.load_nk_table_data(material_data)
            else:
                self.dispersion_type = 'formula 2'
                self.load_nk_formula_data(material_data)

    def load_nk_formula_data(self, data_dict):
        if data_dict['type'] == 'formula 1':
            self.load_nk_formula_1_data(data_dict)
        elif data_dict['type'] == 'formula 2':
            self.load_nk_formula_2_data(data_dict)

    def load_nk_formula_1_data(self, data_dict):
        coeffs = data_dict['coefficients'].split()
        coeffs = [float(x) for x in coeffs]
        A = coeffs[0]
        num_terms = int((len(coeffs) - 1) / 2)
        B_coeffs = [coeffs[2*i+1] for i in range(num_terms)]
        C_coeffs = [coeffs[2*i+2] for i in range(num_terms)]
        def dispersion_formula_er(wavelength):
            L = wavelength
            b_terms = [b * L**2 / (L**2 - c**2) for b, c in zip(B_coeffs, C_coeffs)]
            full_term = 1 + A + np.sum(b_terms)
            return full_term
        def dispersion_formula_n(wavelength):
            return np.sqrt(dispersion_formula_er(wavelength))
        def dispersion_formula_ur(wavelength):
            return 1

        self._er = dispersion_formula_er
        self._ur = dispersion_formula_ur
        self._n = dispersion_formula_n

    def load_nk_formula_2_data(self, data_dict):
        coeffs = data_dict['coefficients'].split()
        coeffs = [float(x) for x in coeffs]
        A, B1, C1, B2, C2 = coeffs
        def dispersion_formula_er(wavelength):
            b1_term = B1 * wavelength **2 / (wavelength**2 - C1)
            b2_term = B2 * wavelength**2 / (wavelength**2 - C2)
            full_term = 1 + A + b1_term + b2_term
            return full_term
        def dispersion_formula_n(wavelength):
            return np.sqrt(dispersion_formula_er(wavelength))
        def dispersion_formula_ur(wavelength):
            return 1

        self._n = dispersion_formula_n
        self._er = dispersion_formula_er
        self._ur = dispersion_formula_ur

    def load_nk_table_data(self, data_dict):
        material_data = data_dict['data']
        nk_data_string = list(filter(None, material_data.split('\n')))
        split_data = [elem.split() for elem in nk_data_string]
        numerical_data = np.array(split_data, dtype=np.float64)
        data_shape = numerical_data.shape
        if data_shape[1] == 3:
            self._n = numerical_data[:,1] + 1j*numerical_data[:,2]
            self.wavelengths = numerical_data[:,0]
        elif data_shape[1] == 2:
            self.wavelengths = numerical_data[:,0]
            self._n = numerical_data[:,1]
        else:
            raise ValueError

        self._er = sq(self._n)
        self._ur = np.ones(self._n.shape)


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

    """
    Outputs the complex permittivity-permeability pair

    :param wavelength: The frequency of interest (in microns)
    :param plasma_frequency: The plasma frequency of the metal of interest
    :param damping_coefficients: A k+1 array of damping coefficients ( in electronvolts), k is number of oscillators.
    :param oscillator_amplitudes: A k+1 array of oscillator amplitudes (unitless). Zeroth element taken to be amplitude of intraband transitions.
    :param offset_coefficients: A k array of offset frequencies (in electronvolts)
    :returns: A two-tuple containing permittivity and permeability (er, ur)
    """
    def drudeLorentzModel(wavelength, plasma_frequency, damping_coefficients,
            oscillator_amplitudes, offset_coefficients):

        frequency = 1.2398419843320028 / wavelength # converts from microns to eV
        er_free = 1 - oscillator_amplitudes[0] * plasma_frequency * plasma_frequency / \
                (frequency * (frequency - 1j*damping_coefficients[0]))
        er_bound = 0
        for i in range(len(damping_coefficients - 1)):
            er_bound += oscillator_amplitudes[i] * plasma_frequency * plasma_frequency / \
                    (offset_coefficients*offset_coefficients - \
                    frequency*frequency + 1j*frequency*damping_coefficients[i])
        er = er_free + er_bound
        ur = 1
        return (er, ur)

    def lookupParameter(self, parameter):
        if self.dispersion_type == 'tabulated':
            return self.lookupNumeric(parameter)
        elif self.dispersion_type == 'formula 2':
            wavelength = self.source.wavelength
            return parameter(wavelength)

    def lookupNumeric(self, parameter):
        """
        Looks up a numeric value of a parameter

        :param parameter: Either _n, _er, or _ur
        """
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
            else: # We need to interpolate the wavelength. The indexOfWavelength is pointing to the *next* value
                slope = (parameter[indexOfWavelength] - parameter[indexOfWavelength-1]) / (self.wavelengths[indexOfWavelength] - self.wavelengths[indexOfWavelength-1]) # wavelength spacing between two points
                deltaWavelength = wavelength - self.wavelengths[indexOfWavelength]
                return_value = parameter[indexOfWavelength] + slope * deltaWavelength

        return return_value
