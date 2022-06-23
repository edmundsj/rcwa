import pandas as pd
import numpy as np
import os
import yaml
from rcwa.utils import nk_location

def nk_to_complex(data):

    if isinstance(data, np.ndarray):
        wavelengths = data[:,0]

        if data.shape[1] == 3:
            nk_complex = data[:,1] + 1j*data[:,2]
        elif data.shape[1] == 2:
            nk_complex = data[:,1]
        else:
            raise ValueError
    elif isinstance(data, pd.DataFrame):
        wavelengths = data.iloc[:,0].values

        if '(nm)' in data.columns[0]:
            wavelengths = wavelengths / 1000

        if data.shape[1] == 3:
            nk_complex = data.iloc[:,1].values + 1j*data.iloc[:,2].values
        elif data.shape[1] == 2:
            nk_complex = data.iloc[:,1].values
        else:
            raise ValueError

    else:
        raise NotImplementedError

    return wavelengths, nk_complex

class CSVLoader:
    def __init__(self, filename):
        self.filename = filename

    def load(self):
        raw_data = pd.read_csv(self.filename)
        wavelengths, n_dispersive = nk_to_complex(raw_data)
        er_dispersive = np.sqrt(n_dispersive)
        ur_dispersive = np.ones(er_dispersive.shape)

        return {'wavelength': wavelengths, 'n': n_dispersive, 'er': er_dispersive, 'ur': ur_dispersive}

class RIDatabaseLoader:

    def __init__(self):
        """
        Loader for RefractiveIndex.info databases
        """
        self.extract_material_database()

    def extract_material_database(self):
        database_filename = os.path.join(nk_location, 'library.yml')
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

    def load(self, filename):
        with open(filename) as fn:
            material_file = yaml.load(fn, Loader=yaml.FullLoader)
            material_data = material_file['DATA'][0]
            if material_data['type'] == 'tabulated nk':
                data = self.load_nk_table_data(material_data)
            else:
                data = self.load_nk_formula_data(material_data)

        return data

    def load_nk_formula_data(self, data_dict):
        if data_dict['type'] == 'formula 1':
            return self.load_nk_formula_1_data(data_dict)
        elif data_dict['type'] == 'formula 2':
            return self.load_nk_formula_2_data(data_dict)
        else:
            raise ValueError(f'Formula type {data_dict["type"]} not supported. Please submit a bug report with this message and the specific material you are trying to use')

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

        return {'er': dispersion_formula_er, 'ur': dispersion_formula_ur, 'n': dispersion_formula_n,
                'dispersion_type': 'formula'}

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

        return {'er': dispersion_formula_er, 'ur': dispersion_formula_ur, 'n': dispersion_formula_n,
                'dispersion_type': 'formula'}


    def load_nk_table_data(self, data_dict):
        material_data = data_dict['data']
        nk_data_string = list(filter(None, material_data.split('\n')))
        split_data = [elem.split() for elem in nk_data_string]
        numerical_data = np.array(split_data, dtype=np.float64)

        wavelengths, n_dispersive = nk_to_complex(numerical_data)

        return {'er': np.square(n_dispersive), 'ur': np.ones(n_dispersive.shape), 'n': n_dispersive,
                'dispersion_type': 'tabulated', 'wavelength': wavelengths}
