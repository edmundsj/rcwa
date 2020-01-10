# This is a parser for our simple netlists.
# TODO:
# 1. Incorporate phase information into Gaussian wave based on the beam waist.
# 2. Switch wavelength to a global property, and only allow people to define a single
#   wavelength per simulation. They can still do multiple wavelengths eventually, but 
#   for now this should be a global setting.


# EXPECTED NETLIST FORMAT
WAVELENGTH_POSITION = 1;
THETA_POSITION = 2;
PHI_POSITION = 3;
TE_POSITION = 4;
TM_POSITION = 5;
PERMITTIVITY_POSITION = 1;
PERMEABILITY_POSITION = 2;
THICKNESS_POSITION = 3;


import re
import numpy as np

class NetlistParser:
    filename = '';
    largest_aperture = 0.0; # The largest aperture size, in mm

    def __init__(self, filename):
        self.filename = filename;

    # Stripts the units after the number. Will want to enhance this so we can
    # add units like 'degrees' and other arbitrary things.
    def stripUnits(self, text):
        # First, find where the digits end and the units begin.
        units_begin = 0;
        i = 0;
        for ch in text:
            if(ch.isdigit() == False and ch != '.'):
                units_begin = i;
                break;
            i += 1;

        # Separate the text into text containing units and text containing digits.
        units_text = text[units_begin:]
        digits_text = text[0:units_begin];

        if(units_begin == 0):
            units_text = '';

        if(units_text == 'mm'):
            return 1e3*float(digits_text);
        elif(units_text == 'um'):
            return float(digits_text);
        elif(units_text == 'nm'):
            return 1e-3*float(digits_text);
        elif(units_text == 'rad'):
            return float(digits_text);
        elif(units_text == 'deg'):
            return np.pi / 180.0 *float(digits_text);
        else:
            return float(text);

    def testLargest(self, current_aperture):
        if(current_aperture > self.largest_aperture):
            return current_aperture;
        else:
            return self.largest_aperture;

    def parseNetlist(self):
        relative_permittivities = [];
        relative_permeabilities = [];
        thicknesses = [];
        sources = []; # excitation sources we are sending at our layer stack. For now, only plane waves are allowed.
        f = open(self.filename, 'r');
        sep = '#'; # The separator for our netlist files. Also comment symbol.
        f1 = f.readlines();
        processed_lines = [];
        sim_directives = [];


        for line in f1:
            comments_removed = line.split(sep, 1)[0];   # Removes all comments
            stripped = comments_removed.strip();        # Removes trailing and preceding whitespace
            commas_removed = stripped.replace(',', ' ');  # Remove all commas
            whitespace_removed = re.sub('\s+', ' ', commas_removed);
            if(stripped != ''):
                processed_lines.append(whitespace_removed);

        # Loop through all the lines that actually contain something interesting
        # Split up each line into 'chunks' (typically numbers and their associated units or 
        # simulation directives
        for line in processed_lines:
            line_chunks = line.split(' ');
            name = line_chunks[0];


            # If we have a plane wave source, we need to know two pieces of information:
            # the intensity of the plane wave, and its tilt. By default (and for now)
            # we will assume the tilt is zero. If it is not, this just changes the phase
            # function we need to apply to the plane wave.
            # This is expecting a grid of x and y values. It will return a numpy array
            # Even if you feed it with a single value.
            if(line[0] == 'W'): # The line contains a plane-wave source of infinite extent.

                wavelength = self.stripUnits(line_chunks[WAVELENGTH_POSITION]);
                theta = self.stripUnits(line_chunks[THETA_POSITION]);
                phi = self.stripUnits(line_chunks[PHI_POSITION]);
                pTE = 1;
                pTM = 1;

                if(len(line_chunks)>PHI_POSITION+1): # We have additional TE and TM polarization information
                    pTE = self.stripUnits(line_chunks[TE_POSITION]);
                    pTM = self.stripUnits(line_chunks[TM_POSITION]);

                # Normalize pTE and pTM so the sum of their squares is one.
                pTEM = np.array([pTE, pTM]);
                pTEM = pTEM / np.linalg.norm(pTEM);

                sources.append((wavelength, theta, phi, pTEM[0], pTEM[1])); # Appends the wavelength, in microns.


            # This needs to be modified to include the Gaussian phase.
            elif(line[0] == 'L'): # The line contains a layer of finite or infinite extent.
                if(len(relative_permittivities) == 0):
                    er = self.stripUnits(line_chunks[PERMITTIVITY_POSITION]);
                    ur = self.stripUnits(line_chunks[PERMEABILITY_POSITION]);

                    relative_permittivities.append(er);
                    relative_permeabilities.append(ur);

                # The line is too thin to contain a thickness, and so must be our final medium.
                elif(len(relative_permittivities) != 0 and len(line_chunks) == 3):
                    er = self.stripUnits(line_chunks[PERMITTIVITY_POSITION]);
                    ur = self.stripUnits(line_chunks[PERMEABILITY_POSITION]);

                    relative_permittivities.append(er);
                    relative_permeabilities.append(ur);

                # The line corresponds to an intermediate layer.
                else:
                    er = self.stripUnits(line_chunks[PERMITTIVITY_POSITION]);
                    ur = self.stripUnits(line_chunks[PERMEABILITY_POSITION]);
                    thickness = self.stripUnits(line_chunks[THICKNESS_POSITION]);

                    relative_permittivities.append(er);
                    relative_permeabilities.append(ur);
                    thicknesses.append(thickness);

            elif(line[0] == '.'): # Simulation directive
                line_chunks[0] = line_chunks[0].upper();
                if(line_chunks[0] == '.NOPHASE'):
                    sim_directives.append('NOPHASE');
                if(line_chunks[0] == '.NO2D'):
                    sim_directives.append('NO2D');
                if(line_chunks[0] == '.NO1D'):
                    sim_directives.append('NO1D');
                if(line_chunks[0] == '.SIZEREL'):
                    sim_directives.append(('SIZEREL', float(line_chunks[1])));
                if(line_chunks[0] == '.SIZEABS'):
                    sim_directives.append(('SIZEABS', self.stripUnits(line_chunks[1])));
                if(line_chunks[0] == '.NX'):
                    sim_directives.append(('NX', int(line_chunks[1])));
                if(line_chunks[0] == '.NY'):
                    sim_directives.append(('NY', int(line_chunks[1])));
            else:
                print("ERROR: Not able to parse line:")
                print(line)
                print(line[0]);
                print("Object type not supported");


        return [relative_permittivities, relative_permeabilities, thicknesses, sources];

