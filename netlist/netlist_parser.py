# This is a parser for our simple netlists.
# TODO:
# 1. Incorporate phase information into Gaussian wave based on the beam waist.
# 2. Switch wavelength to a global property, and only allow people to define a single
#   wavelength per simulation. They can still do multiple wavelengths eventually, but 
#   for now this should be a global setting.

import sys
sys.path.append('core');
sys.path.append('test')

from crystal import *
from layer import *
from source import *

# EXPECTED NETLIST FORMAT
WAVELENGTH_POSITION = 1;
THETA_POSITION = 2;
PHI_POSITION = 3;
TE_POSITION = 4;
TM_POSITION = 5;
PERMITTIVITY_POSITION = 1;
PERMEABILITY_POSITION = 2;
THICKNESS_POSITION = 3;

DBGLVL = 2;

import re
import numpy as np

class NetlistParser:
    filename = '';
    largest_aperture = 0.0; # The largest aperture size, in mm

    def __init__(self, filename):
        self.filename = filename;
        self.useRefractiveIndex = False
        self.sources = []
        self.layers = []
        self.layerStack = LayerStack()
        self.processedLines = []

        f = open(self.filename, 'r');
        self.fileLines = f.readlines();
        f.close()

    def processLines(self):
        for line in self.fileLines:
            sep = '#'; # The separator for our netlist files. Also comment symbol.
            comments_removed = line.split(sep, 1)[0];   # Removes all comments
            stripped = comments_removed.strip();        # Removes trailing and preceding whitespace
            commas_removed = stripped.replace(',', ' ');  # Remove all commas
            whitespace_removed = re.sub('\s+', ' ', commas_removed);
            if(stripped != ''):
                self.processedLines.append(whitespace_removed);

    def extractDirectives(self):
        linesToDelete = []
        for i in range(len(self.processedLines)):
            if self.processedLines[i][0] == '.':
                directiveLine = self.processedLines[i]

                if directiveLine.upper() == '.REFRACTIVEINDEX':
                    self.useRefractiveIndex = True

                linesToDelete.append(i)

        linesToDelete.reverse()
        for index in linesToDelete:
            del self.processedLines[index]

    def extractLayers(self):
        linesToDelete = []
        for i in range(len(self.processedLines)):
            if self.processedLines[i][0] == 'L':
                layerLine = self.processedLines[i]
                er = 0
                ur = 0
                thickness = 0
                lineChunks = self.processedLines[i].split(' ');
                numberChunks = len(lineChunks)
                if self.useRefractiveIndex == True:
                    er = sq(self.stripUnits(lineChunks[PERMITTIVITY_POSITION]))
                    ur = 1
                    if numberChunks >= 3:
                        thickness = self.stripUnits(lineChunks[PERMEABILITY_POSITION])
                else:
                    er = self.stripUnits(lineChunks[PERMITTIVITY_POSITION]);
                    ur = self.stripUnits(lineChunks[PERMEABILITY_POSITION]);
                    if numberChunks >= 2:
                        thickness = self.stripUnits(lineChunks[THICKNESS_POSITION])

                self.layers.append(Layer(er, ur, thickness))
                linesToDelete.append(i)

        linesToDelete.reverse()
        for index in linesToDelete:
            del self.processedLines[index]

        self.layerStack = LayerStack(self.layers)

    # FINALLY - EXTRACT THE SOURCES AND PUT IT ALL TOGETHER. ADD WAVELENGTH SWEEPS.
    def extractSources(self):
        linesToDelete = []
        for i in range(len(self.processedLines)):
            if self.processedLines[i][0] == 'W':
                lineChunks = self.processedLines[i].split(' ');
                numberArguments = len(lineChunks) - 1
                self.startWavelength = self.stripUnits(lineChunks[WAVELENGTH_POSITION]);
                self.stepWavelength = 1.0
                self.stopWavelength = self.startWavelength + self.stepWavelength
                theta = 0*deg
                phi = 0*deg
                pTE = 1
                pTM = 0

                if numberArguments == 3:
                    self.stopWavelength = self.stripUnits(lineChunks[WAVELENGTH_POSITION + 1])
                    self.stepWavelength = self.stripUnits(lineChunks[WAVELENGTH_POSITION + 2])
                elif numberArguments == 4:
                    theta = self.stripUnits(lineChunks[WAVELENGTH_POSITION + 1])
                    pTE = self.stripUnits(lineChunks[WAVELENGTH_POSITION + 2])
                    pTM = self.stripUnits(lineChunks[WAVELENGTH_POSITION + 3])
                elif numberArguments == 5:
                    theta = self.stripUnits(lineChunks[WAVELENGTH_POSITION +1])
                    phi = self.stripUnits(lineChunks[WAVELENGTH_POSITION +2])
                    pTE = self.stripUnits(lineChunks[WAVELENGTH_POSITION +3])
                    pTM = self.stripUnits(lineChunks[WAVELENGTH_POSITION +4])
                elif numberArguments == 6:
                    self.stopWavelength = self.stripUnits(lineChunks[WAVELENGTH_POSITION + 1])
                    self.stepWavelength = self.stripUnits(lineChunks[WAVELENGTH_POSITION + 2])
                    theta = self.stripUnits(lineChunks[WAVELENGTH_POSITION +3])
                    pTE = self.stripUnits(lineChunks[WAVELENGTH_POSITION +4])
                    pTM = self.stripUnits(lineChunks[WAVELENGTH_POSITION +5])
                elif numberArguments == 7:
                    self.stopWavelength = self.stripUnits(lineChunks[WAVELENGTH_POSITION + 1])
                    self.stepWavelength = self.stripUnits(lineChunks[WAVELENGTH_POSITION + 2])
                    theta = self.stripUnits(lineChunks[WAVELENGTH_POSITION +3])
                    phi = self.stripUnits(lineChunks[WAVELENGTH_POSITION +4])
                    pTE = self.stripUnits(lineChunks[WAVELENGTH_POSITION +5])
                    pTM = self.stripUnits(lineChunks[WAVELENGTH_POSITION +6])

                # Normalize pTE and pTM so the sum of their squares is one.
                pTEM = np.array([pTE, pTM]);
                pTEM = pTEM / np.linalg.norm(pTEM);

                source = Source(self.startWavelength, theta, phi, pTEM, self.layerStack.reflectionLayer)
                self.sources.append(source); # Appends the wavelength, in microns.
                linesToDelete.append(i)

        linesToDelete.reverse()
        for index in linesToDelete:
            del self.processedLines[index]

    def parseNetlist(self):
        self.processLines()
        self.extractDirectives()
        self.extractLayers()
        self.extractSources()

    # Stripts the units after the number. Will want to enhance this so we can
    # add units like 'degrees' and other arbitrary things.
    def stripUnits(self, text):
        # First, compress everything so there are no spaces left
        units_begin = None; # By default, assume there are no units.
        text = text.replace(" ", "");
        i = 0;

        # Now, find where the digits end and the units begin.
        contains_complex = False;
        for ch in text:
            if(ch == 'j'):
                contains_complex = True;
            if(((ch >= 'a' and ch <= 'z') or (ch >= 'A' and ch <= 'Z')) and ch != 'j'):
                units_begin = i;
                break;
            i += 1;

        # Separate the text into text containing units and text containing digits.
        units_text = text[units_begin:]
        digits_text = text[0:units_begin];

        if(units_begin == 0):
            units_text = '';

        multiplier = 1;
        if(units_text == 'mm'):
            multiplier = 1e3;
        elif(units_text == 'um'):
            multiplier = 1;
        elif(units_text == 'nm'):
            multiplier = 1e-3;
        elif(units_text == 'rad'):
            multiplier = 1;
        elif(units_text == 'deg'):
            multiplier = np.pi / 180.0;
        else:
            multiplier = 1;

        if(contains_complex == True):
            return multiplier * complex(digits_text);
        else:
            return multiplier * float(digits_text);
