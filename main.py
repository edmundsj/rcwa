# Author: Jordan Edmunds, Ph.D. Student, UC Berkeley
# Contact: jordan.e@berkeley.edu
# Creation Date: 11/01/2019
#
# TODO:
import sys
sys.path.append('./core')
sys.path.append('./netlist')

import numpy as np
import scipy as sp
import scipy.linalg
import sys
from core.matrices import *
from core.source import *
from core.layer import *
from core.solver import *
from netlist.netlist_parser import *
import matplotlib.pyplot as plt
from plotter import Plotter

arguments = len(sys.argv) - 1; # The number of arguments
netlistDirectory = './netlist/predictions/'
netlist1 = netlistDirectory + 'AlN_unmodulated.txt'
netlist2 = netlistDirectory + 'AlN_modulated.txt'
netlist_location = './netlist/predictions/sample_netlist.txt';
if(arguments > 0):
    print(f"Using user defined netlist {sys.argv[1]}")
    netlist_location = sys.argv[1];
print("Parsing netlist... ");

parser1 = NetlistParser(netlist1);
parser1.parseNetlist();
print("Done!")
parser2 = NetlistParser(netlist2)
parser2.parseNetlist();

print("Solving system...")
TMMSolver1 = RCWASolver(parser1.layerStack, parser1.sources[0], (1, 1))
#TMMSolver2 = RCWASolver(parser2.layerStack, parser2.sources[0], (1, 1))
wavelengths = np.arange(parser1.startWavelength, parser1.stopWavelength + parser1.stepWavelength,
        parser1.stepWavelength)

TMMSolver1.Solve(wavelengths=wavelengths)
#TMMSolver2.Solve(wavelengths=wavelengths)
Plotter.plotReflectionSpectra(TMMSolver1.results)
#Plotter.plotReflectionSpectra(TMMSolver2.results)
plt.show()
print("Done!")
