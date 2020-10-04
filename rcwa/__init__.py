__version__ = '0.0.1'
__author__ = 'Jordan Edmunds'
__email__ = 'edmundsj@uci.edu'

import os
import sys

initLocation= os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, initLocation)

import context
import source

from layer import LayerStack, Layer
from source import Source
from solver import Solver
from plotter import Plotter
