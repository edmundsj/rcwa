import sys
sys.path.append('core');
sys.path.append('test')

import unittest
from shorthandTest import *
from matrices import *
from fresnel import *
from convolution import generateConvolutionMatrix
from matrixParser import *

S11Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/S11Layer2.txt")
np.savetxt("test/mathematica/S11Layer2.csv",S11Layer2)

S12Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/S12Layer2.txt")
np.savetxt("test/mathematica/S12Layer2.csv",S12Layer2)

S21Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/S21Layer2.txt")
np.savetxt("test/mathematica/S21Layer2.csv",S21Layer2)

S22Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/S22Layer2.txt")
np.savetxt("test/mathematica/S22Layer2.csv",S22Layer2)

S11Layer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/S11Layer1.txt")
np.savetxt("test/mathematica/S11Layer1.csv",S11Layer1)

S12Layer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/S12Layer1.txt")
np.savetxt("test/mathematica/S12Layer1.csv",S12Layer1)

S21Layer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/S21Layer1.txt")
np.savetxt("test/mathematica/S21Layer1.csv",S21Layer1)

S22Layer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/S22Layer1.txt")
np.savetxt("test/mathematica/S22Layer1.csv",S22Layer1)
