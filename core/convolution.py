# Computes the convolution matrices for multi-dimensional simulations.
import numpy as np
from matrices import *
import math as math

from shorthand import *

# This function is too long.
#def generateConvolutionMatrix(cellData, numberHarmonics):
#    # Now, for the code below to work for any dimension we need to add dimensions to A
#    dataDimension = len(cellData.shape);
#    if(dataDimension == 2):
#        numberHarmonics = (numberHarmonics + (1,))
#    (P, Q, R) = numberHarmonics
#
#    convolutionMatrixSize = P*Q*R;
#    convolutionMatrixShape = (convolutionMatrixSize, convolutionMatrixSize);
#    convolutionMatrix = complexZeros(convolutionMatrixShape)
#
#    cellData = reshapeLowDimensionalData(cellData);
#    (Nx, Ny, Nz) = cellData.shape;
#    zeroHarmonicsLocation = [math.floor(Nx/2), math.floor(Ny/2), math.floor(Nz/2)]
#
#    cellFourierRepresentation = fftn(cellData);
#
#    for rrow in range(R):
#        for qrow in range(Q):
#            for prow in range(P):
#                row = rrow*Q*P + qrow*P + prow;
#                for rcol in range(R):
#                    for qcol in range(Q):
#                        for pcol in range(P):
#                            col = rcol*Q*P + qcol*P + pcol;
#                            # Get the desired harmonics relative to the 0th-order harmonic.
#                            desiredHarmonics = [prow - pcol, qrow - qcol, rrow - rcol]
#
#                            # Get those harmonic locations from the zero harmonic location.
#                            desiredHarmonicsLocation = zeroHarmonicsLocation + desiredHarmonics
#
#                            convolutionMatrix[row][col] = \
#                                cellFourierRepresentation[desiredHarmonicsLocation[0]][desiredHarmonicsLocation[1]][desiredHarmonicsLocation[2]];
#    if convolutionMatrix.shape == (1, 1):
#        convolutionMatrix = convolutionMatrix[0][0]
#
#    return convolutionMatrix;
