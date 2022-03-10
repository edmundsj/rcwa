import numpy as np
from rcwa import Crystal
import matplotlib.pyplot as plt

class Plotter:
    @staticmethod
    def plotRTSpectra(resultsList):
        xData = np.array([result['source'].wavelength for result in resultsList])
        yDataReflection = np.array([result['RTot'] for result in resultsList])
        yDataTransmission = np.array([result['TTot'] for result in resultsList])
        fig, ax = plt.subplots(1, 1)
        ax.plot(xData, yDataReflection)
        ax.plot(xData, yDataTransmission)
        ax.plot(xData, yDataReflection + yDataTransmission)
        ax.legend(['R', 'T', 'Conservation'])
        ax.set_title('R and T vs wavelength')
        ax.set_xlabel('wavelength (um)')
        ax.set_ylabel('R, T')

    @staticmethod
    def plotEllipsometrySpectra(resultsList):
        xData = np.array([result['source'].wavelength for result in resultsList])
        tanPsiData = np.array([result['tanPsi'] for result in resultsList])
        cosDeltaData = np.array([result['cosDelta'] for result in resultsList])
        fig, ax = plt.subplots(2, 1)
        ax[0].plot(xData, tanPsiData)
        ax[1].plot(xData, cosDeltaData)
        ax[0].set_title('Tan(\u03A8) vs wavelength')
        ax[0].set_ylabel('Tan(\u03A8)')
        ax[1].set_title('Cos(\u0394) vs wavelength')
        ax[1].set_xlabel('wavelength (um)')
        ax[1].set_ylabel('Cos(\u0394)')
        fig.suptitle('Ellipsometry Data vs. Wavelength', fontsize=16)

    @staticmethod
    def plotRTEMSpectra(resultsList):
        xData = np.array([result.source.wavelength for result in resultsList])
        RTEData = np.array([np.sq(np.abs(result.rTE)) for result in resultsList])
        RTMData = np.array([np.sq(np.abs(result.rTM)) for result in resultsList])
        fig, ax = plt.subplots(2, 1)
        ax[0].plot(xData, RTEData)
        ax[1].plot(xData, RTMData)
        ax[0].set_title('RTE and RTM vs wavelength')
        ax[0].set_ylabel('RTE')
        ax[1].set_title('RTM vs wavelength')
        ax[1].set_xlabel('wavelength (um)')
        ax[1].set_ylabel('RTM')
        fig.suptitle('RTE / TM vs wavelength', fontsize=16)
