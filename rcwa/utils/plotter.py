import numpy as np
import matplotlib.pyplot as plt

class Plotter:

    @staticmethod
    def plotRTSpectra(result):
        """
        Plots the reflectance, transmittance, and conservation spectra from a set of simulatiotn results

        :param resultsList: List of results from simulation output
        :returns fig, ax: Figure and axes objects for plot
        """
        xData = np.array(result['wavelength'])
        yDataReflection = np.array(result['RTot'])
        yDataTransmission = np.array(result['TTot'])
        fig, ax = plt.subplots(1, 1)
        ax.plot(xData, yDataReflection)
        ax.plot(xData, yDataTransmission)
        ax.plot(xData, yDataReflection + yDataTransmission)
        ax.legend(['R', 'T', 'Conservation'])
        ax.set_title('R and T vs wavelength')
        ax.set_xlabel('wavelength (um)')
        ax.set_ylabel('R, T')
        return fig, ax

    @staticmethod
    def plotEllipsometrySpectra(resultsList):
        """
        Plots the ellipsometry spectra from simulation results

        :param resultsList: List of results from simulation output
        :returns fig, ax: Figure and axes objects for plot
        """
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
        """
        Plots the TE/TM reflectance from simulation results

        :param resultsList: List of results from simulation output
        :returns fig, ax: Figure and axes objects for plot
        """
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
