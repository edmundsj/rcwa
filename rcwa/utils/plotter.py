import numpy as np
import matplotlib.pyplot as plt
import warnings


class Plotter:
    @staticmethod
    def plotRTSpectra(result):
        """
        Plots the reflectance, transmittance, and conservation spectra from a set of simulatiotn results

        :param resultsList: List of results from simulation output
        :returns fig, ax: Figure and axes objects for plot
        """
        warnings.warn('''Plotter is now deprecated and will be removed in a future version. 
                To plot results use results.plot()''')

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
