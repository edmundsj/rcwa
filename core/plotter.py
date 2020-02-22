from results import *
from crystal import *
import matplotlib.pyplot as plt

class Plotter:
    @staticmethod
    def plotRTSpectra(resultsList):
        xData = np.array([result.source.wavelength for result in resultsList])
        yDataReflection = np.array([result.RTot for result in resultsList])
        yDataTransmission = np.array([result.TTot for result in resultsList])
        plt.plot(xData, yDataReflection)
        plt.plot(xData, yDataTransmission)
        plt.plot(xData, yDataReflection + yDataTransmission)
        plt.show()

    def plotReflectionSpectra(resultsList):
        xData = np.array([result.source.wavelength for result in resultsList])
        yDataReflection = np.array([result.RTot for result in resultsList])
        plt.plot(xData, yDataReflection)
