import numpy as np
import math
from softx_functions import energy_conversion


class GratingSpectrometerCalibration():
    """
    Class that contains all the parameters necessary for energy_conversion function.
    """
    def __init__(self):
        self.x0 = 444.0
        self.s = 0.05

        self.r1 = 563.2
        self.alpha0 = 3 * math.pi / 180
        self.y0 = self.r1 * math.tan(self.alpha0)

        self.pos0 = 0.0
        self.currentPos = 85.0

        self.energyMin = 1
        self.energyMax = 350
        self.energyBinWidth = 0.1

        self.d0 = 0.000083

        self.deltaP = np.array([self.currentPos - self.pos0])
        self.image = np.empty(2, dtype=float)
        self.imageIntegrated = np.empty(1, dtype=float)

    def setx0(self, x0):
        self.x0 = x0

    def setS(self, s):
        self.s = s

    def setAlpha0(self, alpha0):
        self.alpha0 = alpha0

    def setPos0(self, pos0):
        self.pos0 = pos0

    def setCurrentPos(self, currentPos):
        self.currentPos = currentPos

    def setDeltaP (self, deltaP):
        self.deltaP = deltaP

    def setEnergyMin(self, energyMin):
        self.energyMin = energyMin

    def setEnergyMax(self, energyMin):
        self.energyMax = energyMin

    def setEnergyBinWidth(self, energyBinWidth):
        self.energyBinWidth = energyBinWidth

    def setd0(self, d0):
        self.d0 = d0

    def setImage(self, image):
        self.image = image

    def setImageIntegrated(self, imageIntegrated):
        self.imageIntegrated = imageIntegrated

    def setGSCalibration(self):
        self.deltaP = np.array([self.currentPos - self.pos0])
        # error if both currentPos and pos0 are 0.0
        self.y0 = self.r1 * math.tan(self.alpha0)
        self.energy, self.spectralIntensityVal, self.spectralImageVal = energy_conversion(
            self.d0, self.alpha0, self.s, self.x0, self.y0, self.deltaP, self.energyMin, self.energyMax,
            self.energyBinWidth, self.image, self.imageIntegrated)
        return self.energy, self.spectralIntensityVal, self.spectralImageVal

    def setGSCalibrationMultiple(self):
        self.y0 = self.r1 * math.tan(self.alpha0)
        self.energy, self.spectralIntensityVal, self.spectralImageVal = energy_conversion(
            self.d0, self.alpha0, self.s, self.x0, self.y0, self.deltaP, self.energyMin, self.energyMax,
            self.energyBinWidth, self.image, self.imageIntegrated)
        return self.energy, self.spectralIntensityVal, self.spectralImageVal
