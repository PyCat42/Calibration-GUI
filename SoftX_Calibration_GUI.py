# -*- coding: utf-8 -*-
"""
Created on Sat Aug 24 2024

Last modification Sep 4 2024

NOTE: program still under modification

@author: Daliborka Hranjec (daliborka.hranjec@gmail.com)

For all further documentation refer to the softx_calibration program.
"""

import os, sys
import matplotlib.pyplot as plt
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5 import uic
from PyQt5.QtCore import Qt, pyqtSlot, pyqtSignal
from PyQt5.Qt import QGraphicsView, QSettings, QGraphicsScene, QPixmap, QWidget, QTableWidgetItem, QPushButton, QBrush, QColor, QKeySequence, QResource
import numpy as np
from math import pi, tan
import pyqtgraph as pg
import re
from datetime import datetime
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from lmfit import minimize, Parameters
import tables
from softx_functions import energy_conversion
from grating_spectrometer_calibration import GratingSpectrometerCalibration

""" --------- SOFTX CALIBRATION PROGRAM FUNCTIONS --------- """

def read_h5_data(str_filepath):
    hdf5_file = tables.open_file(str_filepath, mode='r')
    accumulated_image = hdf5_file.get_node('/data/accumulated_image').read()
    hdf5_file.close()
    return accumulated_image


def gaus(x, a, mean, sigma, baseline):
    return a*np.exp(-(x-mean)**2/(2*sigma**2)) + baseline


def gaus_parameters(to_fit, points):
    n1 = len(to_fit)
    mean_val = sum(points * to_fit) / sum(to_fit)  # weighted average
    sigma_val = np.sqrt(sum((to_fit - np.mean(to_fit)) ** 2 / n1))
    return mean_val, sigma_val


""" --------- PARAMETERS --------- """

# physical constants
c = 299792458e3  # [mm/s]
h_Planck = 6.62606957e-34  # [J s]
eV = 1.602176565e-19  # [J]

# harmonics spacing
wl = 1400e-6  # mm central wavelength
freq = c / wl  # [Hz]
hh_space_eV = h_Planck * freq * 2 / eV  # [eV]

""" --------- CLASSES --------- """


class SoftXCalibrationGUI(QtWidgets.QMainWindow):
    """Main class that sets and runs entire calibration"""
    def __init__(self):
        super().__init__()

        uic.loadUi(os.path.join(os.path.dirname(__file__), 'SoftXCalibrationGUI.ui'), self)
        self.setWindowTitle('SoftX Calibration')
        s = QSettings()

        self._selectLowEnergyRbtn.setChecked(True)
        self.d0 = 0.00083
        self.alpha0 = 3 * pi / 180
        self._selectLowEnergyRbtn.toggled.connect(self.onSelectLowEnergyToggled)
        self._selectHighEnergyRbtn.toggled.connect(self.onSelectHighEnergyToggled)

        self.r1 = 563.2
        self.y0 = self.r1 * tan(self.alpha0)

        self._x0Spn.setMaximum(1500)
        self._useOldx0AndsRbtn.setChecked(True)
        self._first0thPathLine.setEnabled(False)
        self._browseFirst0thBtn.setEnabled(False)
        self._setFirst0thBtn.setEnabled(False)
        self._second0thPathLine.setEnabled(False)
        self._browseSecond0thBtn.setEnabled(False)
        self._setSecond0thBtn.setEnabled(False)
        self._sSpn.setDecimals(4)

        self._useOldx0AndsRbtn.toggled.connect(self.onUseOldx0AndsToggled)
        self._x0Spn.setValue(float(s.value('last x0 value', 0.0)))
        self._sSpn.setValue(float(s.value('last s value', 0.0)))
        self.zero_pos = float(s.value('last zero position', 0.0))
        self._findx0Rbtn.toggled.connect(self.handlingZero)
        self._findx0AndsRbtn.toggled.connect(self.handlingZero)
        self._browseFirst0thBtn.clicked.connect(self.onBrowseFirst0thClicked)
        self._browseSecond0thBtn.clicked.connect(self.onBrowseSecond0thClicked)
        self._setFirst0thBtn.clicked.connect(self.onSetFirst0thClicked)
        self._setSecond0thBtn.clicked.connect(self.onSetSecond0thClicked)

        self._browseForAlpha0Btn.clicked.connect(self.onBrowseForAlpha0Clicked)
        self._prominenceLine.setText("6e5")
        self._checkAlpha0Btn.clicked.connect(self.onCheckAlpha0Clicked)
        self._alpha0CheckLine.setText("0")

        self._energyMinSpn.setMaximum(1000)
        self._energyMinSpn.setValue(1)
        self._energyMaxSpn.setMinimum(self._energyMinSpn.value() + 1)
        self._energyMaxSpn.setMaximum(1000)
        self._energyMaxSpn.setValue(350)
        self._energyBinWidthSpn.setValue(0.1)
        self._energyBinWidthSpn.setSingleStep(0.05)

        self._yMinValueEdit.setText("1e6")
        self._yMaxValueEdit.setText("1e8")

        self._selectSpectraBtn.clicked.connect(self.onSelectSpectraClicked)
        self._selectBackgroundBtn.clicked.connect(self.onSelectBackgroundClicked)
        self._addImagesBtn.clicked.connect(self.onAddImagesClicked)

        self._regionSelectionBtn.clicked.connect(self.onRegionSelectionClicked)
        self._startCalibrationBtn.clicked.connect(self.onStartCalibrationClicked)
        self.str_file_numbers = ""
        self.index_list = []
        self.spectra_file_list = []
        self.image_dict = {}
        self.deltaP = np.array([])
        self.image_list = []
        self.image_integrated_list = []

        self.default_foldername = "SoftX_Calibration " + datetime.today().strftime('%Y-%m-%d %H.%M.%S')
        self._setFolderNameLine.setText(self.default_foldername)
        self._createFolderBtn.clicked.connect(self.onCreateFolderClicked)

        self._startCalibrationBtn.setEnabled(False)
        self._saveEnergyPlotChk.setEnabled(False)
        self._saveFullImageChk.setEnabled(False)
        self._saveFullSpectrumChk.setEnabled(False)
        self._save2DSpectrumChk.setEnabled(False)
        self._saveSaturationChk.setEnabled(False)
        self._saveAllChk.setEnabled(False)
        self._saveEnergyPlotChk.toggled.connect(self.onSaveEnergyPlotToggled)
        self._saveFullImageChk.toggled.connect(self.onSaveFullImageToggled)
        self._saveFullSpectrumChk.toggled.connect(self.onSaveFullSpectrumToggled)
        self._save2DSpectrumChk.toggled.connect(self.onSave2DSpectrumToggled)
        self._saveSaturationChk.toggled.connect(self.onSaveSaturationToggled)
        self._saveAllChk.toggled.connect(self.onSaveAllToggled)

    def onSelectLowEnergyToggled(self, selected):
        if selected:
            self.d0 = 0.000833
            self.alpha0 = 3 * pi / 180

    def onSelectHighEnergyToggled(self, selected):
        if selected:
            self.d0 = 0.000417
            self.alpha0 = 1 * pi / 180

    def handlingZero(self):
        self._browseFirst0thBtn.setEnabled(self._findx0AndsRbtn.isChecked() or self._findx0Rbtn.isChecked())
        self._setFirst0thBtn.setEnabled(self._findx0AndsRbtn.isChecked() or self._findx0Rbtn.isChecked())
        self._first0thPathLine.setEnabled(self._findx0AndsRbtn.isChecked() or self._findx0Rbtn.isChecked())
        self._browseSecond0thBtn.setEnabled(self._findx0AndsRbtn.isChecked())
        self._setSecond0thBtn.setEnabled(self._findx0AndsRbtn.isChecked())
        self._second0thPathLine.setEnabled(self._findx0AndsRbtn.isChecked())

    def onUseOldx0AndsToggled(self):
        s = QSettings()
        self._x0Spn.setValue(float(s.value('last x0 value', 0.0)))
        self._sSpn.setValue(float(s.value('last s value', 0.0)))

    def onBrowseFirst0thClicked(self):
        s = QSettings()
        file_name = QtWidgets.QFileDialog.getOpenFileName(self, 'open file', s.value('first zero file', r'C:/Users'), 'HDF5 files (*.h5)')
        self.zero1_file = file_name[0]
        self._first0thPathLine.setText(self.zero1_file)
        new_file_location_parts = self.zero1_file.split('/')
        new_file_location = '/'.join(new_file_location_parts[:-2])
        s.setValue('first zero file', new_file_location)

        extension = self.zero1_file.split("/")[-1]
        match = re.match('.*pos ([-|+][0-9.]+) mm', extension)
        str_pos = match.groups()

        self.zero_pos1 = float(str_pos[0])
        s = QSettings()
        s.setValue('last zero position', self.zero_pos)
        self.first_zero = read_h5_data(self.zero1_file)
        self.integrated_first_zero = self.first_zero.sum(0)

        # initial guess of the interval
        a = 100
        b = 1200
        to_fit = self.integrated_first_zero[a:b] / 1e6
        points = np.arange(a, b, 1)

        mean_val, sigma_val = gaus_parameters(to_fit, points)

        popt, pcov = curve_fit(gaus, points, to_fit, p0=[1, mean_val, sigma_val, 15])

        self.popUpZeroWindow = ZeroWindow()
        self.popUpZeroWindow.show()

        curve1 = self.popUpZeroWindow.zeroPlot.plot()
        curve1.setData(points, to_fit, symbol='+')

        curve2 = self.popUpZeroWindow.zeroPlot.plot()
        pen = pg.mkPen(color="r")
        curve2.setData(points, gaus(points, *popt), pen=pen)

    def onBrowseSecond0thClicked(self):
        s = QSettings()
        file_name = QtWidgets.QFileDialog.getOpenFileName(self, 'open file',
                                                          s.value('second zero file location', r'C:/Users'),
                                                          'HDF5 files (*.h5)')
        self.zero2_file = file_name[0]
        self._second0thPathLine.setText(self.zero2_file)
        new_file_location_parts = self.zero2_file.split('/')
        new_file_location = '/'.join(new_file_location_parts[:-2])
        s.setValue('second zero file location', new_file_location)

        extension = self.zero2_file.split("/")[-1]
        match = re.match('.*pos ([-|+][0-9.]+) mm', extension)
        str_pos = match.groups()

        self.zero_pos2 = float(str_pos[0])
        self.second_zero = read_h5_data(self.zero2_file)
        self.integrated_second_zero = self.second_zero.sum(0)

        # initial guess of the interval
        a = 100
        b = 1200
        to_fit = self.integrated_second_zero[a:b] / 1e6
        points = np.arange(a, b, 1)

        mean_val, sigma_val = gaus_parameters(to_fit, points)

        popt, pcov = curve_fit(gaus, points, to_fit, p0=[1, mean_val, sigma_val, 15])

        self.popUpZeroWindow = ZeroWindow()
        self.popUpZeroWindow.show()

        curve1 = self.popUpZeroWindow.zeroPlot.plot()
        curve1.setData(points, to_fit, symbol='+')

        curve2 = self.popUpZeroWindow.zeroPlot.plot()
        pen = pg.mkPen(color="r")
        curve2.setData(points, gaus(points, *popt), pen=pen)

    def onSetFirst0thClicked(self):
        a = int(self.popUpZeroWindow.zeroOrderCursor1.getXPos())
        b = int(self.popUpZeroWindow.zeroOrderCursor2.getXPos())

        self.popUpZeroWindow.close()

        to_fit = self.integrated_first_zero[a:b] / 1e6  # taking integrated values in the given limits
        points = np.arange(a, b, 1)  # return evenly spaced values within a given interval (start, stop, step)

        mean_val, sigma_val = gaus_parameters(to_fit, points)

        popt, pcov = curve_fit(gaus, points, to_fit, p0=[1, mean_val, sigma_val, 15])

        self.x0_1 = round(popt[1])
        self.zero_pos = self.zero_pos1 + self.x0_1 * self._sSpn.value()

    def onSetSecond0thClicked(self):
        a = int(self.popUpZeroWindow.zeroOrderCursor1.getXPos())
        b = int(self.popUpZeroWindow.zeroOrderCursor2.getXPos())

        self.popUpZeroWindow.close()

        to_fit = self.integrated_second_zero[a:b] / 1e6  # taking integrated values in the given limits
        points = np.arange(a, b, 1)  # return evenly spaced values within a given interval (start, stop, step)

        mean_val, sigma_val = gaus_parameters(to_fit, points)

        popt, pcov = curve_fit(gaus, points, to_fit, p0=[1, mean_val, sigma_val, 15])

        self.x0_2 = round(popt[1])

        self._sSpn.setValue(abs(self.zero_pos1 - self.zero_pos2) / abs(self.x0_1 - self.x0_2))

        if self.zero_pos1 == 0.0:
            self._x0Spn.setValue(self.x0_1)
        else:
            self._x0Spn.setValue(self.x0_1 + self.zero_pos1 / self._sSpn.value())
        self.zero_pos = self._x0Spn.value() * self._sSpn.value()

        s = QSettings()
        s.setValue('last x0 value', self._x0Spn.value())
        s.setValue('last s value', self._sSpn.value())

    def onBrowseForAlpha0Clicked(self):
        s = QSettings()
        file_name = QtWidgets.QFileDialog.getOpenFileName(self, 'open file', s.value('last spectrum file location', r'C:/Users'), 'HDF5 files (*.h5)')
        file_name[0].split('/')
        file_path = ('/').join(file_name[:-2])
        s.setValue('spectrum file location', file_path)
        self.alphaCheckFile = file_name[0]
        self._alpha0PathLine.setText(self.alphaCheckFile)

        name_parts = self.alphaCheckFile.split("/")
        match = re.match('[0-9.]+_pos ([-|+][0-9.]+) mm', name_parts[-1])
        str_pos = match.groups()
        self.posAlphaCheck = float(str_pos[0])

    def residual(self, params, x, peaks_pixel):
        """
        residual(params, x, peaks_pixel) - copied from softx_calibration.py
        ---------------------
        Args:
        params: Parameters which alpha0 belongs to
        x: peak number
        peaks_pixel: peaks in full spectrum
        ---------------------
        Returns: difference between expected peak value (based on minimal peak) and peak values from the full spectrum
        """
        alpha0 = params['alpha0'] * pi / 180
        peaks_mm = self.r1 * tan(alpha0) + (self.posAlphaCheck - self._x0Spn.value()) + (peaks_pixel - self._x0Spn.value()) * self._sSpn.value()
        peaks_theta = pi / 2 - np.arctan(peaks_mm / self.r1)  # angular position of peaks
        peaks_wl = self.d0 * (np.sin(pi / 2 - alpha0) - np.sin(peaks_theta))  # grating equation

        peaks_eV = h_Planck * c / peaks_wl / eV

        func = np.flip(min(peaks_eV) + hh_space_eV * x)

        return func - peaks_eV

    def onCheckAlpha0Clicked(self):
        # --- FITTING TO FIND ALPHA0 ---

        # 1 - PEAK FINDER IN PIXELS
        # PROMINENCE arg. is the most critical, it gives necessary condition for something to be considered a peak...
        # ...set a value lower than the min peak amplitude.
        alpha_image = read_h5_data(self.alphaCheckFile)
        alpha_image_integrated = alpha_image.sum(0)
        peaks_pixel = find_peaks(alpha_image_integrated, prominence=float(self._prominenceLine.text()))[0]

        # 2 - FITTING
        # NOTE: extremely sensitive to prominence, do NOT fit multiple parameters as it tries to compensate one with each other
        params = Parameters()
        if self._selectLowEnergyRbtn.isChecked():
            params.add('alpha0', value=3, min=2, max=6)  # angle in degrees, LE
        if self._selectHighEnergyRbtn.isChecked():
            params.add('alpha0', value=1, min=0.1, max=2)  # angle in degrees, HE

        x = np.arange(0, len(peaks_pixel))

        out = minimize(self.residual, params, args=(x, peaks_pixel))

        alpha0_check_deg = round(out.params[('alpha0')].value, 4)
        self._alpha0CheckLine.setText(str(alpha0_check_deg))

    def onSelectSpectraClicked(self):
        s = QSettings()
        file_names = QtWidgets.QFileDialog.getOpenFileNames(self, 'open file', s.value('last spectrum file location', r'C:/Users'), 'HDF5 files (*.h5)')
        file_name_parts = file_names[0][0].split('/')
        new_file_location = '/'.join(file_name_parts [:-1])
        s.setValue('last spectrum file location', new_file_location)

        str_file_numbers = ""
        self.current_positions_list = []
        self.current_image_list = []
        self.current_index_list = []
        self.spectra_file_list.append(file_names[0])
        for name in file_names[0]:
            self.current_image_list.append(read_h5_data(name))
            name_parts = name.split("/")
            match = re.match('([0-9.]+)_pos ([-|+][0-9.]+) mm', name_parts[-1])
            if match is None:
                continue
            match_list = match.groups()
            file_number = match_list[0]
            self.current_index_list.append(file_number)
            pos = float(match_list[1])
            self.current_positions_list.append(pos)
            if str_file_numbers == "":
                str_file_numbers = str_file_numbers + file_number
            else:
                str_file_numbers = str_file_numbers + ", " + file_number
        self._spectraPathLine.setText(str_file_numbers)

    def onSelectBackgroundClicked(self):
        s = QSettings()
        file_name = QtWidgets.QFileDialog.getOpenFileName(self, 'open file',
                                                          s.value('last background file', r'C:/Users'),
                                                          'HDF5 files (*.h5)')
        self.background_file = file_name[0]
        s.setValue('last background file', self.background_file)
        self._backgroundPathLine.setText(self.background_file)
        self.background = read_h5_data(self.background_file)
        s.setValue('last background value', self.background)

    def onAddImagesClicked(self):
        s = QSettings()
        self._startCalibrationBtn.setEnabled(False)
        for i in range(0, len(self.current_positions_list)):
            self.image_dict[self.current_positions_list[i]] = self.current_image_list[i] - s.value('last background value')
            if self.str_file_numbers == "":
                self.str_file_numbers = self.str_file_numbers + self.current_index_list[i]
            else:
                self.str_file_numbers = self.str_file_numbers + ", " + self.current_index_list[i]
        self._spectraPathLine.setText("")
        self._addedFilesEdit.setText(self.str_file_numbers)
        self.index_list.append(self.current_index_list)

        self.current_index_list = []
        self.current_image_list = []

    def onRegionSelectionClicked(self):
        self._saveEnergyPlotChk.setEnabled(False)
        self._saveFullImageChk.setEnabled(False)
        self._saveFullSpectrumChk.setEnabled(False)
        self._save2DSpectrumChk.setEnabled(False)
        self._saveSaturationChk.setEnabled(False)
        self._saveAllChk.setEnabled(False)
        self.sorted_dict = dict(sorted(self.image_dict.items()))
        for key in self.sorted_dict.keys():
            self.deltaP = np.append(self.deltaP, abs(key - self.zero_pos))
            self.image_list.append(self.sorted_dict.get(key))
        s = QSettings()
        s.setValue('last region selection img', self.image_list[-1])

        self.popUpSelectRegion = SelectRegion()
        self.popUpSelectRegion.show()
        self._startCalibrationBtn.setEnabled(True)

    def onStartCalibrationClicked(self):
        a = int(self.popUpSelectRegion.selectRegionCursor1.getYPos())
        b = int(self.popUpSelectRegion.selectRegionCursor2.getYPos())
        if a < b:
            self.y_up_sel = b
            self.y_low_sel = a
        else:
            self.y_up_sel = a
            self.y_low_sel = b
        self.popUpSelectRegion.close()
        for i in range(0, len(self.image_list)):
            image_integrated = self.image_list[i][self.y_low_sel:self.y_up_sel, :].sum(0)
            self.image_integrated_list.append(image_integrated)
        self.shifts_pixels = np.append(0, np.abs(self.deltaP[:-1] - self.deltaP[1:]) / self._sSpn.value()).astype(int)
        self.full_image_length = 1280 + sum(self.shifts_pixels)  # length of the full spectrum image in pixels
        self.full_image = np.zeros((1024, self.full_image_length))
        index = self.full_image_length - 1
        i = 0
        while i < len(self.image_list):
            index = index - self.shifts_pixels[i]
            if i == len(self.image_list) - 1:
                self.full_image[:, 0:1280] = np.maximum(self.full_image[:, 0:1280], self.image_list[i])
            else:
                self.full_image[:, index - 1280:index] = np.maximum(self.full_image[:, index - 1280:index],
                                                                    self.image_list[i])
            i = i + 1

        self.full_spectrum = np.zeros((self.full_image_length,))
        index = self.full_image_length - 1
        i = 0
        while i < len(self.image_list):
            index = index - self.shifts_pixels[i]
            if i == len(self.image_list) - 1:
                self.full_spectrum[0:1280] = np.maximum(self.full_spectrum[0:1280], self.image_integrated_list[i])
            else:
                self.full_spectrum[index - 1280:index] = np.maximum(self.full_spectrum[index - 1280:index], self.image_integrated_list[i])
            i = i + 1

        energy_calibration = GratingSpectrometerCalibration()
        energy_calibration.setx0(self._x0Spn.value())
        energy_calibration.setS(self._sSpn.value())
        energy_calibration.setAlpha0(self.alpha0)
        energy_calibration.setDeltaP(self.deltaP)
        energy_calibration.setEnergyMin(self._energyMinSpn.value())
        energy_calibration.setEnergyMax(self._energyMaxSpn.value())
        energy_calibration.setEnergyBinWidth(self._energyBinWidthSpn.value())
        energy_calibration.setd0(self.d0)
        energy_calibration.setImage(self.full_image)
        energy_calibration.setImageIntegrated(self.full_spectrum)
        self.energy, self.spectral_intensity_val, self.spectral_image_val = energy_calibration.setGSCalibrationMultiple()
        self.energy_selection = self.energy[np.where(self.energy <= 300)]
        self.en_sel_len = len(self.energy_selection)
        self.spectral_intensity_sel = self.spectral_intensity_val[:self.en_sel_len]
        self._startCalibrationBtn.setEnabled(True)
        self.default_foldername = "SoftX_Calibration " + datetime.today().strftime('%Y-%m-%d %H.%M.%S')
        self._setFolderNameLine.setText(self.default_foldername)

        self.index_list = []
        self._addedFilesEdit.setText("")
        self.str_file_numbers = ""
        self.spectra_file_list = []
        self.image_dict = {}
        self.last_pos = self.deltaP[-1]
        self.deltaP = np.array([])
        self.image_list = []
        self.image_integrated_list = []

    def onCreateFolderClicked(self):
        s = QSettings()
        self.export_directory_location = QtWidgets.QFileDialog.getExistingDirectory(self, "Find Directory",
                                               s.value('last export directory', r'C:/',),
                                               QtWidgets.QFileDialog.ShowDirsOnly)
        s.setValue('last export directory', self.export_directory_location)
        new_path = os.path.join(self.export_directory_location, self._setFolderNameLine.text())
        default_path = os.path.join(self.export_directory_location, self.default_foldername)
        try:
            os.mkdir(new_path)
            self.root = new_path
        except FileExistsError:
            print("Folder already exists! Files saved in a new folder: ", self.default_foldername)
            os.makedirs(default_path)
            self.root = default_path
        os.startfile(self.root)
        self._saveEnergyPlotChk.setEnabled(True)
        self._saveFullImageChk.setEnabled(True)
        self._saveFullSpectrumChk.setEnabled(True)
        self._save2DSpectrumChk.setEnabled(True)
        self._saveSaturationChk.setEnabled(True)
        self._saveAllChk.setEnabled(True)
        self._saveEnergyPlotChk.setChecked(False)
        self._saveFullImageChk.setChecked(False)
        self._saveFullSpectrumChk.setChecked(False)
        self._save2DSpectrumChk.setChecked(False)
        self._saveSaturationChk.setChecked(False)
        self._saveAllChk.setChecked(False)

    def onSaveEnergyPlotToggled(self):
        if self._saveEnergyPlotChk.isChecked():
            plt.figure("Energy vs counts")
            plt.plot(self.energy, self.spectral_intensity_val)  # counts per eV
            plt.xlabel('Energy [eV]', fontsize=10)
            plt.ylabel('Counts [arb. u.]', fontsize=10)
            y_min = float(self._yMinValueEdit.text())
            y_max = float(self._yMaxValueEdit.text())
            plt.ylim((y_min, y_max))
            plt.yscale('log')
            plt.xticks(fontsize=8)
            plt.yticks(fontsize=8)
            plt.tight_layout()
            path_im = os.path.join(self.root, "energy_vs_counts_im.png")
            plt.savefig(path_im, dpi=300)
            plt.close()
            filename = datetime.today().strftime('%Y-%m-%d %H.%M.%S')
            if self._selectLowEnergyRbtn.isChecked():
                filename = filename + " LE "
            else:
                filename = filename + " HE "
            filename = filename + str(self.last_pos) + "mm"
            path_npz = os.path.join(self.root, filename)
            np.savez(path_npz, self.energy, self.spectral_intensity_val)  # saving calibration in the same format as in the Acquisition GUI
            self._saveEnergyPlotChk.setEnabled(False)

    def onSaveFullImageToggled(self):
        if self._saveFullImageChk.isChecked():
            plt.figure("Full image")
            plt.title("Full image")
            plt.imshow(self.full_image, cmap='gray')

            path_im = path = os.path.join(self.root, "full_image_im.png")
            plt.savefig(path_im, dpi=300)
            path_txt = os.path.join(self.root, "full_image_arr.txt")
            np.savetxt(path_txt, self.full_image)
            self._saveFullImageChk.setEnabled(False)

    def onSaveFullSpectrumToggled(self):
        if self._saveFullSpectrumChk.isChecked():
            plt.figure("Full spectrum")
            plt.title("Full spectrum")
            plt.plot(self.full_spectrum)
            path_im = os.path.join(self.root, "full_spectrum_im.png")
            plt.savefig(path_im, dpi=300)
            plt.close()
            path = os.path.join(self.root, "full_spectrum_arr.txt")
            np.savetxt(path, self.full_spectrum)
            self._saveFullSpectrumChk.setEnabled(False)

    def onSave2DSpectrumToggled(self):
        if self._save2DSpectrumChk.isChecked():
            plt.figure("2D spectrum plot")
            plt.xlabel('Energy [eV]', fontsize=10)
            plt.ylabel('Pixels', fontsize=10)
            plt.pcolormesh(self.energy_selection, np.arange(self.y_low_sel, self.y_up_sel),
                           self.spectral_image_val[self.y_low_sel:self.y_up_sel, :self.en_sel_len], cmap='gray', shading='auto')
            plt.colorbar()
            plt.xticks(fontsize=8)
            plt.yticks(fontsize=8)
            plt.tight_layout()
            path_im = os.path.join(self.root, "2D_spectrum_im.png")
            plt.savefig(path_im, dpi=300)
            plt.close()
            self._save2DSpectrumChk.setEnabled(False)

    def onSaveSaturationToggled(self):
        if self._saveSaturationChk.isChecked():
            num_shots = 200
            mask = np.where(self.full_image >= 1023 * num_shots, 1, 0)
            plt.figure("Saturation check")
            plt.title('Saturation check')
            plt.imshow(self.full_image, cmap='gray')
            plt.imshow(mask, cmap='Purples', alpha=0.5)
            path_im = os.path.join(self.root, "saturation_check_im.png")
            plt.savefig(path_im, dpi=300)
            plt.close()
            self._saveSaturationChk.setEnabled(False)

    def onSaveAllToggled(self):
        if self._saveAllChk.isChecked():
            if self._saveEnergyPlotChk.isEnabled():
                self._saveEnergyPlotChk.setChecked(True)
                self.onSaveEnergyPlotToggled()
                self._saveEnergyPlotChk.setEnabled(False)

            if self._saveFullImageChk.isEnabled():
                self._saveFullImageChk.setChecked(True)
                self.onSaveFullImageToggled()
                self._saveFullImageChk.setEnabled(False)

            if self._saveFullSpectrumChk.isEnabled():
                self._saveFullSpectrumChk.setChecked(True)
                self.onSaveFullSpectrumToggled()
                self._saveFullSpectrumChk.setEnabled(False)

            if self._save2DSpectrumChk.isEnabled():
                self._save2DSpectrumChk.setChecked(True)
                self.onSave2DSpectrumToggled()
                self._save2DSpectrumChk.setEnabled(False)

            if self._saveSaturationChk.isEnabled():
                self._saveSaturationChk.setChecked(True)
                self.onSaveSaturationToggled()
                self._saveSaturationChk.setEnabled(False)

            self._saveAllChk.setEnabled(False)
        else:
            self._saveEnergyPlotChk.setChecked(False)
            self._saveFullImageChk.setChecked(False)
            self._saveFullSpectrumChk.setChecked(False)
            self._save2DSpectrumChk.setChecked(False)
            self._saveSaturationChk.setChecked(False)


class ZeroWindow(QtWidgets.QMainWindow):
    """Defines a pop-up window for fitting the 0th order image"""
    def __init__(self):
        super().__init__()
        self.setGeometry(50, 400, 600, 600)
        self.setWindowTitle("0th order")
        self.layout = QtWidgets.QVBoxLayout()

        self.zeroPlot = pg.PlotWidget(self)
        self.setCentralWidget(self.zeroPlot)
        self.zeroPlot.setXRange(100, 1000, padding=0)

        self.zeroOrderCursor1 = pg.InfiniteLine(movable=True, pos=200)
        self.zeroPlot.addItem(self.zeroOrderCursor1)

        self.zeroOrderCursor2 = pg.InfiniteLine(movable=True, pos=800)
        self.zeroPlot.addItem(self.zeroOrderCursor2)


class SelectRegion(QtWidgets.QMainWindow):
    """Defines a pop-up window for selecting the integration region using the
    spectrum image closest to the 0th order out of the ones that are provided."""
    def __init__(self):
        super().__init__()
        self.setGeometry(50, 400, 640, 512)
        self.setWindowTitle("Integration region selection")
        self.layout = QtWidgets.QVBoxLayout()

        self.selectRegion = pg.PlotWidget(self)
        self.setCentralWidget(self.selectRegion)
        s = QSettings()
        image_arr = s.value('last region selection img')
        image = pg.ImageItem(np.fliplr(image_arr.T))
        self.selectRegion.addItem(image)

        self.selectRegionCursor1 = pg.InfiniteLine(movable=True, pos=200, angle=0)
        self.selectRegion.addItem(self.selectRegionCursor1)

        self.selectRegionCursor2 = pg.InfiniteLine(movable=True, pos=800, angle=0)
        self.selectRegion.addItem(self.selectRegionCursor2)


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    QtWidgets.QApplication.setOrganizationName('DESY CFEL FS-ATTO')
    QtWidgets.QApplication.setOrganizationDomain('atto.cfel.de')
    QtWidgets.QApplication.setApplicationName('SoftX Calibration')
    QResource.registerResource("resources.rcc")
    window = SoftXCalibrationGUI()
    window.show()
    sys.exit(app.exec_())
