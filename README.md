# Calibration-GUI
This GUI takes spectrum images and acquisition parameters, manipulates this data and returns plots. The project was a part of the work on the development of a beamline intended to generate attosecond pulses with photon energies in the water-window region using the high-order harmonic generation process.

# Functionalities overview

Each position within the spectrum determines a certain energy, which is why it is extremely important to determine these positions correctly. Positions in the spectrum are determined with respect to the center of the 0th order (x0). To refer to the position of the spot on the image, one can use pixels, but the position of the image itself is given in millimeters. In order to determine this conversion between these two ways of expressing the position, two images of the same 0th order but different acquisition positions must be provided. If the user wants to both find the conversion factor and change the 0th order that is used for the reference, this can be done by selecting the ”Find x0 and s” button. If only the change of the reference 0th order is needed, this can be done by selecting the ”Find x0” button. 

![](Images/Fitting%20Zero.png)

The user first clicks the ”Browse” button, which opens the dialog that enables the user to choose an .h5 file containing the 0th order data. This data is then loaded, and the intensity of pixels for certain positions is extracted. The Gaussian fit of this data is then attempted in the range from 100 to 1000 pixels. Besides the class that defines the main GUI window, there is a separate class called ZeroWindow that defines a new pop-up window. Said window contains a plot of the 0th order data and two cursors for selecting the new fitting range. Upon choosing the 0th order file, a new instance of the ZeroWindow class is created.

One of the two gratings can be chosen. In order to check if grating input angle α0 is correct, one can use fitting of the intensity peaks. Function that checks the α0 angle value first finds the peaks on the acquired images, then counts them, and calculates the distance between them. Based on this distance, the α0 angle can be determined. The procedure is highly sensitive to the prominence parameter of the fitting function since this parameter determines the pixel intensity level above which peaks are to be looked for. If the signal obtained during the acquisition is weak in comparison to the background, peaks might not be detected.

Camera pixels are equally spaced to record the spectrum, so they will correspond to unequally spaced points in photon energy. For this reason, data must undergo a process of nonlinear resampling in order to find a conversion from pixel positions to energy positions. The main goal in this process is to keep integrated intensity the same in both cases. This part of the programme is adapted into a separate function that returns the value of spaced energy values and the corresponding pixel count values for each energy section. In the SoftX Calibration main window, the user can select the minimum and maximum value for the energy conversion, as well as the energy steps or energy bin widths.

Upon clicking on the ”Browse” button for the selection of the spectrum files, the user can choose .h5 files containing the spectrum data. A multiple-file section is possible, and the numbers of the chosen files are written on the edit line. Upon clicking on the ”Browse” button for the background, one can also choose the .h5 background file. It is important that all the images chosen in the first step have the same background that is chosen in the second step. All the files will be included in the program after clicking the ”Add” button. Adding more images that do not have the same background can be done after the first set of images is added using the ”Add” button by repeating the same steps. All the images are extracted from the .h5 files, and then the corresponding background is subtracted from them.

The specific region in which the provided spectra will be integrated is selected by the use of the ”Region Selection” button. Another class contained in the SoftX Calibration software is the SelectRegion class, which contains a plotted reference image and two horizontal cursors. The image taken for the reference is the last image added. Upon clicking the ”Region Selection” button, an element of the SelectRegion class is created. After this action, the ”Add” button is disabled and the ”Start” button is enabled.

<p align="center">
  <img src=Images/Integration%20region.png>
</p>

Finally, the user can choose which plots to save to a new folder with a given default name that contains the date and time of the execution. The user can also choose another name for this folder, but in the case that the folder with the chosen name already exists, the default name is taken. Upon clicking the ”Create” button, it is possible to select the place in which the new folder will be created. After the selection, a newly created folder is opened. There is also the option to save various plots to this folder. Available plots are:

- energy plot - energy [eV] VS pixel counts [arb. u.], saved both as .png file and .npz file that contains the output of the energy conversion function,

- full image - combined spectrum image, saved both as .png, and as .txt file,

- full spectrum - combined spectrum image integrated in the chosen boundaries, also saved both as .png, and as .txt file,

- 2D spectrum - plot that shows the distribution of the energies up to 300eV with the position values in the selected range, saved as .png file,

- saturation check - combined spectrum image with a mask that shows the highest intensity places in purple, saved as .png file.

There is also the option to save all the above-mentioned files by checking the ”save all” checkbox. Once a certain plot is saved, the corresponding check button is disabled until the calibration is restarted in order to avoid saving multiple identical plots.



**Note:** File 010 is example of a spectrum file, 011 of a background and 024 and 025 of 0th order.
