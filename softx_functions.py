import numpy as np
from math import pi, tan

"""
energy_conversion(image_integrated, image, y0, x0, s, deltaP, energy_bin_width, alpha0, d)
-----------------------
Args:
(float) d: grating period
(float) alpha0: input grating angle in radians (nominal is 3 --> shift correction seeing Al filter)
(float) s: mm/px conversion
(float) x0: center of the 0th order in pixels
(float) y0: 0th position with respect to the tangent to the grating
(float array) deltaP: consecutive positions of spectrum images that are taken into account with respect to 0th order in descending order
(float) energy_bin_width: width of the energy bins in eV
(array) image: full image acquired by overlapping of the individual images (1024x(1280 + image_shifts.sum))
(array) image_integrated: full spectrum acquired by adding the partial spectra (which are in turn acquired by integrating images in the range given by user)
-----------------------
Returns:
(array) energy: energy axis with the center of each energy bin; used on x-axis in plots
(array) spectral_intensity_val: counts per eV
(2d array) spectral_image_val: energy vs pixels
"""


def energy_conversion(d, alpha0, s, x0, y0, deltaP, energy_min, energy_max, energy_bin_width, image, image_integrated):
    # parameters (global)
    c = 299792458e3  # [mm/s]
    h_Planck = 6.62606957e-34  # [J s]
    eV = 1.602176565e-19  # [J]
    r1 = 563.2  # [mm] distance of the focal plane

    n = len(image_integrated)
    x = np.arange(1, n + 1)  # [pixels]
    # conversion for each pixel --> array, element-wise
    y = y0 + deltaP[-1] + (x - x0) * s  # - 9.4  # [mm]

    # edges separating energy bins
    energy_edges = np.arange(energy_min, energy_max, energy_bin_width)  # [eV]

    thetaF = np.arcsin(np.sin(pi / 2 - alpha0) - h_Planck * c / energy_edges / eV / d)
    first_valid = np.argmax(thetaF >= 0)
    if not thetaF[first_valid] >= 0:
        raise ValueError('No valid angles for chosen energy range.')
    thetaF = thetaF[first_valid:]
    energy_edges = energy_edges[first_valid:]

    alphaF = pi / 2 - thetaF
    yF = r1 * np.tan(alphaF)

    spectral_intensity = np.full(len(yF) - 1, np.NaN)
    image_energy = np.zeros([np.size(image, 0), len(yF) - 1])
    avg_pxs_per_bin = 0
    number_of_values = 0
    edges = []
    pixels_number_per_bin = []
    for i in range(len(yF) - 1):
        # Compare to y-array and find start and end indices to integrate
        start_distance_from_left_image_edge = yF[i + 1] - y[0]  # [mm]
        end_distance_from_left_image_edge = yF[i] - y[0]  # [mm]
        start = start_distance_from_left_image_edge / s  # [px]
        end = end_distance_from_left_image_edge / s  # [px]
        start_floor = int(start)
        end_floor = int(end)
        if end < start:
            raise ValueError('Unexpected slope of function, got end smaller than start.')
        if start < 0:
            break
        if abs(end) > len(image_integrated):
            continue
        # TODO: Think about physical units and variable pixel size. Should we multiply .sum() by s or not?
        # instead of s, shouldn't be the number of pixels in the bin (i.e. start-end)?
        # s is considered already in the definition of y and
        # in addition this loop relies on counting the pixels (start and end ended to be indices and so pixels numbers)
        # spectral_intensity[i] = big_rot_integrated[start:end].sum() / energy_bin_width  # counts per eV (sort of density)
        spectral_intensity[i] = (image_integrated[start_floor:end_floor].sum()
                                 + (end - end_floor) * image_integrated[end_floor]
                                 - (start - start_floor) * image_integrated[start_floor]) / energy_bin_width
        image_energy[:, i] = (image[:, start_floor:end_floor].sum(1)
                              + (end - end_floor) * image[:, end_floor]
                              - (start - start_floor) * image[:, start_floor]) / energy_bin_width
        number_of_values = number_of_values + 1
        pixels_number_per_bin.append(end - start)
        edges.append(start)  # edges in pixels (i.e. indices of the axis)
        avg_pxs_per_bin = avg_pxs_per_bin + (end - start)

        nan_indices = np.isnan(spectral_intensity)
        spectral_intensity_val = spectral_intensity[~ nan_indices]
        energy = energy_edges[0:-1] + 0.5 * energy_bin_width # energy axis with the center of each bin
        energy = energy[~ nan_indices]
        spectral_image_val = image_energy[:, ~ nan_indices]

    return energy, spectral_intensity_val, spectral_image_val
