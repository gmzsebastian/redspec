import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy import optimize
import glob
from numpy.polynomial.legendre import Legendre
from pyraf import iraf
import os
from matplotlib.widgets import Slider, Button
import subprocess
import argparse


def check_existence(file_name, function, verbose=True):
    '''
    Check if some files with file_name already exist.
    If they exist return True, if they don't return False.
    Print the name of the function too.

    Parameters
    -------------
    file_name: name of the files to search.
    function: Name of the function, only for plotting purposes
    verbose: Print the name of the function?

    Output
    -------------
    True if the file already exists, False if it doesn't

    '''

    # Check that the files don't exist
    exists = subprocess.Popen(
        'ls ' + file_name, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = exists.communicate()

    # If files were returned:
    if out != b'':
        if verbose:
            print("%s -- %s already exists, skipping." % (function, file_name))
        return True
    else:
        return False


def gaussian(x, h, c, w):
    '''
    A Simple gaussian function
    '''
    return h * np.exp(-(x - c)**2 / (2 * w**2))

# Define an error function to fit the gaussian


def errfunc_gauss(p, x, y): return (gaussian(x, *p) - y)**2


def normalize(x, mina, maxa):
    '''
    Normalization function to make the maximum = 1 and the minimum = -1
    '''
    return (2.0 * x - (maxa + mina)) / (maxa - mina)


def wavelength_solution(directory, objecto, cenwave, pix_scale, max_wave=9999999, min_wave=0, width_guess=4, upper_sigma=2.0, lower_sigma=2.0, order=3, max_delta=10, fit_range=9, max_width=5, arc_name='arc'):
    '''
    Calculate the wavelength solution to a spectra. The function needs a very good initial guess of the
    slope and intercept of a linear solution to the wavelength correction. It will also check if the
    file already exists.

    Parameters
    -------------
    directory: Directory where all the science files and the lamp file are located
    objecto  : Name of the object to correct for
    cenwave  : Best guess of the central wavelength in Angstroms
    pix_scale: Best guess of the pixel scale in Angstroms / pixel
    max_wave : Maximum wavelength to fit in Angstroms
    min_wave : Minimum wavelength to fit in Angstroms
    width_guess: Guess of the width of the lines in Angstroms
    upper_sigma: Used for sigma clipping of solution fit
    lower_sigma: Used for sigma clipping of solution fit
    order : Order of the Legendre fit, can be 2, 3, 4
    max_delta: Maximum separation in Angstroms between guess and result of center wavelength
    fit_range: Fit plus/minus these many Angstroms around the central wavelength
    max_width: Maximum width in Angstroms of the resulting fit of a line



    Output
    -------------
    Function will append the name of the lamp to its respective science file.
    And output a database file with the parameters of the Legendre fit.
    '''

    # Find the files with the lamp and open the first one
    input_lamp = glob.glob(directory + '/%s*BiasFlatOut.fits' % arc_name)[0]
    Lamp = fits.open(input_lamp)

    # Extract the X data from the lamp
    counts_y = Lamp[0].data[0][0]
    pixel_x = np.arange(len(counts_y)) + 1

    # Apply the approximate wavelength correction
    x_shift = cenwave - len(pixel_x) * pix_scale / 2.0
    approx_wave = pixel_x * pix_scale + x_shift

    # Import the known wavelengths
    CuArNeLines_all = np.genfromtxt("HeNeAr.dat")

    # Select only the lines that are within the range of the spectrum
    CuArNeLines_most = CuArNeLines_all[np.where(
        (min(approx_wave) < CuArNeLines_all) & (max(approx_wave) > CuArNeLines_all))]

    # Crop the lines to only the range specified by max_wave and min_wave
    CuArNeLines = CuArNeLines_most[np.where(
        (CuArNeLines_most < max_wave) & (CuArNeLines_most > min_wave))]

    # Empty variables for the parameters of the gaussian fits
    good_centers = []
    good_truths = []
    good_widths = []

    # For every line in the CuArNe file fit a gaussian around the center wavelength.
    for Line in CuArNeLines:

        # Slice the data to only fit a region of 9 angstroms to each side of the line
        Xdata = approx_wave[np.where(
            (approx_wave > Line - fit_range) & (approx_wave < Line + fit_range))]
        Ydata = counts_y[np.where(
            (approx_wave > Line - fit_range) & (approx_wave < Line + fit_range))]

        # Make an initial guess to the gaussian to fit
        h_guess = np.max(Ydata) - np.min(Ydata)     # Line Height
        c_guess = Line                              # Line Center
        w_guess = width_guess                       # Line Width
        o_guess = np.min(Ydata)
        FirstGuess = [h_guess, c_guess, w_guess]

        # Do a quick least squares fit to the data
        Optimized, success = optimize.leastsq(
            errfunc_gauss, FirstGuess[:], args=(Xdata, Ydata))

        # Rename the output parameters to something readable
        output_height = Optimized[0]
        output_center = Optimized[1]
        output_width = Optimized[2]
        delta_wavelength = np.abs(Line - output_center)

        # Accept the line only if:
        # It's less than max_width Angstroms wide
        # It's less than max_delta Angstroms from the initial guess
        # It's height is positive
        if (delta_wavelength < max_delta) and (output_width < max_width) and (output_height > 0):
            good_centers = np.append(good_centers, output_center)
            good_truths = np.append(good_truths, Line)
            good_widths = np.append(good_widths, output_width)

        # Plot the individual line gaussian function
        Xdata_out = np.linspace(min(Xdata), max(Xdata), 100)
        plt.plot(Xdata_out, gaussian(Xdata_out, *Optimized),
                 linewidth=2, color='g')
    print(str(len(good_centers)) + " Lines Found")

    # Plot the gaussian fits on top of the approximate solution for diagnosis
    plt.plot(approx_wave, counts_y, linewidth=1, color='b')
    for i in range(len(good_truths)):
        plt.axvline(x=good_truths[i], color='k', linestyle='--', label='True')
    for i in range(len(good_centers)):
        plt.axvline(x=good_centers[i], color='g',
                    linestyle='-', label='Center')
    plt.show()

    # Fit a 2nd order polynomial to the centers (center of lines in spectra determined from fits) and the
    # good_truths (lines in the input HeNeAr file that were deemed acceptable)
    # Fit the Line centers to true centers
    Parameters = np.polyfit(good_centers, good_truths, 2)
    # Function with best fit parameters
    Function = np.poly1d(Parameters)
    # Resulting function with true centers
    Model = Function(good_centers)
    Output = good_truths - Model                       # Data - Model

    # Do Sigma Clipping: Mark as Bad the points that are higher than some sigma.
    hi_simga = np.average(Output) + upper_sigma * np.std(Output)
    lo_sigma = np.average(Output) - lower_sigma * np.std(Output)
    bad_lines = np.where((Output >= hi_simga) | (Output <= lo_sigma))
    good_lines = np.where((Output < hi_simga) & (Output > lo_sigma))

    # Remove the clipped points from the list of good lines
    best_centers = good_centers[good_lines]

    # Find minimum and maximum of the centers in the original pixel array
    back_to_pixels = (best_centers - (cenwave - len(pixel_x) *
                                      pix_scale / 2.0)) / pix_scale
    #back_to_pixels = (best_centers - x_shift) / pix_scale
    Xmin = min(back_to_pixels)
    Xmax = max(back_to_pixels)

    # Normalize the good centers so that they go from -1 to 1
    Center_min = min(best_centers)
    Center_max = max(best_centers)
    model_x = normalize(best_centers, Center_min, Center_max)

    # Do the fit again with an Nth order Legendre polynomial
    Parameters_Good = np.polynomial.legendre.legfit(
        model_x, good_truths[good_lines], order)  # Fit the line centers to the true centers
    # Resulting function with true centers
    Model_Good = Legendre(Parameters_Good)(model_x)
    # Model - Data
    Output_Good = Model_Good - good_truths[good_lines]

    # Plot results to make sure everything looks ok
    plt.subplots_adjust(hspace=0.0)
    plt.subplot(211)
    plt.xlim(min(good_centers), max(good_centers))
    plt.scatter(good_centers, good_truths - good_centers, color='b')
    plt.scatter(good_centers[bad_lines], good_truths[bad_lines] -
                good_centers[bad_lines], color='r')
    plt.plot(best_centers, Model_Good - best_centers, color='g')

    plt.subplot(212)
    plt.axhline(color='g')
    plt.xlabel("Wavelength")
    plt.ylabel("Residuals")
    plt.xlim(min(good_centers), max(good_centers))
    plt.scatter(best_centers, good_truths[good_lines] - Model_Good, color='b')
    plt.show()

    # Calculate the error and resolution on the wavelength correction
    Error = np.std(Output_Good)
    true_widths = good_widths[good_lines] * 2 * np.sqrt(2 * np.log(2))
    Resolution = np.average(true_widths, weights=1 / Output_Good ** 2)
    Resolution_std = np.std(true_widths)
    print("Error = " + str(Error))
    print("Resolution = " + str(Resolution))

    # Plot Resolution of Lines
    #plt.errorbar(np.arange(len(true_widths)), true_widths, Output_Good, fmt = '.')
    #plt.axhline(y = Resolution, color = 'k')
    #plt.axhline(y = Resolution + Resolution_std, color = 'k', linestyle = '--')
    #plt.axhline(y = Resolution - Resolution_std, color = 'k', linestyle = '--')
    # plt.show()

    # Do the Final Correction on the entire spectra
    Lamp_Wave_Out = Legendre(Parameters_Good)(
        normalize(approx_wave, Center_min, Center_max))

    # Plot the output before and after the correction
    plt.plot(approx_wave, counts_y, color='k', alpha=0.2, linestyle='--')
    plt.plot(Lamp_Wave_Out, counts_y)
    for i in CuArNeLines_all:
        plt.axvline(x=i, color='r')
    plt.show()

    # Save the output to the file that will be used by IRAF
    if order == 2:
        Line0 = "# std = " + str(Error) + "\n"
        Line1 = "begin identify Lamp.ms - Ap 1" + "\n"
        Line2 = "coefficients 7" + "\n"
        Line3 = "2.0" + "\n"   # 2 = Legendre Polynomial
        Line4 = "3.0" + "\n"   # 3 = Order of the fit
        Line5 = str(Xmin) + "\n"
        Line6 = str(Xmax) + "\n"
        Line7 = str(Parameters_Good[0]) + "\n"
        Line8 = str(Parameters_Good[1]) + "\n"
        Line9 = str(Parameters_Good[2]) + "\n"
        Output = Line0 + Line1 + Line2 + Line3 + \
            Line4 + Line5 + Line6 + Line7 + Line8 + Line9
    elif order == 3:
        Line0 = "# std = " + str(Error) + "\n"
        Line1 = "begin identify Lamp.ms - Ap 1" + "\n"
        Line2 = "coefficients 8" + "\n"
        Line3 = "2.0" + "\n"   # 2 = Legendre Polynomial
        Line4 = "4.0" + "\n"   # 4 = Order of the fit
        Line5 = str(Xmin) + "\n"
        Line6 = str(Xmax) + "\n"
        Line7 = str(Parameters_Good[0]) + "\n"
        Line8 = str(Parameters_Good[1]) + "\n"
        Line9 = str(Parameters_Good[2]) + "\n"
        Line10 = str(Parameters_Good[3]) + "\n"
        Output = Line0 + Line1 + Line2 + Line3 + Line4 + \
            Line5 + Line6 + Line7 + Line8 + Line9 + Line10
    elif order == 4:
        Line0 = "# std = " + str(Error) + "\n"
        Line1 = "begin identify Lamp.ms - Ap 1" + "\n"
        Line2 = "coefficients 9" + "\n"
        Line3 = "2.0" + "\n"   # 2 = Legendre Polynomial
        Line4 = "5.0" + "\n"   # 5 = Order of the fit
        Line5 = str(Xmin) + "\n"
        Line6 = str(Xmax) + "\n"
        Line7 = str(Parameters_Good[0]) + "\n"
        Line8 = str(Parameters_Good[1]) + "\n"
        Line9 = str(Parameters_Good[2]) + "\n"
        Line10 = str(Parameters_Good[3]) + "\n"
        Line11 = str(Parameters_Good[4]) + "\n"
        Output = Line0 + Line1 + Line2 + Line3 + Line4 + Line5 + \
            Line6 + Line7 + Line8 + Line9 + Line10 + Line11
    else:
        print("Order Must be 2, 3, or 4")

    # Find where the actual image name begins for saving
    filename = input_lamp[input_lamp.find('/%s' % arc_name) + 1:-5]

    # Save the output to the IRAF database
    outputname = 'database/id' + filename
    f = open(outputname, 'w')
    f.write(Output)
    f.close()
    print("Saved " + outputname)

    # Append the lamp wavelength link to the science target
    # Do that for each image
    Inputs = glob.glob('%s/%s*SkyOut.fits' % (directory, objecto))

    # Check that the files don't exist
    if check_existence('%s/%s*OutWave.fits' % (directory, objecto), 'wavelength_solution'):
        return

    for name in Inputs:
        print(name)
        iraf.hedit(images=name, fields='REFSPEC1',
                   value=filename, add='yes', update='yes')

        # Do wavelength correction for each image
        iraf.noao.dispcor(input=name,
                          output=name[:-5] + "Wave.fits",
                          linearize='no',
                          database='database',
                          w1='INDEF',
                          w2='INDEF',
                          dw='INDEF',
                          nw='INDEF')


def test_solution(directory, cenwave, pix_scale, bright_lines, width_guess=4, max_delta=10, fit_range=9, max_width=5, arc_name='HeNeAr'):
    '''
    Find an approximate solution to the wavelength by plotting and sliding
    the central wavelength and pixel scale.

    Parameters
    ---------------
    directory : what directory is the lamp in.
    cenwave   : Apprimxate central wavelength in Angstroms
    pix_scale : Approximate pixel scale in Angstroms / pixel
    bright_lines : List of bright HeNeAr lines

    Returns
    --------------
    best_pix_scale and best_cenwave found from the slider
    '''

    # Find the files with the lamp and open the first one
    print(directory + '/%s*BiasFlat.fits' % arc_name)
    input_lamp = glob.glob(directory + '/%s*BiasFlat.fits' % arc_name)[0]
    Lamp = fits.open(input_lamp)

    # Extract the X data from the lamp
    counts_y = Lamp[0].data[600]
    pixel_x = (np.arange(len(counts_y)) + 1)

    # Apply the approximate wavelength correction
    x_shift = cenwave - len(pixel_x) * pix_scale / 2.0
    approx_wave = pixel_x * pix_scale + x_shift

    # Brightest lines in a HeNeAr lamp
    CuArNeLines_all = np.array(bright_lines)

    # Function to plot the data
    def plot_data(pixel_x, pix_scale, cenwave):
        x_shift = cenwave - len(pixel_x) * pix_scale / 2.0
        return pixel_x * pix_scale + x_shift

    fig, ax = plt.subplots()
    plt.xlabel("Wavelength")
    plt.ylabel("Counts")
    # Plot bright HeNeAr lines
    for i in CuArNeLines_all:
        plt.axvline(x=i, color='r', linestyle='--', linewidth=1, alpha=0.6)
    plt.subplots_adjust(bottom=0.28)
    l, = plt.plot(approx_wave, counts_y, color='C0', lw=1)
    plt.axis([min(approx_wave), max(approx_wave), 0 -
              np.std(counts_y) * 0.2, max(counts_y) * 1.01])

    ax_scale = plt.axes([0.2, 0.15, 0.65, 0.03])
    ax_wave = plt.axes([0.2, 0.10, 0.65, 0.03])

    s_scale = Slider(ax_scale,  'Pixel Scale', pix_scale / 2.5,
                     pix_scale * 2.5, valinit=pix_scale, valfmt='%1.4f')
    s_wave = Slider(ax_wave, 'Central Wave',   cenwave - 1500,
                    cenwave + 1500, valinit=cenwave,   valfmt='%1.4f')

    # Change the value of the plot
    def update(val):
        new_pix_scale = s_scale.val
        new_cenwave = s_wave.val
        l.set_xdata(plot_data(pixel_x, new_pix_scale, new_cenwave))
    s_scale.on_changed(update)
    s_wave.on_changed(update)

    # Reset parameters
    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset')

    # Reset Button
    def reset(event):
        s_scale.reset()
        s_wave .reset()
    button.on_clicked(reset)

    plt.show()

    print('Selected values: pixel scale = %s; central wavelength = %s' %
          (s_scale.val, s_wave.val))

    # Apply the approximate wavelength correction
    x_shift_new = s_wave.val - len(pixel_x) * s_scale.val / 2.0
    approx_wave_new = pixel_x * s_scale.val + x_shift_new

    # Empty variables for the parameters of the gaussian fits
    good_centers = []
    good_truths = []
    good_widths = []

    # For every line in the CuArNe file fit a gaussian around the center wavelength.
    for Line in CuArNeLines_all:

        # Slice the data to only fit a region of 9 angstroms to each side of the line
        Xdata = approx_wave_new[np.where(
            (approx_wave_new > Line - fit_range) & (approx_wave_new < Line + fit_range))]
        Ydata = counts_y[np.where(
            (approx_wave_new > Line - fit_range) & (approx_wave_new < Line + fit_range))]

        # Make an initial guess to the gaussian to fit
        h_guess = np.max(Ydata) - np.min(Ydata)     # Line Height
        c_guess = Line                              # Line Center
        w_guess = width_guess                       # Line Width
        o_guess = np.min(Ydata)
        FirstGuess = [h_guess, c_guess, w_guess]

        # Do a quick least squares fit to the data
        Optimized, success = optimize.leastsq(
            errfunc_gauss, FirstGuess[:], args=(Xdata, Ydata))

        # Rename the output parameters to something readable
        output_height = Optimized[0]
        output_center = Optimized[1]
        output_width = Optimized[2]
        delta_wavelength = np.abs(Line - output_center)

        # Accept the line only if:
        # It's less than max_width Angstroms wide
        # It's less than max_delta Angstroms from the initial guess
        # It's height is positive
        if (delta_wavelength < max_delta) and (output_width < max_width) and (output_height > 0):
            good_centers = np.append(good_centers, output_center)
            good_truths = np.append(good_truths, Line)
            good_widths = np.append(good_widths, output_width)

        # Plot the individual line gaussian function
        Xdata_out = np.linspace(min(Xdata), max(Xdata), 100)
        plt.plot(Xdata_out, gaussian(Xdata_out, *Optimized),
                 linewidth=2, color='g')
    print(str(len(good_centers)) + " Lines Found")

    # Plot the gaussian fits on top of the approximate solution for diagnosis
    plt.plot(approx_wave_new, counts_y, linewidth=1, color='b')
    for i in range(len(good_truths)):
        plt.axvline(x=good_truths[i], color='k', linestyle='--', label='True')
    for i in range(len(good_centers)):
        plt.axvline(x=good_centers[i], color='g',
                    linestyle='-', label='Center')
    plt.show()

    # Find minimum and maximum of the centers in the original pixel array
    back_to_pixels = (good_centers - (s_wave.val - len(pixel_x) *
                                      s_scale.val / 2.0)) / s_scale.val

    # Fit that result with a linear function to find the best solution
    m_fit, b_fit = np.polyfit(back_to_pixels, good_truths, 1)
    cenwave_fit = b_fit + len(pixel_x) * m_fit / 2.0

    print('Best values: pixel scale = %s; central wavelength = %s' %
          (m_fit, cenwave_fit))
    # return m_fit, cenwave_fit
    return s_wave.val, s_scale.val


parser = argparse.ArgumentParser()
parser.add_argument("input_dir", help="Directory of raw data", type=str)
args = parser.parse_args()

obj_list = ["ASASSN-19bt", "AT2019aov", "AT2019aqv", "ZTF18acbvkwl", "LTT3218"]
for obj in obj_list:
    best_cenwave, best_pix_scale = test_solution(args.input_dir + '/' + obj,
                                                 6098,
                                                 1.89,
                                                 bright_lines=[4471.4770, 5015.6750, 6402.2460, 7032.4127, 7438.899, 8377.6070], arc_name='HeNeAr')

    wavelength_solution(args.input_dir + '/' + obj, obj,
                        best_cenwave, best_pix_scale, arc_name='HeNeAr')
