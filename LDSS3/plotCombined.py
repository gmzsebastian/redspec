
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from pyraf import iraf
from pyraf.iraf import imred
from pyraf.iraf import twodspec
from pyraf.iraf import longslit
import glob
from astropy.time import Time
import subprocess
from pyraf.iraf import onedspec
import os
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


def create_individual(directory, objecto, suffix=''):
    '''
    Take all the reduced .fits files in the directory and
    break them up into their individual .fits files and
    convert them to .txt files
    Parameters
    -------------
    directory: Directory where all the files are located.
    objecto: Name of the files to process
    suffix: Optional subset of files to plot, i.e. different standards.
            If not specified, all the files will be plotted.
    Output
    -------------
    Set of .fits files with
    - Clean Spectra
    - Raw Spectra
    - Errorbars
    One .txt file with wavelength and clean spectra
    '''

    # Check that the files do not exist
    if check_existence('%s/%s_combined_single*.fits' % (directory, objecto), 'create_individual'):
        return

    Files = glob.glob("%s/%s_combined.fits" % (directory, objecto))

    # For each one of thoese, divide them into 3 text files
    for i in range(len(Files)):
        print(Files[i])

        # Divide the file into 3 fits files
        Output = Files[i][:-5] + "_single.fits"
        longslit.scopy(input=Files[i],   # List of input spectra
                       output=Output,     # List of output spectra
                       w1='INDEF',    # Starting wavelength
                       w2='INDEF',    # Ending wavelength
                       apertures='',         # List of input apertures or columns
                       bands='',         # List of input bands or lines
                       beams='',         # List of beams or echelle orders
                       # Input aperture modulus (0 = none)
                       apmodulus='0',
                       # Output spectra format (onedspec or multispec)
                       format='onedspec',
                       renumber='no',       # Renumber output apertures?
                       offset=0,          # Output aperture number offset
                       clobber='no',       # Modify existing output images?
                       merge='no',       # Merge with existing output images?
                       rebin='no',       # Rebin to exact wavelength region?
                       verbose='no',       # Print operations?
                       mode='ql'        # IRAF mode
                       )

        # Create a text file from the first fits file  with the cleaned data
        onedspec.wspectext(input=Output[:-4] + 'fits.0001.fits',
                           output=Output[:-5] + '_opt.txt',
                           header='no',
                           wformat="%0.6f",
                           mode='ql')


def generate_output(directory, objecto, suffix=''):
    '''
    Generate a three column ASCII file with wavelength, flux, and sigma for
    every specified file. Also generate the same output with the optimally
    extracted and raw spectra.
    Parameters
    -------------
    directory: Directory where all the science files are located.
    objecto: Name of the files to process
    suffix: Optional subset of files to plot, i.e. different standards.
            If not specified, all the files will be plotted.

    Output
    -------------
    Three column files with (flux, wavelength, sigma)
    '''

    Inputs_optimal = glob.glob('%s/%s_combined_single_opt.txt' %
                               (directory, objecto))

    for i in range(len(Inputs_optimal)):
        # Import data
        wavelength, optimal = np.genfromtxt(Inputs_optimal[i], unpack=True)

        # Save the optimal extraction and raw spectra to different files
        output_optimal = (wavelength, optimal)

        outfile_optimal = Inputs_optimal[i][:-14] + 'optimal_final.txt'
        np.savetxt(outfile_optimal, np.transpose(output_optimal),
                   delimiter=' ', fmt='%9.9E')
        print("Saved " + outfile_optimal)


def plot_object(directory, objecto, suffix=''):
    '''
    Plot all the spectra by reading in the three column files.
    Parameters
    -------------
    directory: Directory where all the science files are located.
    objecto: Name of the files to process
    suffix: Optional subset of files to plot, i.e. different standards.
            If not specified, all the files will be plotted.

    Output
    -------------
    Plot of spectra
    '''

    # Import all the wavelength and flux corrected text files
    optimal_files = glob.glob("%s/%s_combined_optimal_final.txt" %
                              (directory, objecto))
    size = len(optimal_files)

    # Create Colors
    color = plt.cm.plasma(np.linspace(0, 1, size + 2))

    # Make one plot with everything
    for i in range(len(optimal_files)):
        wavelength, optimal_flux = np.genfromtxt(
            optimal_files[i], unpack=True)

        plt.plot(wavelength, optimal_flux / np.average(optimal_flux),
                 alpha=0.7, color='r', label='Optimal Flux')
        plt.legend(loc='upper left')
        plt.xlabel("Wavelength")
        plt.ylabel("Normalized Flux")
        plt.title(objecto)
        #plt.ylim(-1, 5)
        plt.xlim(3800, 9600)
        plt.savefig(optimal_files[i].replace('.txt', '.jpg'), dpi=150)
        plt.clf()


def all_in_one(directory, objecto, suffix=''):
    os.system('mv %s A' % directory)
    create_individual('A', objecto, suffix)
    generate_output('A', objecto, suffix)
    plot_object('A', objecto, suffix)
    os.system('mv A %s' % directory)


parser = argparse.ArgumentParser()
parser.add_argument("input_dir", help="Directory of raw data", type=str)
args = parser.parse_args()

obj_list = ["AT2019yx"]

for obj in obj_list:
    all_in_one(args.input_dir + '/' + obj, obj)