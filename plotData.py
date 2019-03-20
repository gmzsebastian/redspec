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

def check_existence(file_name, function, verbose = True):
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
    exists   = subprocess.Popen('ls ' + file_name, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = exists.communicate()

    # If files were returned:
    if out != b'':
        if verbose:
            print("%s -- %s already exists, skipping."%(function, file_name))
        return True
    else:
        return False

def create_individual(directory, objecto, flux_corrected = True, suffix = ''):
    '''
    Take all the reduced .fits files in the directory and 
    break them up into their individual .fits files and 
    convert them to .txt files

    Parameters
    -------------
    directory: Directory where all the files are located.
    objecto: Name of the files to process
    flux_corrected: Plot the files that are flux corrected?
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
    if check_existence('%s/%s*_single*.fits'%(directory, objecto), 'create_individual'):
        return

    # Import all the wavelength and flux corrected files
    if flux_corrected:
        Files = glob.glob("%s/%s*WaveStd*%s*.fits"%(directory, objecto, suffix))
    # If not, import the files that are wavelength corrected, but not flux corrected.
    else:
        Files = glob.glob("%s/%s*Wave*.fits"%(directory, objecto))

    # For each one of thoese, divide them into 3 text files
    for i in range(len(Files)):
        print(Files[i])

        # Divide the file into 3 fits files
        Output = Files[i][:-5] + "_single.fits"
        longslit.scopy(input     = Files[i],   # List of input spectra
                       output    = Output,     # List of output spectra
                       w1        = 'INDEF',    # Starting wavelength
                       w2        = 'INDEF',    # Ending wavelength
                       apertures = '',         # List of input apertures or columns
                       bands     = '',         # List of input bands or lines
                       beams     = '',         # List of beams or echelle orders
                       apmodulus = '0',        # Input aperture modulus (0 = none)
                       format    = 'onedspec', # Output spectra format (onedspec or multispec)
                       renumber  = 'no',       # Renumber output apertures?
                       offset    = 0,          # Output aperture number offset
                       clobber   = 'no',       # Modify existing output images?
                       merge     = 'no',       # Merge with existing output images?
                       rebin     = 'no',       # Rebin to exact wavelength region?
                       verbose   = 'no',       # Print operations?
                       mode      = 'ql'        # IRAF mode
                       )

        # Create a text file from the first fits file  with the cleaned data
        onedspec.wspectext(input=Output[:-4] + 'fits.0001.fits', output= Output[:-5] + '_opt.txt', header='no', wformat="%0.6f", mode = 'ql')
        onedspec.wspectext(input=Output[:-4] + 'fits.1001.fits', output= Output[:-5] + '_raw.txt', header='no', wformat="%0.6f", mode = 'ql')
        onedspec.wspectext(input=Output[:-4] + 'fits.2001.fits', output= Output[:-5] + '_sig.txt', header='no', wformat="%0.6f", mode = 'ql')

def generate_output(directory, objecto, flux_corrected = True, suffix = ''):
    '''
    Generate a three column ASCII file with wavelength, flux, and sigma for
    every specified file. Also generate the same output with the optimally 
    extracted and raw spectra.

    Parameters
    -------------
    directory: Directory where all the science files are located.
    objecto: Name of the files to process
    flux_corrected: Plot the files that are flux corrected?
    suffix: Optional subset of files to plot, i.e. different standards.
            If not specified, all the files will be plotted.
    
    Output
    -------------
    Three column files with (flux, wavelength, sigma)
    '''

    Inputs_optimal = glob.glob('%s/%s*%s*_single_opt.txt'%(directory, objecto, suffix))
    Inputs_raw     = glob.glob('%s/%s*%s*_single_raw.txt'%(directory, objecto, suffix))
    Inputs_sigma   = glob.glob('%s/%s*%s*_single_sig.txt'%(directory, objecto, suffix))

    if len(np.unique((len(Inputs_optimal), len(Inputs_raw), len(Inputs_sigma)))) != 1:
        print('Number of files not equal, something went wrong')
        return

    for i in range(len(Inputs_optimal)):
        # Import data
        wavelength, optimal = np.genfromtxt(Inputs_optimal[i], unpack = True)
        wavelength, raw     = np.genfromtxt(Inputs_raw[i]    , unpack = True)
        wavelength, sigma   = np.genfromtxt(Inputs_sigma[i]  , unpack = True)

        # Save the optimal extraction and raw spectra to different files
        output_optimal = (wavelength, optimal, sigma)
        output_raw     = (wavelength, raw, sigma)

        outfile_optimal = Inputs_optimal[i][:-14] + 'optimal_final.txt'
        np.savetxt(outfile_optimal, np.transpose(output_optimal), fmt='%9.9E')
        print("Saved " + outfile_optimal)

        outfile_raw = Inputs_optimal[i][:-14] + 'raw_final.txt'
        np.savetxt(outfile_raw, np.transpose(output_raw), fmt='%9.9E')
        print("Saved " + outfile_raw)

def plot_object(directory, objecto, suffix = ''):
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
    optimal_files = glob.glob("%s/*%s*%s*_optimal_final.txt"%(directory, objecto, suffix))
    raw_files     = glob.glob("%s/*%s*%s*_raw_final.txt"%(directory, objecto, suffix))
    size          = len(optimal_files)

    # Create Colors
    color=plt.cm.plasma(np.linspace(0,1,size + 2))

    # Make one plot with everything
    for i in range(len(optimal_files)):
        wavelength, optimal_flux, sigma = np.genfromtxt(optimal_files[i], unpack = True)
        wavelength, raw_flux, sigma     = np.genfromtxt(raw_files[i], unpack = True)

        plt.plot(wavelength, optimal_flux / np.average(optimal_flux), alpha = 0.7, color = 'r', label = 'Optimal Flux')
        plt.plot(wavelength, raw_flux / np.average(raw_flux),  linestyle = '--', alpha = 0.7, color = 'b', label = 'Raw Flux')
        plt.legend(loc = 'upper left')
        plt.xlabel("Wavelength")
        plt.ylabel("Normalized Flux")
        plt.title(objecto)
        plt.ylim(-1, 5)
        plt.xlim(min(wavelength), max(wavelength))
        plt.savefig(optimal_files[i][:-18] + '.jpg', dpi = 150)
        plt.clf()

    # Make an average plot
    if size == 1:
        wavelength, flux, sigma = np.genfromtxt(optimal_files[0], unpack = True)
    elif size > 1:
        for i in range(len(optimal_files)):
            wavelength_out, flux_out, sigma_out = np.genfromtxt(optimal_files[i], unpack = True)
            try:
                wavelength_array = np.vstack((wavelength_array, wavelength_out))
                flux_array       = np.vstack((flux_array      , flux_out      ))
                sigma_array      = np.vstack((sigma_array     , sigma_out     ))
            except:
                wavelength_array = np.copy(wavelength_out)
                flux_array       = np.copy(flux_out)
                sigma_array      = np.copy(sigma_out)

        wavelength = np.nanmedian(wavelength_array, axis = 0)
        flux       = np.nanmedian(flux_array,       axis = 0)
        sigma      = np.nanmedian(sigma_array,      axis = 0)

    # Plot Averaged Spectra
    plt.plot(wavelength, flux / np.average(flux), linewidth = 0.5)
    plt.xlabel("Wavelength")
    plt.ylabel("Normalized Flux")
    plt.title(objecto)
    plt.ylim(-1, 5)
    plt.xlim(min(wavelength), max(wavelength))
    plt.savefig('%s/%s_average%s.jpg'%(directory, objecto, suffix), dpi = 150)
    plt.clf()

    # Save Averaged Spectra
    array = (wavelength, flux, sigma)
    np.savetxt('%s/%s_average%s.txt'%(directory, objecto, suffix), np.transpose(array))

def molly_parameter(directory, objecto, flux_corrected = True, suffix = ''):
    '''
    Create two molly files with the information about each spectrum.

    Telescope options:
    Palomar 200in
    Wilson
    Campanas
    Lemmon
    WHT
    INT
    JKT
    UKIRT
    Kitt Peak
    AAT
    CTIO
    McDonald
    MMT
    VLT
    ANU 2.3m
    SAAO 1.9m
    NTT
    'elsewhere'

    Parameters
    -------------
    directory: Directory where all the science files are located.
    objecto: Name of the files to process
    flux_corrected: Plot the files that are flux corrected?
    suffix: Optional subset of files to plot, i.e. different standards.
            If not specified, all the files will be plotted.
    
    Output
    -------------
    headerfile and listfile of molly files
    '''

    # Import all the wavelength and flux corrected files
    if flux_corrected:
        fits_files = glob.glob("%s/%s*WaveStd*%s*.fits"%(directory, objecto, suffix))
    # If not, import the files that are wavelength corrected, but not flux corrected.
    else:
        fits_files = glob.glob("%s/%s*Wave*.fits"%(directory, objecto))

    # Text Files
    optimal_files = glob.glob("%s/*%s*%s*_optimal_final.txt"%(directory, objecto, suffix))
    raw_files     = glob.glob("%s/*%s*%s*_raw_final.txt"%(directory, objecto, suffix))

    # Header of header output file
    header_data = "Object     UTC     Date     RA     DEC     Dwell     Airmass     Equinox     JD     Day     Month     Year \n"
    header_data += "  C         D        C       D      D        D          D            D        D      I        I         I \n"

    # Empty list for list output file
    list_data = ""

    for i in range(len(optimal_files)):
        # Edges of the filename to find the filename
        a = optimal_files[i].find('/')

        # Name of the row name that molly needs
        filename   = optimal_files[i][a+1:]
        molly_name = "lasc %s %s 1 2 3 A M 0.05 \n"%(filename, i + 1)

        # Add to list of all files
        list_data += molly_name

        ##### Header output file #####
        # For every variable and every file, extract the value listed here
        input_file = fits.open(fits_files[i])
        objecto  = input_file[0].header['OBJECT']
        ut_time  = input_file[0].header['UT']
        date_obs = input_file[0].header['DATE-OBS']
        ra_raw   = input_file[0].header['RA']
        dec_raw  = input_file[0].header['DEC']
        exptime  = input_file[0].header['EXPTIME']
        airmass  = input_file[0].header['AIRMASS']
        equinox  = 2000.0
        juliand  = Time(date_obs + " " + ut_time).jd

        # Get the variables into the right format
        # Time
        hr = float(ut_time[0:2])
        mi = float(ut_time[3:5])
        se = float(ut_time[6:])
        ut_out = hr + mi / 60.0 + se / 3600.0

        # Date
        year  = date_obs[0:4]
        month = date_obs[5:7]
        day   = date_obs[8:10]
        date_out = str(day) + "/" + str(month) + "/" + str(year)

        # Right Ascension
        if ra_raw[0] == '+':
            rahr = float(ra_raw[1:3])
            rami = float(ra_raw[4:6])
            rase = float(ra_raw[7:])
            ra_out = rahr + rami / 60.0 + rase / 3600.0
        if ra_raw[0] == '-':
            rahr = float(ra_raw[1:3])
            rami = float(ra_raw[4:6])
            rase = float(ra_raw[7:])
            ra_out = np.negative(rahr + rami / 60.0 + rase / 3600.0)

        # Declination
        if dec_raw[0] == '+':
            dechr = float(dec_raw[1:3])
            decmi = float(dec_raw[4:6])
            decse = float(dec_raw[7:])
            dec_out = dechr + decmi / 60.0 + decse / 3600.0
        if dec_raw[0] == '-':
            dechr = float(dec_raw[1:3])
            decmi = float(dec_raw[4:6])
            decse = float(dec_raw[7:])
            dec_out = np.negative(dechr + decmi / 60.0 + decse / 3600.0)

        # One line for each file with all the variables in the right format
        molly_header = str(objecto) + " " + str(ut_out) + " " + str(date_out) + " " + str(ra_out) + " " + str(dec_out) + " " + str(exptime) + " " + str(airmass) + " " + str(equinox) + " " + str(juliand) + " " + str(day) + " " + str(month) + " " + str(year) + "\n"
        # Append to the big list of all headers
        header_data += molly_header

    # Save header output file
    headerfile_name = 'headerfile_%s%s.txt'%(objecto, suffix)
    headerfile = open(directory + '/' + headerfile_name, 'w')
    headerfile.write(header_data)
    print('Saved %s/headerfile_%s%s.txt'%(directory, objecto, suffix))

    # Save list output file
    listfile_name = "listfile_" + objecto + suffix + ".txt"
    listfile = open(directory + '/' + listfile_name, "w")
    listfile.write(list_data)
    print("Saved " + "listfile_" + objecto + suffix + ".txt")

    # Create molly instructions file
    size = len(optimal_files)
    instructions = 'mxpix 4500 sure\n\n@%s.txt\n\nedit 1 %s\nfile\n%s\nq\n\nhfix 1 %s MMT\n\nvbin 1 %s 101\n\n\n\n\n'%(listfile_name, size, headerfile_name, size, size)

    for i in range(len(optimal_files)):
        a = optimal_files[i].find('/')
        # Name of the row name that molly needs
        filename   = optimal_files[i][a+1:-9] + '_molly.txt'
        number     = i + 101
        instructions += 'wasc %s %s ANGSTROMS MJY\n'%(filename, number)

    instructions_name = "instructions_" + objecto + suffix + ".txt"
    instructions_file = open(directory + '/' + instructions_name, "w")
    instructions_file.write(instructions)
    print("Saved " + "instructions_" + objecto + suffix + ".txt")

def all_in_one(directory, objecto, flux_corrected = True, suffix = ''):
    os.system('mv %s A'%directory)
    create_individual('A', objecto, flux_corrected, suffix)
    generate_output('A', objecto, flux_corrected, suffix)
    plot_object('A', objecto, suffix)
    os.system('mv A %s'%directory)

def example():
    '''
    Example of how to use this script
    '''
    # Create individual files and plot them
    all_in_one('AT2018lfe', 'AT2018lfe')

    # If you wish to create a file to be read in by Molly, do:
    molly_parameter('AT2018lfe', 'AT2018lfe')
    # WARNING, this function has not been tested in whiiiiile, might break a lot
    