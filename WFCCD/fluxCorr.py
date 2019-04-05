from pyraf.iraf import twodspec
from pyraf.iraf import standard
from pyraf.iraf import sensfunc
from pyraf.iraf import calibrate
import glob
import os
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


def create_standard_sens(directory, objecto, iraf_name, iraf_directory='iidscal', suffix=''):
    '''
    Calculate the flux calibration for the standard star and save. The list of standards
    can be found in: http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?onedstds

    Parameters
    -------------
    directory: Directory where all the reduced standard files are located.
               i.e. Objectname/Science
    iraf_name: Name of the standard in the IRAF tables
    object_name: The object you intend to correct
    iraf_directory: Directory where the extinction file is. Options:
        blackbody
        bstdscal
        ctiocal
        ctionewcal
        iidscal
        irscal
        oke1990
        redcal
        spec16cal
        spec50cal
        spechayescal
    suffix: Optional suffix to be appended at the end of the file name

    Output
    -------------
    A file that will be used by SENSFUNC.
    The calibratied sens file
    The flux corrected images of the target.
    '''

    # Use the first image of the standard to calculate the calibration
    standard_file = glob.glob(
        "%s/%s*BiasFlatSkyOutWave.fits" % (directory, objecto))[0]
    standard_output = directory + '/' + iraf_name + '_flux' + suffix

    # Extinction information
    extinction_file = '/Users/pcowpert/anaconda3/envs/iraf27/iraf/noao/'\
        'lib/onedstds/%s/%s.dat' % (iraf_directory, iraf_name)
    caldir_file = '/Users/pcowpert/anaconda3/envs/iraf27/iraf/noao/'\
        'lib/onedstds/%s/' % iraf_directory

    # Check that the files don't exist
    if not check_existence(standard_output, 'create_standard_sens'):
        standard(input=standard_file,  # Input name of the image
                 output=standard_output,  # Output flux file used by Sensfunc
                 samestar='yes',  # Same star in all apertures?
                 beam_switch='no',  # Beam switch spectra?
                 apertures='',  # Aperture selection list
                 # Bandpass widths, if INDEF use the default values in the standard calibration file
                 bandwidth='INDEF',
                 # Bandpass separation, if INDEF use the default values in the standard calibration file
                 bandsep='INDEF',
                 fnuzero='3.68E-20',  # Absolute flux zero point
                 # Extinction file
                 extinction=extinction_file,
                 # Directory containing calibration data
                 caldir=caldir_file,
                 observatory=')_.observatory',  # Observatory for data
                 interact='yes',  # Graphic interaction to define new bandpasses
                 graphics='stdgraph',  # Graphics output device
                 cursor='',  # Graphics cursor input
                 star_name=iraf_name,  # Star name in calibration list
                 airmass='',  # Airmass
                 exptime='',  # Exposure time (seconds)
                 mag='',  # Magnitude of Star
                 magband='U',  # Magnitude type
                 teff='',  # Effective temperature or spectral type
                 answer='y',  # no, yes, NO, YES, NO!, YES!
                 mode='ql',  # IRAF mode
                 )

    sens_output = directory + '/sens_' + iraf_name + suffix

    if not check_existence(sens_output + '.fits', 'create_standard_sens'):
        # Create the sensitivity file
        sensfunc(  # Input standard star data file
            standards=standard_output,
            # Output root sensitivty function imagename
            sensitivity=sens_output,
            apertures='',   # Aperture selection list
            # Ignore apertures and make one sensitivity function?
            ignoreaps='yes',
            logfile='logfile',  # Output log for statistics information
            extinction=extinction_file,  # Extinction file
            newextinction='',  # Output revised extinction file
            observatory=')_.observatory',  # Observatory for data
            function='legendre',  # Fitting function
            order='8',  # Order to fit
            interactive='yes',  # Determine sensitivity function interactively?
            graphs='sr',  # Graphs per frame
            marks='plus cross box',  # Data mark types
            colors='2 1 3 4',  # Colors
            cursor='',  # Graphics cursor input
            device='stdgraph',  # Graphics output device
            answer='yes',  # no, yes, NO, YES
            mode='ql',  # IRAF mode
        )


def iraf_standard(directory, objecto, sensfile):
    '''
    Do the flux calibration on the targets in the directory using the
    sensfile

    Parameters
    -------------
    directory: Directory where all the science files are located.
               i.e. Objectname/Science
    objecto: Name of the object in the directory to flux calibrate
    sensfile: Full name of the sensfile to calibrate the data

    Output
    -------------
    Flux calibrated data
    '''

    # For each of the science targets, correct the flux
    list_of_objects = glob.glob(
        '%s/%s*BiasFlatSkyOutWave.fits' % (directory, objecto))

    # Create file suffix
    suffix = sensfile[sensfile.find('sens_'):sensfile.find('.fits')]

    for file_name in list_of_objects:

        # Check that the files don't exist
        output_file = file_name.replace('.fits', "Std_%s.fits" % suffix)
        if check_existence(output_file, 'iraf_standard'):
            continue

        calibrate(input=file_name,  # Input spectra to calibrate
                  # Output calibrated spectra
                  output=output_file,
                  extinct='no',  # Apply extinction correction
                  flux='yes',  # Apply flux calibration
                  extinction='',  # Extinction file
                  observatory=')_.observatory',  # Observatory for data
                  ignoreaps='yes',  # Ignore aperture numbers in flux calibration
                  # Image root name for sensitivity spectra.
                  sensitivity=sensfile,
                  fnu='no',  # Create spectra having units of FNU?
                  airmass='',
                  exptime='',
                  mode='ql'  # IRAF mode
                  )


parser = argparse.ArgumentParser()
parser.add_argument("input_dir", help="Directory of raw data", type=str)
args = parser.parse_args()

create_standard_sens(args.input_dir + '/LTT3218', 'LTT3218', 'l3218',
                     iraf_directory='ctionewcal')

obj_list = ["ASASSN-19bt", "AT2019aov", "AT2019aqv", "ZTF18acbvkwl"]

for obj in obj_list:
    iraf_standard(args.input_dir + '/' + obj, obj,
                  args.input_dir + '/LTT3218/sens_l3218.fits')
