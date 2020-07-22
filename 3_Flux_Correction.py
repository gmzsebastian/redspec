from pyraf.iraf import twodspec
from pyraf.iraf import standard
from pyraf.iraf import sensfunc
from pyraf.iraf import calibrate
import glob
import os
import subprocess

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

def create_standard_sens(directory, objecto, iraf_name, iraf_directory = 'iidscal', suffix = ''):
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
    standard_file = glob.glob("%s/%s*BiasFlatSkyOutWave.fits"%(directory, objecto))[0]

    # Check that the files don't exist
    if check_existence(iraf_name+'_flux'+suffix+'.fits', 'create_standard_sens'):
        return

    # Run the standard
    standard(input       = standard_file,             # Input name of the image
             output      = iraf_name+'_flux'+suffix,  # Output flux file used by Sensfunc
             samestar    = 'yes',                     # Same star in all apertures?
             beam_switch = 'no',                      # Beam switch spectra?
             apertures   = '',                        # Aperture selection list
             bandwidth   = 'INDEF',                   # Bandpass widths, if INDEF use the default values in the standard calibration file
             bandsep     = 'INDEF',                   # Bandpass separation, if INDEF use the default values in the standard calibration file
             fnuzero     = '3.68E-20',                # Absolute flux zero point
             extinction  = '/iraf/iraf/noao/lib/onedstds/%s/%s.dat'%(iraf_directory, iraf_name), # Extinction file
             caldir      = '/iraf/iraf/noao/lib/onedstds/%s/'%iraf_directory,                    # Directory containing calibration data
             observatory = ')_.observatory',          # Observatory for data
             interact    = 'yes',                     # Graphic interaction to define new bandpasses
             graphics    = 'stdgraph',                # Graphics output device
             cursor      = '',                        # Graphics cursor input
             star_name   = iraf_name,                 # Star name in calibration list
             airmass     = '',                        # Airmass
             exptime     = '',                        # Exposure time (seconds)
             mag         = '',                        # Magnitude of Star
             magband     = 'U',                       # Magnitude type
             teff        = '',                        # Effective temperature or spectral type
             answer      = 'y',                       # no, yes, NO, YES, NO!, YES!
             mode        = 'ql',                      # IRAF mode
             )

    # Create the sensitivity file
    sensfunc(standards     = iraf_name+'_flux'+suffix,  # Input standard star data file
             sensitivity   = 'sens_'+iraf_name+suffix,  # Output root sensitivty function imagename
             apertures     = '',                        # Aperture selection list
             ignoreaps     = 'yes',                     # Ignore apertures and make one sensitivity function?
             logfile       = 'logfile',                 # Output log for statistics information
             extinction    = '',                        # Extinction file
             newextinction = '',                        # Output revised extinction file
             observatory   = ')_.observatory',          # Observatory for data
             function      = 'legendre' ,               # Fitting function
             order         = '6',                       # Order to fit
             interactive   = 'yes',                     # Determine sensitivity function interactively?
             graphs        = 'sr',                      # Graphs per frame
             marks         = 'plus cross box',          # Data mark types
             colors        = '2 1 3 4',                 # Colors
             cursor        = '',                        # Graphics cursor input
             device        = 'stdgraph',                # Graphics output device
             answer        = 'yes',                     # no, yes, NO, YES
             mode          = 'ql',                      # IRAF mode
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

    # Check that the files don't exist
    if check_existence("%s/%s*OutWaveStd.fits"%(directory, objecto), 'iraf_standard'):
        return

    # Create file suffix
    suffix = sensfile[5:-5]

    # For each of the science targets, correct the flux
    list_of_objects = glob.glob('%s/%s*BiasFlatSkyOutWave.fits'%(directory, objecto))
    for file_name in list_of_objects:
        calibrate(input       = file_name, # Input spectra to calibrate
                  output      = file_name[:-5] + "Std_%s.fits"%suffix, # Output calibrated spectra
                  extinct     = 'no', # Apply extinction correction
                  flux        = 'yes', # Apply flux calibration
                  extinction  = '', # Extinction file
                  observatory   = ')_.observatory', # Observatory for data
                  ignoreaps   = 'yes', # Ignore aperture numbers in flux calibration
                  sensitivity = sensfile, # Image root name for sensitivity spectra.
                  fnu         = 'no', # Create spectra having units of FNU?
                  airmass     = '', #
                  exptime     = '', #
                  mode        = 'ql' # IRAF mode
                  )

def example():
    '''
    Example of how to use this script
    '''
    # Create the standard sens file for a standard star
    create_standard_sens('LTT3218' , 'LTT3218' , 'l3218', iraf_directory = 'ctionewcal')

    # Correct the target with the standard
    iraf_standard('AT2018lfe', 'AT2018lfe', 'sens_l2415.fits')

#create_standard_sens('Feige110' , 'spec' , 'feige110', iraf_directory = 'iidscal')

#iraf_standard('AT2019itq', 'spec', 'sens_feige110.fits')

