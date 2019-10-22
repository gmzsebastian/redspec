from astropy.io import fits
from pyraf import iraf
from pyraf.iraf import imred
from pyraf.iraf import twodspec
from pyraf.iraf import apextract
from pyraf.iraf import apall
from pyraf.iraf import longslit
from pyraf.iraf import noao
from pyraf.iraf import ccdred
from pyraf.iraf import zerocombine
from pyraf.iraf import flatcombine
from pyraf.iraf import response
from pyraf.iraf import ccdproc
from pyraf.iraf import imarith
from pyraf.iraf import background
import glob
import os
import subprocess
import numpy as np

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

def create_lists_science(directory, objecto, extension):
    '''
    Create a series of lists of the science files to be read in by
    IRAF. Each list will contain every fits file in the specified
    directory.

    Parameters
    -------------
    directory: Directory where all the science files are located.
               i.e. Objectname/Science

    Output
    -------------
    5 lists:
    - List of objects
    - List + _BiasFlat
    - List + Sky
    - List + Out
    - List + Wave
    '''

    # Create lists for science files
    os.system('ls %s/%s*.fits > %s/Allfiles'%(directory, objecto, directory))
    os.system("sed 's/.fits/_BiasFlat.fits/g' %s/Allfiles > %s/AllfilesBiasFlat"%(directory, directory))
    os.system("sed 's/.fits/Sky.fits/g' %s/AllfilesBiasFlat > %s/AllfilesBiasFlatSky"%(directory, directory))
    os.system("sed 's/.fits/Out.fits/g' %s/AllfilesBiasFlatSky > %s/AllfilesBiasFlatSkyOut"%(directory, directory))
    os.system("sed 's/.fits/Wave.fits/g' %s/AllfilesBiasFlatSkyOut > %s/AllfilesBiasFlatSkyOutWave"%(directory, directory))

    if extension != '':
        os.system("sed 's/.fits/.fits[%s]/g' %s/Allfiles > %s/Allfiles%s"%(extension, directory, directory, extension))

def create_lists_flats(directory, objecto, extension):
    '''
    Create a series of lists of the flat files to be read in by
    IRAF. Each list will contain every fits file in the specified
    directory.

    Parameters
    -------------
    directory: Directory where all the bias files are located.
               i.e. Objectname/Flats

    Output
    -------------
    5 lists:
    - List of flat files
    - List + _Dark
    '''

    # Create lists for flat files
    os.system('ls %s/%s*.fits > %s/flatlist'%(directory, objecto, directory))
    os.system("sed 's/.fits/_Bias.fits/g' %s/flatlist > %s/flatbiaslist"%(directory, directory))

    if extension != '':
        os.system("sed 's/.fits/.fits[%s]/g' %s/flatlist > %s/flatlist%s"%(extension, directory, directory, extension))

def create_lists_lamp(directory, objecto, extension):
    '''
    Create a series of lists of the lamp files to be read in by
    IRAF. Each list will contain every fits file in the specified
    directory.

    Parameters
    -------------
    directory: Directory where all the lamp files are located.
               i.e. Objectname/Lamp

    Output
    -------------
    5 lists:
    - List of lamp files
    - List + _BiasFlat
    - List + Out
    '''

    # Create lists for lamp files
    os.system('ls %s/%s*.fits > %s/lamplist'%(directory, objecto, directory))
    os.system("sed 's/.fits/_BiasFlat.fits/g' %s/lamplist > %s/lamplistbiasflat"%(directory, directory))
    os.system("sed 's/.fits/Out.fits/g' %s/lamplistbiasflat > %s/lamplistbiasflatout"%(directory, directory))

    if extension != '':
        os.system("sed 's/.fits/.fits[%s]/g' %s/lamplist > %s/lamplist%s"%(extension, directory, directory, extension))

def iraf_zerocombine(directory, extension = '', datasec_key = 'DATASEC'):
    '''
    Function will create a list of bias files out of every .fits file in the
    specified directory. It will then combine all of those files into a single
    Bias.fits file using IRAF's zerocombine.

    If a file called directory/Bias.fits already exists then don't do anything.

    Parameters
    -------------
    directory: Directory where all the bias files are located.
               i.e. 'Bias'

    Output
    -------------
    One master Bias.fits file
    '''

    # Create list of bias files
    os.system("ls " + directory + "/*.fits > " + directory + "/biaslist")
    bias_list = 'biaslist'

    # Test to see if the Bias.fits file already exists
    if check_existence(directory + '/Bias.fits', 'iraf_zerocombine'):
        return

    # Add extension if specified
    if extension != '':
        os.system("sed 's/.fits/.fits[%s]/g' %s/biaslist > %s/biaslist%s"%(extension, directory, directory, extension))
        bias_list = 'biaslist' + str(extension)

    zerocombine(input   = '@%s/%s'%(directory, bias_list), # List of zero images to combine, start with @ for a list      
                output  = '%s/Bias.fits'%directory,  # Output bias level name
                combine = 'median',                  # Type of combine operation (average or median)
                reject  = 'avsigclip',               # Type of rejection algorithm
                                                     #      none - No rejection
                                                     #    minmax - Reject the nlow and nhigh pixels
                                                     #   ccdclip - Reject pixels using CCD noise parameters
                                                     #  crreject - Reject only positive pixels using CCD noise parameters
                                                     #   sigclip - Reject pixels using a sigma clipping algorithm
                                                     # avsigclip - Reject pixels using an averaged sigma clipping algorithm
                                                     #     pclip - Reject pixels using sigma based on percentiles
                ccdtype = '',                        # Tpe of CCD image to combine, if '' all images will be used
                process = 'no',                      # Process images before combining?
                delete  = 'no',                      # Delete input images after combining?
                clobber = 'no',                      # Clobber existing output image?
                scale   = 'none',                    # Scale the images? for example by exposure time
                                                     # options: none, mode, median, mean, exposure
                statsec = '',                        # Image section for computing statistics
                nlow    = '0',                       # Number of low pixels to reject
                nhigh   = '1',                       # Number of high pixels to reject
                nkeep   = '1',                       # If positive, minimum to keep.
                                                     # If negative, maximum to reject
                mclip   = 'yes',                     # Use median in sigma clipping algorithm?
                lsigma  = '3.0',                     # Lower sigma clipping factor
                hsigma  = '3.0',                     # Upper sigma clipping factor
                rdnoise = 'rdnoise',                 # CCD redout noise (electrons)
                gain    = 'gain',                    # CCD gain (electrons / DN)
                snoise  = '0.0',                     # Sensitiviy noise
                pclip   = '-0.5',                    # Percentile clipping algorithm parameter. 
                                                     # If positive: Specifies number of pixels above or below the median to use for computing the clipping sigma.
                                                     # If negative: Specifies the fraction of pixels above or below the median to use.
                                                     # The default of -0.5 selects approximately the quartile point.
                blank   = '0.0',                     # Value if there are no pixels
                mode    = 'ql'                       # IRAF mode
                )

    input_bias = np.genfromtxt('%s/%s'%(directory, bias_list), dtype = 'str')[0]
    if extension != '':
        datasec_in = fits.getval(input_bias[:-3], datasec_key,ext = extension)
    else:
        datasec_in = fits.getval(input_bias, datasec_key,ext = extension)
    fits.setval('%s/Bias.fits'%directory, datasec_key, value = datasec_in)

    print('Saved Bias.fits to %s/ \n'%directory)

def iraf_ccdproc(directory, file_type, bias_file, flat_file = '', extension = ''):
    '''
    Do the bias and flats correction for a series of fits files.
    If the files are flats, then it doesn't do a flat correction.

    Parameters
    -------------
    directory: Directory where all the files to be corrected are located.
               i.e. Objectname/Science
    bias_file: Directory with all the bias files
               i.e. bias/Bias.fits
    flat_dir: Directory with the Flat_norm file
              i.e. Objectname/Flats
    file_type: Type of file being corrected. Either 'science' or 'lamp'
               if a flat normalization will happen. Or 'flat' if there
               won't be a flat correction. 

    Output
    -------------
    A set of bias and flat corrected fits files in the 'directory'
    '''

    # Check that the files don't exist
    if check_existence('%s/%s*_BiasFlat.fits'%(directory, file_type), 'iraf_ccdproc'):
        return

    # If correcting flat files, don't do flat corrections
    if (file_type == 'flat') or (file_type == 'flats') or (file_type == 'BClamp'):
        images_in  = '@%s/flatlist%s'%(directory, extension)
        output_in  = '@%s/flatbiaslist'%directory
        flat       = ''
        flatcor_in = 'no'

    # If correcting the lamp files, do flat correction
    elif file_type in ['lamp', 'arc', 'HeNeAr', 'henear', 'arcs', 'Arc', 'HeArNe', 'COMP']:
        images_in  = '@%s/lamplist%s'%(directory, extension)
        output_in  = '@%s/lamplistbiasflat'%directory
        flat       = flat_file
        flatcor_in = 'yes'

    # If correcting sceicne files, do flat correction
    else:
        images_in  = '@%s/Allfiles%s'%(directory, extension)
        output_in  = '@%s/AllfilesBiasFlat'%directory
        flat       = flat_file
        flatcor_in = 'yes'

    # Do the correction
    ccdproc(images      = images_in,   # List of CCD images to correct
            output      = output_in,   # Name of output images
            ccdtype     = '',          # CCD image type to correct, if empty all images will be used
            max_cache   = '0',         # Maximum image cachine memory in MB
            noproc      = 'no',        # List processing steps only?
            fixpix      = 'no',        # Fix bad CCD lines and columns?
            overscan    = 'no',        # Apply overscan strip correction?
            trim        = 'no',        # Trim the image?
            zerocor     = 'yes',       # Apply zero level correction?
            darkcor     = 'no',        # Apply dark count correction?
            flatcor     = flatcor_in,  # Apply flat field correction?
            illumcor    = 'no',        # Apply illumination correction?
            fringecor   = 'no',        # Apply fringe correction?
            readcor     = 'no',        # Convert zero level image to readout correction?
                                       # If yes, zero images are averaged across readout axis
            scancor     = 'no',        # Convert flat field image to scan correction?
                                       # If yes, the form of scan mode correction is specified by scantype
            readaxis    = 'line',      # Read out axis (Column or line)
            fixfile     = '',          # File describing the bad lines and columns
            biassec     = '',          # Overscan strip image section
            trimsec     = '',          # Trim data section
            zero        = bias_file,   # Zero level calibration image
            dark        = '',          # Dark count calibration image
            flat        = flat,        # Flat field images
            illum       = '',          # Illumination correction images
            fringe      = '',          # Fringe correction images
            minreplace  = '1.0',       # Minimum flat filed value, values below this value
                                       # are replaced by this value. To avoid divide by 0
            scantype    = 'shortscan', # Scan type (Shortscan or longscan), longscan is a mode
                                       # where the CCD is clocked and read out continuously.
            nscan       = '1',         # Number of short scan lines
            interactive = 'no',        # Fit overscan interactively?
            function    = 'legendre',  # Fitting function
            order       = '1',         # Number of polynomial terms or spline pieces
            sample      = '*',         # Sample points to fit
            naverage    = '1',         # Number of sample points to combine
            niterate    = '3',         # Number of rejection iterations
            low_reject  = '3.0',       # Low sigma rejection factor
            high_reject = '2.0',       # High sigma rejection factor
            grow        = '0.0',       # Rejection growing radius
            mode        = 'ql',        # IRAF mode
            )

def iraf_flatcombine(directory, response_sample, fit_order):
    '''
    Merge all the flat files into one Flat.fits and then calculate the 
    response function for that file a generate a Flat_norm.fits file. 

    Parameters
    -------------
    directory: Directory where all the flat files are located.
               i.e. Objectname/Flats
    response_sample: X-axis from which to sample and fit the response
                     function
    fit_order: Order of the fit to the response function.

    Output
    -------------
    - A Flat.fits file with the median of all the flat files
    - Flat_norm.fits that is the response function of this file.
    '''

    # Check that the files don't exist
    if check_existence('%s/Flat_norm.fits'%directory, 'iraf_flatcombine'):
        return

    # Run flatcombine to create the average flat file
    flatcombine(input   = '@%s/flatbiaslist'%directory,   # List of flat field images to combine
                output  = '%s/Flat.fits'%directory,       # Output name for the flat field
                combine = 'median',                       # Type of combine operation (average or median)
                reject  = 'avsigclip',                    # Type of rejection algorithm
                                                          #      none - No rejection
                                                          #    minmax - Reject the nlow and nhigh pixels
                                                          #   ccdclip - Reject pixels using CCD noise parameters
                                                          #  crreject - Reject only positive pixels using CCD noise parameters
                                                          #   sigclip - Reject pixels using a sigma clipping algorithm
                                                          # avsigclip - Reject pixels using an averaged sigma clipping algorithm
                                                          #     pclip - Reject pixels using sigma based on percentiles
                ccdtype = '',                             # Tpe of CCD image to combine, if '' all images will be used
                process = 'no',                           # Process images before combining?
                subsets = 'no',                           # Combine different images into groups with different subset parameters.
                delete  = 'no',                           # Delete input images after combining?
                clobber = 'no',                           # Clobber existing output image?
                scale   = 'none',                         # Scale the images? for example by exposure time
                                                          # options: none, mode, median, mean, exposure
                statsec = '',                             # Image section for computing statistics
                nlow    = '0',                            # Number of low pixels to reject
                nhigh   = '1',                            # Number of high pixels to reject
                nkeep   = '1',                            # If positive, minimum to keep.
                                                          # If negative, maximum to reject
                mclip   = 'yes',                          # Use median in sigma clipping algorithm?
                lsigma  = '3.0',                          # Lower sigma clipping factor
                hsigma  = '3.0',                          # Upper sigma clipping factor
                rdnoise = 'rdnoise',                      # CCD redout noise (electrons)
                gain    = 'gain',                         # CCD gain (electrons / DN)
                snoise  = '0.0',                          # Sensitiviy noise
                pclip   = '-0.5',                         # Percentile clipping algorithm parameter. 
                                                          # If positive: Specifies number of pixels above or below the median to use for computing the clipping sigma.
                                                          # If negative: Specifies the fraction of pixels above or below the median to use.
                                                          # The default of -0.5 selects approximately the quartile point.
                blank   = '1.0',                          # Value if there are no pixels
                mode    = 'ql'                            # IRAF mode
                )

    # Determine response calibrations
    response(calibration   = '%s/Flat.fits'%directory,        # Images to use in determining response calibrations.
             normalization = '%s/Flat.fits'%directory,        # Images to use in determining the normalization spectrum, usually the same as calibration.
             response      = '%s/Flat_norm.fits'%directory,   # Output image name
             interactive   = 'yes',                           # Graph and fit interactively? yes or no
             threshold     = 'INDEF',                         # Set the response to 1 when the normalization is below this threshold.
                                                              # If INDEF then don't do anything.
             sample        = response_sample,                 # Samples of points to use in fitting the average function.
                                                              # If * all points will be used
             naverage      = '1',                             # Number of sample points to average or median before fitting the function.
                                                              # Positive is average, negative is median.
             function      = 'legendre',                      # Function to fit the average, either 'spline1', 'spline3', 'legendre', or 'chebyshev'
             order         = str(int(fit_order)),             # Order of the fitting function, or number of splines.
             low_reject    = '0.0',                           # rejection limits below the fit in units of sigma
             high_reject   = '0.0',                           # rejection limits above the fit in units of sigma.
             niterate      = '1',                             # Number of rejection iterations.
             grow          = '0.0',                           # Reject additional points within this distance of 
                                                              # points exceeding the rejection threshold
             graphics      = 'stdgraph',                      # Graphics for IRAF to use
             cursor        = '',                              # Cursor?
             mode          = 'ql'                             # IRAF mode
            )

def check_0(filelist):
    '''
    Check that the files in the directory don't have 0 values.
    If they do have 0 values, replace them with the value of the
    nearest pixel.

    Parameters
    -------------
    filelist: glob style command with the list of files to check

    Output
    -------------
    Function will overwrite the existing files.
    '''
    import numpy as np

    # Read in the files
    files = glob.glob(filelist)

    for j in range(len(files)):

        # Extract the data
        file = fits.open(files[j])
        data = file[0].data

        # Is there 0's in the data?
        continue_function = True
        if len(np.where(data == 0.0)[0]) == 0:
            print(files[j] + " has no bad pixels, skipping" + '\n')
            continue_function = False

        if continue_function:
            # Check if there is a row of 0's somewhere that can 
            # be easily replaced
            zeros = np.where(data == 0.0)

            # If all of the zero's are along the y-axis
            if len(np.unique(zeros[0])) == 1:
                # If the zeros are along the highest y-value
                # Subtract 1 from y, otherwise add.
                if zeros[0][0] + 1 == len(data):
                    zeros_mod = tuple((zeros[0]-1, zeros[1]))
                else:
                    zeros_mod = tuple((zeros[0]+1, zeros[1]))

            # If all of the zero's are along the x-axis
            elif len(np.unique(zeros[1])) == 1:
                # If the zeros are along the highest x-value
                # Subtract 1 from x, otherwise add.
                if zeros[1][0] + 1 == len(data.T):
                    zeros_mod = tuple((zeros[0], zeros[1]-1))
                else:
                    zeros_mod = tuple((zeros[0], zeros[1]+1))

            # Otherwise just make it nan
            else:
                zeros_mod = 'nan'

            # If it worked, replace the row of 0's
            if zeros_mod != 'nan':
                data[zeros] = data[zeros_mod]

            # If that didn't work but there are still 0's remove them
            # one by one.
            while len(np.where(data == 0.0)[0]) != 0:
                Total = 0
                zeros_columns = np.where(data == 0.0)[1]
                # Check if line of the data has 0
                for i in zeros_columns:
                    how_many = len(np.where(data.T[i] == 0)[0])
                    # If there are 0's in that line
                    if how_many > 0:
                        bad  = np.where(data.T[i] == 0)
                        try:
                            good = bad[0] + 1
                            data.T[i][bad] = data.T[i][good]
                        except:
                            good = bad[0] - 1
                            data.T[i][bad] = data.T[i][good]
                        Total += how_many

                print('\n' + "*** Replaced %s bad pixels in %s ***"%(Total, files[j]) + '\n' )
                file.writeto(files[j], overwrite = True)
            else:
                print('\n' + "*** Replaced %s bad pixels in %s ***"%(len(zeros_mod[0]), files[j]) + '\n' )
                file.writeto(files[j], overwrite = True)

    del np

def iraf_background(directory, objecto, fit_sample, fit_order = 20):
    '''
    Fit and substract a line of commond background. You will 
    need to manually type the columns to sample. (i.e. 200 2000)

    Parameters
    -------------
    directory: Directory where all the science files are located.
               i.e. Objectname/Science

    Output
    -------------
    Background corrected files in the same directory as the files
    that are being correced.
    '''

    # Check that the files don't exist
    if check_existence('%s/%s*BiasFlatSky.fits'%(directory, objecto), 'iraf_background'):
        return

    # Do background subtraction
    background(input       = '@%s/AllfilesBiasFlat'%(directory),    # Input images names
               output      = '@%s/AllfilesBiasFlatSky'%(directory), # Output images names
               axis        = '2',                                   # Axis along which background is fit and substracted
               interactive = 'yes',                                 # Run the fit interactively?
               sample      =  fit_sample,                           # What points to sample in the fit
               naverage    = '1',                                   # Number of points in sample averaging
               function    = 'legendre',                            # Fitting function
               order       =  str(fit_order),                       # Order of fitting function
               low_reject  = '3.0',                                 # Low sigma rejection factor
               high_reject = '2.0',                                 # High sigma rejection factor
               niterate    = '10',                                  # Number of rejection iterations
               grow        = '0.0',                                 # Rejection growing radius
               graphics    = 'stdgraph',                            # Graphics for IRAF to use
               cursor      = '',                                    # Cursor?
               mode        = 'ql',                                  # IRAF mode
               )

def iraf_apall(directory, objecto, file_type, trace_order = 3, gain_name = 'GAIN', noise_name = 'enoise'):
    '''
    Extract a 1D spectra from a 2D image, if file_type is 'lamp' then
    the code will do everything automatically, reading in the apertures
    from the science_dir file. If the file_type is 'science' then the
    interactive mode will be active to fit for the correct aperture.

    Parameters
    -------------
    directory: Directory where all the science or lamp files are located.
               i.e. Objectname/Science
    science_dir: Directory where all the science files are located.
               i.e. Objectname/Science
    file_type: either 'science' or 'lamp' depending on what type of
               file you are reducing.

    Output
    -------------
    List of 1D spectra that have been cleaned for cosmic rays and 
    undergone optimal extraction, the files have 3 arrays:
    - Clean Spectra
    - Raw Spectra
    - Sigma
    '''

    # Check that the files don't exist
    if check_existence('%s/%s*Out.fits'%(directory, objecto), 'iraf_apall'):
        return

    # If the files are lamp files, there was no sky subtraction done
    if file_type in ['lamp', 'arc', 'HeNeAr', 'henear', 'arcs', 'Arc', 'HeArNe']:
        input_in       = '@%s/lamplistbiasflat'%directory
        output_in      = '@%s/lamplistbiasflatout'%directory
        interactive_in = 'no'

        # Open a reference science image, the first one
        with open('%s/AllfilesBiasFlatSky'%directory) as f:
            read_data = f.readline()
        f.close()
        if read_data[-1:] == '\n':
            references_in = read_data[:-1]
        else:
            references_in = read_data

    # If the files are science files
    else:
        # Definte the necessary names
        input_in = '@%s/AllfilesBiasFlatSky'%directory
        output_in = '@%s/AllfilesBiasFlatSkyOut'%directory
        interactive_in = 'yes'
        references_in = ''

    # Run Apall to extract the spectra
    apall(input         = input_in,        # List of input images
          output        = output_in,       # List of output spectra
          apertures     = '',              # Apertures
          format        = 'multispec',     # Extracted spectra format
          references    = references_in,   # List of aperture reference images
          profiles      = '',              # List of aperture profile images
          interactive   = interactive_in,  # Run task interactively?
          find          = interactive_in,  # Find apertures?
          recenter      = interactive_in,  # Recenter apertures?
          resize        = interactive_in,  # Resize apertures?
          edit          = interactive_in,  # Edit apertures? 
          trace         = interactive_in,  # Trace apertures? 
          fittrace      = interactive_in,  # Fit the traced points interactively?
          extract       = 'yes',           # Extract spectra? 
          extras        = 'yes',           # Extract sky, sigma, etc? 
          review        = 'yes',           # Review extractions?
          line          = 'INDEF',         # Dispersion line, INDEF selects the middle line
          nsum          = '500',           # Number of dispersion lines to sum or median
          # Aperture Parameters
          lower         = '-5.0',          # Lower aperture limit relative to center
          upper         = '5.0',           # Upper aperture limit relative to center
          apidtable     = '',              # Aperture ID table
          # Background Parameters
          b_function    = 'legendre',      # Background function
          b_order       = '1',             # Background function order 
          b_sample      = '-50:-25,25:50', # Background sample regions
          b_naverage    = '-3',            # Background average or median
          b_niterate    = '5',             # Background rejection iterations
          b_low_reject  = '3',             # Background lower rejection sigma
          b_high_reject = '3',             # Background upper rejection sigma
          b_grow        = '0',             # Background rejection growing radius
          # Aperture Centering Parameters
          width         = '5',             # Profile centering width
          radius        = '10',            # Profile centering radius
          threshold     = '0',             # Detection threshold for profile centering
          # Automatic Finding and ordering Parameters
          nfind         = '1',             # Number of apertures to be found automatically
          minsep        = '5',             # Minimum separation between spectra
          maxsep        = '10000.0',       # Maximum separation between spectra
          order         = 'increasing',    # Order of apertures
          # Recentering Parameters
          aprecenter    = '',              # Apertures for recentering calculation
          npeaks        = 'INDEF',         # Select brightest peaks, if INDEF all apertures will be selected.
          shift         = 'yes',           # Use average shift instead of recentering?
          # Resizing Parameters
          llimit        = 'INDEF',         # Lower aperture limit relative to center,
                                           # INDEF places the limit at the edge of the image
          ulimit        = 'INDEF',         # Upper aperture limit relative to center
          ylevel        = '0.1',           # Fraction of peak or intensity for automatic width
          peak          = 'yes',           # Is ylevel a fraction of the peak?
          bkg           = 'yes',           # Substract background in automatic width?
          r_grow        = '0.0',           # Grow limits by this factor
          avglimits     = 'no',            # Average limits over all apertures?
          # Tracing Parameters
          t_nsum        = '40',            # Number of dispersion lines to sum
          t_step        = '40',            # Tracing step
          t_nlost       = '3',             # Number of consequtive times profile is lost before quitting
          t_function    = 'legendre',      # Trace fitting function
          t_order       = str(trace_order),# Trace fitting function order
          t_sample      = '*',             # Trace sample regions
          t_naverage    = '1',             # Trace average or median
          t_niterate    = '5',             # Trace rejection iterations
          t_low_reject  = '3.0',           # Trace lower rejection sigma
          t_high_reject = '3.0',           # Trace upper rejection sigma
          t_grow        = '0.0',           # Trace rejection growing radius
          # Extraction Parameters
          background    = 'none',          # Background to subtract (none, average, median, minimum, or fit)
          skybox        = '1',             # Box car smoothing length of sky
          weights       = 'variance',      # Extraction weights (none or variance = optimal extraction)
          pfit          = 'fit1d',         # Profile fitting type (fit1d or fit2d)
          clean         = 'yes',           # Detect and replace bad pixels?
          saturation    = 'INDEF',         # Saturation level
          readnoise     = noise_name,      # Read out noise sigma
          gain          = gain_name,       # Photon gain
          lsigma        = '4.0',           # Lower region threshold
          usigma        = '4.0',           # Upper rejection threshold
          nsubaps       = '1',             # Number of subapertures per aperture
          mode          = 'ql',            # IRAF mode
          )

def individual_flats(directory, objecto, bias_file = 'bias/Bias.fits', response_sample = '10:2670', fit_order = 20, extension = ''):
    '''
    If there was only one set of flats to be used for multiple objects, then use this function.

    Parameters
    -------------
    directory: Directory where all the bias files are located.
               i.e. Objectname/Flats
    objecto: Name of the flat file to correct
    response_sample: X-axis from which to sample and fit the response
                     function
    fit_order: Order of the fit to the response function.
    '''

    # Create lists of Flats
    create_lists_flats(directory, objecto, extension)
    # Bias correct the flats
    iraf_ccdproc(directory, 'flat', bias_file = bias_file, extension = extension)
    # Create master Flat file
    iraf_flatcombine(directory,  response_sample = response_sample, fit_order = fit_order)
    # Check that the master Flat file has no 0's
    check_0(directory + '/Flat_norm.fits')

def reduce_data(directory, objecto, do_individual_flats = True, bias_file = 'bias/Bias.fits', flat_file = '', arc_name = 'arc', flat_name = 'flat', extension = '', gain_name = 'GAIN', noise_name = 'enoise'):
    '''
    Big function that encompases all of the other data reduction functions.
    Must be run in sudo, maybe it needs to run in Python2. If the file 
    structure is as follows:
    Science: object_name/Science
    Flats:   object_name/Flats
    Lamp:    object_name/Lamp
    Bias:    Bias/

    Parameters
    -------------
    object_name: Name of the folder where the targets are.

    '''
    
    # Create lists of all files
    create_lists_science(directory, objecto, extension)
    create_lists_lamp(directory, arc_name, extension)

    if do_individual_flats:
        # Create lists of Flats
        create_lists_flats(directory, flat_name, extension)
        # Bias correct the flats
        iraf_ccdproc(directory, 'flat', bias_file = bias_file, extension = extension)
        # Create master Flat file
        iraf_flatcombine(directory, response_sample = '*', fit_order = 60)
        # Check that the master Flat file has no 0's
        check_0(directory + '/Flat_norm.fits')
        # Bias and Flat correct the science and lamps
        iraf_ccdproc(directory, objecto,  bias_file = bias_file, flat_file = directory + '/Flat_norm.fits', extension = extension)
        iraf_ccdproc(directory, arc_name, bias_file = bias_file, flat_file = directory + '/Flat_norm.fits', extension = extension)
    else:
        # Bias and Flat correct the science and lamps
        iraf_ccdproc(directory, objecto,  bias_file = bias_file, flat_file = flat_file, extension = extension)
        iraf_ccdproc(directory, arc_name, bias_file = bias_file, flat_file = flat_file, extension = extension)

    # Substract the background from the images
    iraf_background(directory, objecto, fit_sample = '*', fit_order = 5)

    # Extract the Source spectra
    iraf_apall(directory, objecto, 'science', trace_order = 3, gain_name = gain_name, noise_name = noise_name) 

    # Extract the lamp spectra
    iraf_apall(directory, arc_name, 'lamp', trace_order = 3, gain_name = gain_name, noise_name = noise_name)

def example():
    '''
    Example of how to use this script
    '''
    # Create the master Bias.fits file
    iraf_zerocombine('bias')

    # Reduce data where every target has its own flat
    reduce_data('AT2018lfe'        , 'AT2018lfe'    , arc_name='HeNeAr' , flat_name = 'Qh')
    reduce_data('AT2019go'         , 'AT2019go'     , arc_name='HeNeAr' , flat_name = 'Qh')

    # Or Reduce data where you only have one set of flats by first reducing the flats
    individual_flats('flats_1.0' , 'BClamp')
    individual_flats('flats_1.25', 'BClamp')

    # And then specifying where to find them
    reduce_data('AT2018lfe', 'AT2018lfe' , do_individual_flats = False, flat_file = 'flats_1.0/Flat_norm.fits' )
    reduce_data('AT2019go' , 'AT2019go'  , do_individual_flats = False, flat_file = 'flats_1.25/Flat_norm.fits')

# Create the master Bias.fits file
iraf_zerocombine('bias', extension = 0)

# IMACS
#reduce_data('AT2018hyz', 'AT2018hyz', arc_name='HeNeAr' , flat_name = 'Qh')

# Binospec
#reduce_data('Feige_110', 'spectra', arc_name = 'HeNeAr', flat_name = 'flat', extension = 1) 

# FAST
#individual_flats('FLAT', 'FLAT')
#reduce_data('BDp284211', 'BDp284211', arc_name='COMP', flat_name = 'Qh', do_individual_flats = False, flat_file = 'FLAT/Flat_norm.fits', gain_name = 'GAIN', noise_name = 'RDNOISE')

# WHT
#iraf_zerocombine('bias_Redarm', extension = 1, datasec_key = 'RTDATSEC')
#iraf_zerocombine('bias_Bluearm', extension = 1, datasec_key = 'RTDATSEC')
#individual_flats('Flat_Redarm' , 'field', extension = 1, bias_file = 'bias_Redarm/Bias.fits')
#individual_flats('Flat_Bluearm', 'field', extension = 1, bias_file = 'bias_Redarm/Bias.fits')
#reduce_data('SP2148+286_Redarm', 'SP2148+286', arc_name='arc', flat_name = 'Qh', do_individual_flats = False, flat_file = 'Flat_Redarm/Flat_norm.fits', gain_name = 'GAIN', noise_name = 'READNOIS', bias_file = 'bias_Redarm/Bias.fits', extension = 1)

# LDSS3
#iraf_zerocombine('bias', extension = 0)
#reduce_data('AT2019itq', 'spec', arc_name = 'HeNeAr', flat_name = 'flat', gain_name = 'EGAIN', noise_name = 'ENOISE') 
#reduce_data('ltt7987'  , 'spec', arc_name = 'HeNeAr', flat_name = 'flat', gain_name = 'EGAIN', noise_name = 'ENOISE') 
