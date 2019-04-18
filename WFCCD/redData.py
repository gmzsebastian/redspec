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
from pyraf.iraf import imcombine
from pyraf.iraf import response
from pyraf.iraf import ccdmask
from pyraf.iraf import fixpix
from pyraf.iraf import ccdproc
from pyraf.iraf import imarith
from pyraf.iraf import background
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
    exists = subprocess.Popen('ls ' + file_name, shell=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = exists.communicate()

    # If files were returned:
    if out != b'':
        if verbose:
            print("%s -- %s already exists, skipping." % (function,
                                                          file_name))
        return True
    else:
        return False


def create_lists_science(directory, objecto):
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
    6 lists:
    - List of objects
    - List + _BiasFlat
    - List + Sky
    - List + Out
    - List + Wave
    '''

    # Create lists for science files
    os.system('ls %s/%s*.fits > %s/Allfiles' % (directory,
                                                objecto, directory))
    os.system("sed 's/.fits/_BiasFlat.fits/g' %s/Allfiles > "
              "%s/AllfilesBiasFlat" % (directory, directory))
    os.system("sed 's/.fits/Sky.fits/g' %s/AllfilesBiasFlat > "
              "%s/AllfilesBiasFlatSky" % (directory, directory))
    os.system("sed 's/.fits/Out.fits/g' %s/AllfilesBiasFlatSky > "
              "%s/AllfilesBiasFlatSkyOut" % (directory, directory))
    os.system("sed 's/.fits/Wave.fits/g' %s/AllfilesBiasFlatSkyOut > "
              "%s/AllfilesBiasFlatSkyOutWave" % (directory, directory))
    os.system("sed 's/.fits/Trim.fits/g' %s/AllfilesBiasFlatSkyOutWave > "
              "%s/AllfilesBiasFlatSkyOutWaveTrim" % (directory, directory))


def create_lists_flats(directory, objecto='flat'):
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
    os.system('ls %s/%s*.fits > %s/flatlist' % (directory, objecto, directory))
    os.system("sed 's/.fits/_Bias.fits/g' %s/flatlist > "
              "%s/flatbiaslist" % (directory, directory))


def create_lists_lamp(directory, objecto='arc'):
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
    os.system('ls %s/%s*.fits > %s/lamplist' % (directory, objecto, directory))
    os.system("sed 's/.fits/_BiasFlat.fits/g' %s/lamplist > %s/lamplistbiasflat" %
              (directory, directory))
    os.system("sed 's/.fits/Out.fits/g' %s/lamplistbiasflat > %s/lamplistbiasflatout" %
              (directory, directory))


def create_comb_list(directory, objecto):
    '''
    Create a series of lists of the input files for coadding with imcombine

    Parameters
    -------------
    directory: Directory where all the lamp files are located.
    objecto: Object name

    Output
    -------------
    1 list: Images to combine.
    '''

    # Create lists for lamp files
    os.system('ls %s/%s*.fits > %s/Comblist' % (directory, objecto, directory))


def iraf_zerocombine(directory):
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

    # Test to see if the Bias.fits file already exists
    if check_existence(directory + '/Bias.fits', 'iraf_zerocombine'):
        return

    zerocombine(input='@%s/biaslist' % directory,  # List of zero images to combine, start with @ for a list
                output='%s/Bias.fits' % directory,  # Output bias level name
                # Type of combine operation (average or median)
                combine='median',
                reject='avsigclip',  # Type of rejection algorithm
                #      none - No rejection
                #    minmax - Reject the nlow and nhigh pixels
                #   ccdclip - Reject pixels using CCD noise parameters
                #  crreject - Reject only positive pixels using CCD noise parameters
                #   sigclip - Reject pixels using a sigma clipping algorithm
                # avsigclip - Reject pixels using an averaged sigma clipping algorithm
                #     pclip - Reject pixels using sigma based on percentiles
                # Tpe of CCD image to combine, if '' all images will be used
                ccdtype='',
                process='no',  # Process images before combining?
                delete='no',  # Delete input images after combining?
                clobber='no',  # Clobber existing output image?
                scale='none',  # Scale the images? for example by exposure time
                # options: none, mode, median, mean, exposure
                statsec='',  # Image section for computing statistics
                nlow='0',  # Number of low pixels to reject
                nhigh='1',  # Number of high pixels to reject
                # If positive, minimum to keep.
                nkeep='1',
                # If negative, maximum to reject
                mclip='yes',  # Use median in sigma clipping algorithm?
                lsigma='3.0',  # Lower sigma clipping factor
                hsigma='3.0',  # Upper sigma clipping factor
                # CCD redout noise (electrons)
                rdnoise='ENOISE',
                gain='EGAIN',  # CCD gain (electrons / DN)
                snoise='0.0',  # Sensitiviy noise
                # Percentile clipping algorithm parameter.
                pclip='-0.5',
                # If positive: Specifies number of pixels above or below the median to use for computing the clipping sigma.
                # If negative: Specifies the fraction of pixels above or below the median to use.
                # The default of -0.5 selects approximately the quartile point.
                blank='0.0',  # Value if there are no pixels
                mode='ql'  # IRAF mode
                )

    print('Saved Bias.fits to %s/ \n' % directory)


def iraf_ccdproc(directory, file_type, bias_file, flat_file='', mask_file=''):
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
    if check_existence('%s/%s*_BiasFlat.fits' % (directory, file_type), 'iraf_ccdproc'):
        return

    # If correcting flat files, don't do flat corrections
    if (file_type == 'flat') or (file_type == 'flats') or (file_type == 'BClamp'):
        images_in = '@%s/flatlist' % directory
        output_in = '@%s/flatbiaslist' % directory
        flat = ''
        flatcor_in = 'no'
        mask = ''
        fixpix_in = 'no'

    # If correcting the lamp files, do flat correction
    elif file_type in ['lamp', 'arc', 'HeNeAr', 'henear', 'arcs']:
        images_in = '@%s/lamplist' % directory
        output_in = '@%s/lamplistbiasflat' % directory
        flat = flat_file
        flatcor_in = 'yes'
        mask = mask_file
        fixpix_in = 'yes'

    # If correcting sceicne files, do flat correction
    else:
        images_in = '@%s/Allfiles' % directory
        output_in = '@%s/AllfilesBiasFlat' % directory
        flat = flat_file
        flatcor_in = 'yes'
        mask = mask_file
        fixpix_in = 'yes'

    # Do the correction
    ccdproc(images=images_in,  # List of CCD images to correct
            output=output_in,  # Name of output images
            ccdtype='',  # CCD image type to correct, if empty all images will be used
            max_cache='0',  # Maximum image cachine memory in MB
            noproc='no',  # List processing steps only?
            fixpix=fixpix_in,  # Fix bad CCD lines and columns?
            overscan='no',  # Apply overscan strip correction?
            trim='yes',  # Trim the image?
            zerocor='yes',  # Apply zero level correction?
            darkcor='no',  # Apply dark count correction?
            flatcor=flatcor_in,  # Apply flat field correction?
            illumcor='no',  # Apply illumination correction?
            fringecor='no',  # Apply fringe correction?
            readcor='no',  # Convert zero level image to readout correction?
            # If yes, zero images are averaged across readout axis
            scancor='no',   # Convert flat field image to scan correction?
            # If yes, the form of scan mode correction is specified by scantype
            readaxis='line',  # Read out axis (Column or line)
            fixfile=mask,  # File describing the bad lines and columns
            biassec='',  # Overscan strip image section
            trimsec='[50:4030,1500:2600]',  # Trim data section
            zero=bias_file,  # Zero level calibration image
            dark='',  # Dark count calibration image
            flat=flat,  # Flat field images
            illum='',  # Illumination correction images
            fringe='',  # Fringe correction images
            minreplace='1.0',  # Minimum flat filed value, values below this value
            # are replaced by this value. To avoid divide by 0
            # Scan type (Shortscan or longscan), longscan is a mode
            scantype='shortscan',
            # where the CCD is clocked and read out continuously.
            nscan='1',  # Number of short scan lines
            interactive='no',  # Fit overscan interactively?
            function='legendre',  # Fitting function
            order='1',  # Number of polynomial terms or spline pieces
            sample='*',  # Sample points to fit
            naverage='1',  # Number of sample points to combine
            niterate='3',  # Number of rejection iterations
            low_reject='3.0',  # Low sigma rejection factor
            high_reject='3.0',  # High sigma rejection factor
            grow='0.0',  # Rejection growing radius
            mode='ql',  # IRAF mode
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
    if check_existence('%s/Flat_norm.fits' % directory, 'iraf_flatcombine'):
        return

    # Run flatcombine to create the average flat file
    flatcombine(input='@%s/flatbiaslist' % directory,  # List of flat field images to combine
                output='%s/Flat.fits' % directory,  # Output name for the flat field
                # Type of combine operation (average or median)
                combine='median',
                reject='avsigclip',   # Type of rejection algorithm
                #      none - No rejection
                #    minmax - Reject the nlow and nhigh pixels
                #   ccdclip - Reject pixels using CCD noise parameters
                #  crreject - Reject only positive pixels using CCD noise parameters
                #   sigclip - Reject pixels using a sigma clipping algorithm
                # avsigclip - Reject pixels using an averaged sigma clipping algorithm
                #     pclip - Reject pixels using sigma based on percentiles
                # Tpe of CCD image to combine, if '' all images will be used
                ccdtype='',
                process='no',  # Process images before combining?
                # Combine different images into groups with different subset parameters.
                subsets='no',
                delete='no',  # Delete input images after combining?
                clobber='no',  # Clobber existing output image?
                scale='none',  # Scale the images? for example by exposure time
                # options: none, mode, median, mean, exposure
                statsec='',  # Image section for computing statistics
                nlow='0',  # Number of low pixels to reject
                nhigh='1',  # Number of high pixels to reject
                # If positive, minimum to keep.
                nkeep='1',
                # If negative, maximum to reject
                mclip='yes',  # Use median in sigma clipping algorithm?
                lsigma='3.0',  # Lower sigma clipping factor
                hsigma='3.0',  # Upper sigma clipping factor
                # CCD redout noise (electrons)
                rdnoise='ENOISE',
                # CCD gain (electrons / DN)
                gain='EGAIN',
                snoise='0.0',  # Sensitiviy noise
                # Percentile clipping algorithm parameter.
                pclip='-0.5',
                # If positive: Specifies number of pixels above or below the median to use for computing the clipping sigma.
                # If negative: Specifies the fraction of pixels above or below the median to use.
                # The default of -0.5 selects approximately the quartile point.
                blank='1.0',  # Value if there are no pixels
                mode='ql'  # IRAF mode
                )

    # Determine response calibrations
    response(calibration='%s/Flat.fits' % directory,  # Images to use in determining response calibrations.
             # Images to use in determining the normalization spectrum, usually the same as calibration.
             normalization='%s/Flat.fits' % directory,
             response='%s/Flat_norm.fits' % directory,  # Output image name
             # Graph and fit interactively? yes or no
             interactive='yes',
             # Set the response to 1 when the normalization is below this threshold.
             threshold='INDEF',
             # If INDEF then don't do anything.
             # Samples of points to use in fitting the average function.
             sample='*',
             # If * all points will be used
             # Number of sample points to average or median before fitting the function.
             naverage='1',
             # Positive is average, negative is median.
             # Function to fit the average, either 'spline1', 'spline3', 'legendre', or 'chebyshev'
             function='legendre',
             # Order of the fitting function, or number of splines.
             order=str(int(fit_order)),
             # rejection limits below the fit in units of sigma
             low_reject='0.0',
             # rejection limits above the fit in units of sigma.
             high_reject='0.0',
             # Number of rejection iterations.
             niterate='1',
             grow='0.0',  # Reject additional points within this distance of
             # points exceeding the rejection threshold
             graphics='stdgraph',  # Graphics for IRAF to use
             cursor='',  # Cursor?
             mode='ql'  # IRAF mode
             )


def iraf_ccdmask(directory):

    flat_in = directory + '/Flat_norm.fits'
    mask_out = flat_in.replace('.fits', '_mask')

    if not check_existence(mask_out, 'iraf_ccdmask'):
        ccdmask(image=flat_in,  # Input image. Ideally a flat field
                mask=mask_out,  # Output pixel mask
                ncmed=10,  # Column size of moving median rectangle (min 3)
                nlmed=10,  # Line size of moving median rectangle (min 3)
                ncsig=20,  # Column size of search box (Box min 100 pixels)
                nlsig=20,  # Line size of search box (Box min 100 pixels)
                lsigma=6,  # Sigma to select below average
                hsigma=6,  # Sigma to select above average
                ngood=7,  # Gap size of pixels to flag as bad
                linterp=2,  # Mask code for interpolation along lines
                cinterp=3,  # Mask code for interpolation along columns
                eqinterp=2,  # Mask code for interpolation along both.
                )


def iraf_imcombine(directory, objecto):

    imcombine(  # List of input images to combine.
        input='@%s/Comblist' % directory,
        # Output combined image(s).
        output='%s' % objecto,
        # Optional output multiextension FITS file(s).
        headers="",
        # Optional output bad pixel mask(s) with good values of 0 and bad values of 1.
        bpmasks="",
        # Optional output mask file(s) identifying rejected or excluded pixels.
        rejmask="",
        # Optional output pixel mask(s) giving the number of input pixels rejected or excluded from the input images.
        nrejmasks="",
        # Optional output exposure mask(s) giving the sum of the exposure values of the input images with non - zero weights that contributed to that pixel.
        expmasks="",
        # Optional output sigma image(s).
        sigma='',
        logfile="logfile",  # Optional output log file.
        # (average | median | sum) Type of combining operation performed on the final set of pixels(after offsetting, masking, thresholding, and rejection).
        combine="average",
        # (none | minmax | ccdclip | crreject | sigclip | avsigclip | pclip) Type of rejection operation performed on the pixels remaining after offsetting, masking and thresholding.
        reject="avsigclip",
        # Project(combine) across the highest dimension of the input images?
        project="no",
        # (none | short | ushort | integer | long | real | double) Output image pixel datatype.
        outtype="real",
        # Output region limits specified as pairs of whitespace separated values.
        outlimits="",
        offsets="none",  # (none|wcs|world|physical|grid|<filename>)
        # (none|goodvalue|badvalue|goodbits|badbits|!<keyword>) Type of pixel masking to use.
        masktype="none",
        maskvalue="0",  # Mask value used with the masktype parameter.
        # Output value to be used when there are no pixels.
        blank='0.',
        # (none|mode|median|mean|exposure|@ < file > |!< keyword > ) Multiplicative image scaling to be applied.
        scale="exposure",
        # (none|mode|median|mean|@ < file > |!< keyword > ) Additive zero level image shifts to be applied.
        zero="none",
        # (none|mode|median|mean|exposure|@ < file > |!< keyword > ) Weights to be applied during the final averaging.
        weight="exposure",
        # Section of images to use in computing image statistics for scaling and weighting.
        statsec="[50:4030,1500:2600]",
        # Image header keyword to be used with the exposure scaling and weighting options.
        expname="EXPTIME",

        ## Algorithm Parameters ##

        lthreshold='INDEF',
        # Low and high thresholds to be applied to the input pixels.
        hthreshold='INDEF',
        nlow='1',
        # (minmax) The number of low and high pixels to be rejected by the "minmax" algorithm.
        nhigh='1',
        # The minimum number of pixels to retain or the maximum number to reject when using the clipping algorithms(ccdclip, crreject, sigclip, avsigclip, or pclip).
        nkeep='1',
        mclip='yes',  # (ccdclip, crreject, sigclip, avsigcliip) Use the median as the estimate for the true intensity rather than the average with high and low values excluded in the "ccdclip", "crreject", "sigclip", and "avsigclip" algorithms?
        lsigma='3.',
        # (ccdclip, crreject, sigclip, avsigclip, pclip) Low and high sigma clipping factors for the "ccdclip", "crreject", "sigclip", "avsigclip", and "pclip" algorithms.
        hsigma='3.',
        rdnoise="ENOISE",
        # (ccdclip, crreject, sigclip, avsigclip) This parameter determines when poisson corrections are made to the computation of a sigma for images with different scale factors.
        sigscale='0.1',
        # (pclip) Percentile clipping algorithm parameter.
        pclip='-0.5',
        # Radius in pixels for additional pixel to be rejected in an image with a rejected pixel from one of the rejection algorithms.
        grow='0.',
        mode='ql',  # IRAF mode
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
                    zeros_mod = tuple((zeros[0] - 1, zeros[1]))
                else:
                    zeros_mod = tuple((zeros[0] + 1, zeros[1]))

            # If all of the zero's are along the x-axis
            elif len(np.unique(zeros[1])) == 1:
                # If the zeros are along the highest x-value
                # Subtract 1 from x, otherwise add.
                if zeros[1][0] + 1 == len(data.T):
                    zeros_mod = tuple((zeros[0], zeros[1] - 1))
                else:
                    zeros_mod = tuple((zeros[0], zeros[1] + 1))

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
                        bad = np.where(data.T[i] == 0)
                        try:
                            good = bad[0] + 1
                            data.T[i][bad] = data.T[i][good]
                        except:
                            good = bad[0] - 1
                            data.T[i][bad] = data.T[i][good]
                        Total += how_many

                print('\n' + "*** Replaced %s bad pixels in %s ***" %
                      (Total, files[j]) + '\n')
                file.writeto(files[j], overwrite=True)
            else:
                print('\n' + "*** Replaced %s bad pixels in %s ***" %
                      (len(zeros_mod[0]), files[j]) + '\n')
                file.writeto(files[j], overwrite=True)

    del np


def iraf_background(directory, objecto, fit_sample, fit_order=20):
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
    if check_existence('%s/%s*BiasFlatSky.fits' % (directory, objecto), 'iraf_background'):
        return

    # Do background subtraction
    background(input='@%s/AllfilesBiasFlat' % (directory),  # Input images names
               # Output images names
               output='@%s/AllfilesBiasFlatSky' % (directory),
               # Axis along which background is fit and substracted
               axis='2',
               interactive='yes',  # Run the fit interactively?
               sample=fit_sample,  # What points to sample in the fit
               # Number of points in sample averaging
               naverage='1',
               function='legendre',  # Fitting function
               # Order of fitting function
               order=str(fit_order),
               low_reject='3.0',  # Low sigma rejection factor
               high_reject='3.0',  # High sigma rejection factor
               niterate='10',  # Number of rejection iterations
               grow='0.0',  # Rejection growing radius
               graphics='stdgraph',  # Graphics for IRAF to use
               cursor='',  # Cursor?
               mode='ql',  # IRAF mode
               )


def iraf_apall(directory, objecto, file_type, trace_order=3):
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
    if check_existence('%s/%s*Out.fits' % (directory, objecto), 'iraf_apall'):
        return

    # If the files are lamp files, there was no sky subtraction done
    if file_type in ['lamp', 'arc', 'HeNeAr', 'henear', 'arcs']:
        input_in = '@%s/lamplistbiasflat' % directory
        output_in = '@%s/lamplistbiasflatout' % directory
        interactive_in = 'no'

        # Open a reference science image, the first one
        with open('%s/AllfilesBiasFlatSky' % directory) as f:
            read_data = f.readline()
        f.close()
        if read_data[-1:] == '\n':
            references_in = read_data[:-1]
        else:
            references_in = read_data

    # If the files are science files
    else:
        # Definte the necessary names
        input_in = '@%s/AllfilesBiasFlatSky' % directory
        output_in = '@%s/AllfilesBiasFlatSkyOut' % directory
        interactive_in = 'yes'
        references_in = ''

    # Run Apall to extract the spectra
    apall(input=input_in,  # List of input images
          output=output_in,  # List of output spectra
          apertures='',  # Apertures
          format='multispec',  # Extracted spectra format
          references=references_in,  # List of aperture reference images
          profiles='',  # List of aperture profile images
          interactive=interactive_in,  # Run task interactively?
          find=interactive_in,  # Find apertures?
          recenter=interactive_in,  # Recenter apertures?
          resize=interactive_in,  # Resize apertures?
          edit=interactive_in,  # Edit apertures?
          trace=interactive_in,  # Trace apertures?
          fittrace=interactive_in,  # Fit the traced points interactively?
          extract='yes',  # Extract spectra?
          extras='yes',  # Extract sky, sigma, etc?
          review='yes',  # Review extractions?
          line='INDEF',  # Dispersion line, INDEF selects the middle line
          nsum='500',  # Number of dispersion lines to sum or median
          # Aperture Parameters
          lower='-5.0',  # Lower aperture limit relative to center
          upper='5.0',  # Upper aperture limit relative to center
          apidtable='',  # Aperture ID table
          # Background Parameters
          b_function='legendre',  # Background function
          b_order='1',  # Background function order
          b_sample='-50:-25,25:50',  # Background sample regions
          b_naverage='-3',  # Background average or median
          b_niterate='5',  # Background rejection iterations
          b_low_reject='3',  # Background lower rejection sigma
          b_high_reject='3',  # Background upper rejection sigma
          b_grow='0',  # Background rejection growing radius
          # Aperture Centering Parameters
          width='5',  # Profile centering width
          radius='10',  # Profile centering radius
          threshold='0',  # Detection threshold for profile centering
          # Automatic Finding and ordering Parameters
          nfind='1',  # Number of apertures to be found automatically
          minsep='5',  # Minimum separation between spectra
          maxsep='10000.0',  # Maximum separation between spectra
          order='increasing',  # Order of apertures
          # Recentering Parameters
          aprecenter='',  # Apertures for recentering calculation
          # Select brightest peaks, if INDEF all apertures will be selected.
          npeaks='INDEF',
          shift='yes',  # Use average shift instead of recentering?
          # Resizing Parameters
          llimit='INDEF',  # Lower aperture limit relative to center,
          # INDEF places the limit at the edge of the image
          ulimit='INDEF',  # Upper aperture limit relative to center
          ylevel='0.1',  # Fraction of peak or intensity for automatic width
          peak='yes',  # Is ylevel a fraction of the peak?
          bkg='yes',  # Subtract background in automatic width?
          r_grow='0.0',  # Grow limits by this factor
          avglimits='no',  # Average limits over all apertures?
          # Tracing Parameters
          t_nsum='40',  # Number of dispersion lines to sum
          t_step='40',  # Tracing step
          t_nlost='3',  # Number of consequtive times profile is lost before quitting
          t_function='legendre',  # Trace fitting function
          t_order=str(trace_order),  # Trace fitting function order
          t_sample='*',  # Trace sample regions
          t_naverage='1',  # Trace average or median
          t_niterate='5',  # Trace rejection iterations
          t_low_reject='3.0',  # Trace lower rejection sigma
          t_high_reject='3.0',  # Trace upper rejection sigma
          t_grow='0.0',  # Trace rejection growing radius
          # Extraction Parameters
          # Background to subtract (none, average, median, minimum, or fit)
          background='median',
          skybox='1',  # Box car smoothing length of sky
          # Extraction weights (none or variance = optimal extraction)
          weights='variance',
          pfit='fit1d',  # Profile fitting type (fit1d or fit2d)
          clean='yes',  # Detect and replace bad pixels?
          saturation='INDEF',  # Saturation level
          readnoise='ENOISE',  # Read out noise sigma
          gain='EGAIN',  # Photon gain
          lsigma='4.0',  # Lower region threshold
          usigma='4.0',  # Upper rejection threshold
          nsubaps='1',  # Number of subapertures per aperture
          mode='ql',  # IRAF mode
          )


def individual_flats(directory, objecto, bias_file,
                     response_sample='*', fit_order=20):
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

    # Test to see if the Bias.fits file already exists
    if check_existence(directory + '/Flat_norm.fits', 'individual_flats'):
        return

    # Create lists of Flats
    create_lists_flats(directory, objecto=objecto)
    # Bias correct the flats
    iraf_ccdproc(directory, 'flat', bias_file=bias_file)
    # Create master Flat file
    iraf_flatcombine(
        directory,  response_sample=response_sample, fit_order=fit_order)
    # Check that the master Flat file has no 0's
    check_0(directory + '/Flat_norm.fits')
    iraf_ccdmask(directory)


def combine_images(directory, objecto):

    if not os.path.exists(directory + '/single_images'):
        os.system('mkdir %s/single_images' % directory)

    create_comb_list(directory, objecto)

    output_image = "%s/%s_Comb.fits" % (directory, objecto)

    if check_existence(output_image, "combine_images"):
        return

    iraf_imcombine(directory, output_image)

    file_list = glob.glob("%s/%s_*.fits" % (directory, objecto))

    for f in file_list:
        if (not "_Comb" in f):
            os.system('mv %s %s/single_images' % (f, directory))


def reduce_data(directory, objecto,
                do_individual_flats=True,
                science_flat='', arc_flat='',
                science_mask='', arc_mask='',
                bias_file='',
                arc_name='arc', flat_name='flat'):
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
    create_lists_science(directory, objecto)
    create_lists_lamp(directory, objecto=arc_name)

    if do_individual_flats:
        # Create lists of Flats
        create_lists_flats(directory, objecto=flat_name)
        # Bias correct the flats
        iraf_ccdproc(directory, 'flat', bias_file=bias_file)
        # Create master Flat file
        iraf_flatcombine(directory, response_sample='*', fit_order=40)
        # Check that the master Flat file has no 0's
        check_0(directory + '/Flat_norm.fits')
        # Create Bad Pixel Mask
        iraf_ccdmask(directory)
        # Bias and Flat correct the science and lamps
        iraf_ccdproc(directory, objecto,  bias_file=bias_file,
                     flat_file=directory + '/Flat_norm.fits',
                     mask_file=directory + '/Flat_norm_mask.pl')
        iraf_ccdproc(directory, arc_name, bias_file=bias_file,
                     flat_file=directory + '/Flat_norm.fits',
                     mask_file=directory + '/Flat_norm_mask.pl')
    else:
        # Bias and Flat correct the science and lamps
        iraf_ccdproc(directory, objecto,
                     bias_file=bias_file,
                     flat_file=science_flat,
                     mask_file=science_mask)
        iraf_ccdproc(directory, arc_name,
                     bias_file=bias_file,
                     flat_file=arc_flat,
                     mask_file=arc_mask)

    # Substract the background from the images
    iraf_background(directory, objecto, fit_sample='*', fit_order=10)

    # Extract the Source spectra
    iraf_apall(directory, objecto, 'science', trace_order=3)

    # Extract the lamp spectra
    iraf_apall(directory, arc_name, 'lamp', trace_order=3)


parser = argparse.ArgumentParser()
parser.add_argument("input_dir", help="Directory of raw data", type=str)
args = parser.parse_args()

# Create the master Bias.fits file
iraf_zerocombine(args.input_dir + '/bias')

individual_flats(args.input_dir + '/flat_750', 'flat',
                 bias_file=args.input_dir + '/bias/Bias.fits')
individual_flats(args.input_dir + '/flat_150', 'flat',
                 bias_file=args.input_dir + '/bias/Bias.fits')

# Determine bad pixel mask
# Process Science Targets
obj_list = ["AT2018iao", "AT2019ahk", "AT2019aov", "AT2019aqv"]

for obj in obj_list:
    combine_images(args.input_dir + '/' + obj, obj)
    reduce_data(args.input_dir + '/' + obj, obj,
                do_individual_flats=False,
                science_flat=args.input_dir + '/flat_150/Flat_norm.fits',
                arc_flat=args.input_dir + '/flat_150/Flat_norm.fits',
                science_mask=args.input_dir + '/flat_150/Flat_norm_mask.pl',
                arc_mask=args.input_dir + '/flat_150/Flat_norm_mask.pl',
                bias_file=args.input_dir + '/bias/Bias.fits',
                arc_name='HeNeAr')

# Process Standards as they have different arc and science flats.
std_list = ["LTT3218"]
for obj in std_list:
    combine_images(args.input_dir + '/' + obj, obj)
    reduce_data(args.input_dir + '/' + obj, obj,
                do_individual_flats=False,
                science_flat=args.input_dir + '/flat_750/Flat_norm.fits',
                arc_flat=args.input_dir + '/flat_150/Flat_norm.fits',
                science_mask=args.input_dir + '/flat_750/Flat_norm_mask.pl',
                arc_mask=args.input_dir + '/flat_150/Flat_norm_mask.pl',
                bias_file=args.input_dir + '/bias/Bias.fits',
                arc_name='HeNeAr')
