from pyraf.iraf import scombine
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
            print('%s -- %s already exists, skipping.' % (function, file_name))
        return True
    else:
        return False


def create_spec_list(directory, objecto):
    '''
    Create a series of lists of the spectra files to be combined. This assumes they are already flux calibrated.

    Parameters
    -------------
    directory: Directory where all the spectra files
    objecto: Name of object to be placed in list
    '''

    # Create lists for spectra files
    os.system('ls %s/%s*BiasFlatSkyOutWaveStd*.fits > %s/speclist' %
              (directory, objecto, directory))


def iraf_scombine(directory, objecto):

    # Run scombine to coadd spectra
    scombine(input='@%s/speclist' % (directory),
             # List of output images to be created containing the combined spectra.
             output='%s/%s_combined.fits' % (directory, objecto),
             # List of output images to be created containing the number of spectra combined.
             noutput='',
             # File name for recording log information about the combining operation.
             logfile='logfile',
             apertures='',  # List of apertures to be selected for combining.
             # (all|images|apertures) Option for grouping input spectra for combining (after selection by aperture) from one or more input images.
             group='apertures',
             # (average|median|sum) Option for combining pixels at the same dispersion coordinate.
             combine='average',
             # Type of rejection operation performed on the pixels which overlap at each dispersion coordinate.
             reject='avsigclip',
             first='yes',  # Use the first input spectrum of each set to be combined to define the dispersion coordinates for combining and output?
             w1='INDEF',
             w2='INDEF',
             dw='INDEF',
             nw='INDEF',
             log='no',
             #  The output linear or log linear wavelength scale if the dispersion of the first spectrum is not used. INDEF values are filled in from the maximum wavelength range and minimum dispersion of the spectra to be combined.
             # (none|mode|median|mean|exposure|@<file>|!<keyword>) Multiplicative image scaling to be applied.
             scale='none',
             # (none|mode|median|mean|@<file>|!<keyword>) Additive zero level image shifts to be applied.
             zero='none',
             # (none|mode|median|mean|exposure|@<file>|!<keyword>) Weights to be applied during the final averaging.
             weight='exposure',
             # Wavelength sample regions to use in computing spectrum statistics for scaling and weighting.
             sample='',

             ## Algorithm Parameters ##

             lthreshold='INDEF',
             # Low and high thresholds to be applied to the input pixels.
             hthreshold='INDEF',
             nlow='1',
             # The number of low and high pixels to be rejected by the 'minmax' algorithm.
             nhigh='1',
             # The minimum number of pixels to retain or the maximum number to reject when using the clipping algorithms (ccdclip, crreject, sigclip, avsigclip, or pclip).
             nkeep='1',
             mclip='yes',  # (ccdclip, crreject, sigclip, avsigcliip) Use the median as the estimate for the true intensity rather than the average with high and low values excluded in the 'ccdclip', 'crreject', 'sigclip', and 'avsigclip' algorithms?
             lsigma='3.',
             # (ccdclip, crreject, sigclip, avsigclip, pclip) Low and high sigma clipping factors for the 'ccdclip', 'crreject', 'sigclip', 'avsigclip', and 'pclip' algorithms.
             hsigma='3.',
             rdnoise='ENOISE',
             gain='EGAIN',
             # (ccdclip, crreject) Effective CCD readout noise in electrons, gain in electrons/DN, and sensitivity noise as a fraction.
             snoise='0.',
             # (ccdclip, crreject, sigclip, avsigclip) This parameter determines when poisson corrections are made to the computation of a sigma for images with different scale factors.
             sigscale='0.1',
             pclip='-0.5',  # (pclip) Percentile clipping algorithm parameter.
             # Number of pixels to either side of a rejected pixel to also be rejected.
             grow='0',
             # Value to use when there are no input pixels to combine for an output pixel.
             blank='0.0',
             mode='ql',  # IRAF mode
             )


def combine_spectra(directory, objecto):
    os.system('mv %s A' % directory)
    create_spec_list('A', objecto)

    if check_existence('A/%s/%s_combined.fits' % (objecto, objecto),
                       'combine_spectra'):
        return

    iraf_scombine('A', objecto)
    os.system('mv A %s' % directory)


parser = argparse.ArgumentParser()
parser.add_argument('input_dir', help='Directory of raw data', type=str)
args = parser.parse_args()

obj_list = ["AT2018ldn", "AT2019cbb", "AT2019cet", "AT2019cqg"]

for obj in obj_list:
    combine_spectra(args.input_dir + '/' + obj, obj)
