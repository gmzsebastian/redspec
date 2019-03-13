import numpy as np
from astropy.io import fits
import glob
from astropy.table import Table
import os

def prepare_data(file_directory = 'raw_data/*.fits', crop = False, rotate = False, flip = False):
    '''
    Copy the raw science image into directories and rename them to something useful.
    Also crop, rotate, or flip them if specified (not yet implemented)

    Parameters
    ---------------------
    file_directory: Directory where to search for files in glob format, i.e. 'data/*.fits'
    crop          : Crop the images?
    rotate        : Rotate the images?
    flip          : Flip the images?

    Returns
    ---------------------
    Nothing, saves the .fits files to their correct directories.
    '''

    # Import file names
    Files = glob.glob(file_directory)

    # Make sure there are files
    if len(Files) == 0:
        print('No Files Found With %s'%file_directory)
        return

    # Prepare each file into the right format
    for i in range(len(Files)):
        # Open File
        print(Files[i])
        File = fits.open(Files[i])

        # Get File name
        filename = Files[i][Files[i].find('/')+1:Files[i].find('.fits')]

        # Get the type of file (arc, flat, object, etc.)
        file_type = File[0].header['IMAGETYP']

        # If it is not an object, copy into the current directory
        if file_type == 'zero':
            if len(glob.glob('bias')) == 0:
                os.system('mkdir bias')
            os.system('cp %s %s/%s_%s.fits'%(Files[i], 'bias', File[0].header['OBJECT'], filename))
        elif file_type != 'object':
            os.system('cp %s %s_%s.fits'%(Files[i], File[0].header['OBJECT'], filename))
        else:
            # Print relevant information
            objecto   = File[0].header['OBJECT']
            slit      = File[0].header['APERTURE']
            disperser = File[0].header['DISPERSE']
            cenwave   = File[0].header['CENWAVE']
            print('Object \t Slit \t\t Disperser \t Central Wave')
            print('%s \t %s \t %s \t %s'%(objecto, slit, disperser, cenwave))
            print('\n')

            directory_name = '%s_%s_%s_%s'%(objecto, slit, disperser, cenwave)
            # Make directory with relevant information
            if len(glob.glob(directory_name)) == 0:
                os.system("mkdir %s"%directory_name)

            # Copy the science file in the directory and rename it to something useful
            os.system('cp %s %s/%s_%s.fits'%(Files[i], directory_name, File[0].header['OBJECT'], filename))

def extract_fits_info(file_directory, variable_names, data_index = 0, return_counts = True):
    '''
    Extract relevant information from the header of several fits files.

    Parameters
    ---------------------
    variable_names: Names of the header items you wish to extract, i.e. ['OBJECT', 'EXPTIME', 'RA']
    file_directory: Directory where to search for files in glob format, i.e. 'data/*.fits'
    data_index    : In what index of a multi-dimensional fits file is the data stored, i.e. 0
    return_counts : Include a column with the value of the hottest pixel in the image? 

    Returns
    ---------------------
    Nothing, saves a file to the current directory in an Astropy Table format.
    '''

    # Import file names
    Files = glob.glob(file_directory)

    # Make sure there are files
    if len(Files) == 0:
        print('No Files Found With %s'%file_directory)
        return

    # Variables to extract
    variables = np.array(variable_names)

    # Add the filename and value of maximum pixel if specified
    variables = np.append('FILENAME', variables)
    if return_counts:
        variables = np.append(variables, 'MAX_COUNTS')

    # Create table to put that data in
    file_list = Table(names = variables, dtype=['S']*len(variables))

    # Get the date of the first file for the file name
    File0 = fits.open(Files[0])
    Date0 = str(File0[0].header['DATE-OBS'])

    # Extract data for each file
    for i in range(len(Files)):
        # Open File
        print(Files[i])
        File = fits.open(Files[i])
        
        # Get File name
        filename = Files[i][Files[i].find('/')+1:Files[i].find('.fits')]

        # Get value of maximum pixel if specified
        if return_counts:
            max_counts = np.max(File[data_index].data)
            max_value  = len(variables)
        else:
            max_counts = '--'
            max_value  = len(variables) - 1

        # Create array with only the file name
        parameters = np.array([filename])

        # For each variable after the file name, extract it
        for j in range(1, max_value):
            try:
                if variables[j] == 'MAX_COUNTS':
                    parameters = np.append(parameters, max_counts)
                else:
                    parameters = np.append(parameters, File[0].header[variables[j]])
            # If the variable is not found, add '--'
            except:
                print('%s not found'%variables[j])
                parameters = np.append(parameters, '--')

        file_list.add_row(parameters)

    # Save output table
    file_list.write("Nightlog_%s.txt"%Date0, format='ascii.tab')

extract_fits_info('raw_data/*.fits', ['OBJECT', 'EXPTIME', 'RA', 'DEC', 'DATE-OBS', 'UT', 'AIRMASS', 'APERTURE', 'DISPERSE', 'CENWAVE'])

prepare_data()
