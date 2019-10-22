import numpy as np
from astropy.io import fits
import glob
from astropy.table import Table
import os

def prepare_data(file_directory = 'raw_data/*.fits', crop = False, flip = False, variables = ['DISPERSR', 'FILTER'], disperser_name = 'x270', header_extension = 0, data_extension = 0):
    '''
    Copy the raw science image into directories and rename them to something useful.
    Also crop, rotate, or flip them if specified (not yet implemented)

    Parameters
    ---------------------
    file_directory : Directory where to search for files in glob format, i.e. 'data/*.fits'
    crop           : Crop the images?
    flip           : Flip the images?
    disperser_name : Only read in images with this disperser

    Returns
    ---------------------
    Nothing, saves the .fits files to their correct directories.
    '''

    # Import file names
    Files = sorted(glob.glob(file_directory))

    # Make sure there are files
    if len(Files) == 0:
        print('No Files Found With %s'%file_directory)
        return

    # Prepare each file into the right format
    for i in range(len(Files)):
        # Open File
        print(Files[i])
        File = fits.open(Files[i], ignore_missing_end=True)

        # Figure out what the file is
        spec_vs_phot = File[1].header['MASK']
        # Is it spectra?
        if 'Longslit' in spec_vs_phot:
            image_type = 'spectra'
        elif 'imaging' in spec_vs_phot:
            image_type = 'image'
        else:
            image_type = 'other'

        # Is it a flat?
        screen = File[header_extension].header['SCRN']
        if (screen == 'deployed') and (image_type == 'spectra'):
            image_type = 'flat'

        # Is it actually a lamp?
        lamp = File[header_extension].header['HENEAR']
        if (lamp == 'on') and (image_type == 'flat'):
            image_type = 'HeNeAr'

        # Is it a bias?
        file_kind = File[header_extension].header['IMAGETYP']
        if file_kind == 'bias':
            image_type = 'bias'

        # Get file and object name
        filename    = Files[i][Files[i].find('/')+1:Files[i].find('.fits')]
        object_name = File[header_extension].header['OBJECT']

        # Get the disperser and overwrite if not consistent
        disperser = File[header_extension].header['DISPERS1']
        if (disperser != disperser_name) and (image_type != 'bias'):
            image_type = 'other'

        print(image_type, ' - ', object_name)

        # If the object is a bias frame, copy that into the bias folder
        if image_type == 'bias':
            if len(glob.glob('bias')) == 0:
                os.system('mkdir bias')
            input_image_name = '%s/%s_%s.fits'%('bias', 'bias', filename)

        # If the object is a flat, spectra, or lamp copy that into it's respective object folder name
        if image_type in ['flat', 'spectra', 'HeNeAr']:
            # Make directory with relevant information
            if len(glob.glob(object_name)) == 0:
                os.system("mkdir %s"%object_name)
            input_image_name = '%s/%s_%s.fits'%(object_name, image_type, filename)

        if image_type in ['bias', 'flat', 'spectra', 'HeNeAr']:
            # Copy the science file in the directory and rename it to something useful
            output_image_name = input_image_name[:input_image_name.find('sci')] + input_image_name[input_image_name.find('g_201')+10:]
            os.system('cp %s %s'%(Files[i], output_image_name))

            ############ Images are Copied, now Overwrite and Modify ############
            if crop:
                # Crop Limits for Binospec
                xmin, xmax =    1, 4095
                ymin, ymax = 1700, 2400

                # Crop the data
                fits_file = fits.open(output_image_name, ignore_missing_end=True)
                fits_file[data_extension].data = fits_file[data_extension].data[ymin:ymax,xmin:xmax]
                fits_file[data_extension].writeto(output_image_name, overwrite = True, output_verify = 'ignore')
                #fits.setval(output_image_name,  'BIASSEC',  value='[%s:%s,%s:%s]'%(xmin - xmin + 1, xmax - xmin, ymin - ymin + 1, ymax - ymin), ext = extension)
                fits.setval(output_image_name,  'DATASEC',  value='[%s:%s,%s:%s]'%(xmin - xmin + 1, xmax - xmin, ymin - ymin + 1, ymax - ymin), ext = header_extension)
                fits.setval(output_image_name,  'DISPAXIS', value='1', ext = header_extension)
                fits.setval(output_image_name,  'enoise'  , value='1', ext = header_extension)
                print('Cropped ' + output_image_name)

            if flip:
                # Flip the data
                fits_file = fits.open(output_image_name, ignore_missing_end=True)
                fits_file[1].data = np.flip(fits_file[1].data, axis = 1)
                fits_file.writeto(output_image_name, overwrite = True)
                print('Flipped ' + output_image_name)

def extract_fits_info(file_directory, variable_names, data_index = 0, header_index = 0, return_counts = True):
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
    Files = sorted(glob.glob(file_directory))

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
    File0 = fits.open(Files[0], ignore_missing_end=True)
    Date0 = str(File0[header_index].header['DATE-OBS'])

    # Extract data for each file
    for i in range(len(Files)):
        # Open File
        print(Files[i])
        File = fits.open(Files[i], ignore_missing_end=True)
        
        # Get File name
        filename = Files[i][Files[i].find('/')+1:Files[i].find('.fits')]

        # Get value of maximum pixel if specified
        if return_counts:
            max_counts = np.nanmax(File[data_index].data)
            max_value  = len(variables)
        else:
            max_counts = '--'
            max_value  = len(variables)

        # Create array with only the file name
        parameters = np.array([filename])

        # For each variable after the file name, extract it
        for j in range(1, max_value):
            try:
                if variables[j] == 'MAX_COUNTS':
                    parameters = np.append(parameters, max_counts)
                else:
                    parameters = np.append(parameters, File[header_index].header[variables[j]])
            # If the variable is not found, add '--'
            except:
                print('%s not found'%variables[j])
                parameters = np.append(parameters, '--')

        file_list.add_row(parameters)

    # Save output table
    file_list.write("Nightlog_%s.txt"%Date0, format='ascii.fixed_width', delimiter=None)

#### Binospec ####
#extract_fits_info('Binospec*/raw_data/*.fits', ['OBJECT', 'IMAGETYP', 'SCRN', 'EXPTIME', 'RA', 'DEC', 'DATE-OBS', 'FILTER', 'MASK', 'DISPERS1', 'DISPERS2', 'HENEAR', 'MJD', 'AIRMASS', 'EXPMODE', 'PI'], header_index = 1, data_index = 1, return_counts = False)
prepare_data(variables = [''], crop = True, flip = False, header_extension = 1, data_extension = 1)
#prepare_data(variables = [''], crop = True, flip = True, header_extension = 1, data_extension = 1)


