import numpy as np
from astropy.io import fits
import glob
from astropy.table import Table
import os
import argparse


def prepare_data(file_directory,
                 crop=False, rotate=False, flip=False,
                 variables=['APERTURE', 'FILTER']):
    '''
    Copy the raw science image into directories and
    rename them to something useful.
    Also crop, rotate, or flip them if specified (not yet implemented)

    Parameters
    ---------------------
    file_directory: Directory where to search for files
    in glob format, i.e. 'data/*.fits'
    crop          : Crop the images?
    rotate        : Rotate the images?
    flip          : Flip the images?

    Returns
    ---------------------
    Nothing, saves the .fits files to their correct directories.
    '''

    # Import file names
    Files = sorted(glob.glob(file_directory + '/raw_data/*.fits'))

    # Make sure there are files
    if len(Files) == 0:
        print('No Files Found With %s' % file_directory)
        return

    # Prepare each file into the right format
    for i in range(len(Files)):
        # Open File
        print(Files[i])
        File = fits.open(Files[i])
        file_kind = File[0].header['FILTER']
        file_type = File[0].header['EXPTYPE']

        if (file_kind in ['grism_blue']) or (file_type in ['bias', 'Bias']):
            # Get File name
            filename = Files[i][Files[i].find(
                'ccd') + 3:Files[i].find('.fits')]
            full_name = File[0].header['OBJECT']

            try:
                name_break = full_name.find(' ')
            except AttributeError:
                full_name = File[0].header['EXPTYPE']
                name_break = full_name.find(' ')

            if name_break != -1:
                object_name = full_name[:name_break]
                type_name = full_name[1 + name_break:]
            elif full_name == 'arc':
                # Get Object name from previous exposure
                with fits.open(Files[i - 1]) as tmp:
                    tmp_name = tmp[0].header['OBJECT']
                object_name = tmp_name
                type_name = 'HeNeAr'
            elif full_name == 'dflat':
                object_name = 'flat'
                type_name = 'flat'
            else:
                object_name = full_name
                type_name = full_name

            # If the object is a bias frame, copy that inot the bias folder
            if file_type in ['zero', 'Bias', 'Zero', 'bias']:
                if len(glob.glob('bias')) == 0:
                    os.system('mkdir %s/bias' % file_directory)
                os.system('cp %s %s/%s_%s.fits' % (Files[i],
                                                   file_directory + '/bias',
                                                   object_name, filename))
            elif file_type in ['Flat', 'flat', 'object', 'Object']:
                print_name = object_name
                directory_name = file_directory + '/' + object_name

                if variables[0] != '':
                    for k in range(len(variables)):
                        value = File[0].header[variables[k]].replace(' ', '')
                        print_name += ' \t ' + value
                        directory_name += '_' + value

                print(print_name)
                print('\n')

                # Make directory with relevant information
                if len(glob.glob(directory_name)) == 0:
                    os.system("mkdir %s" % directory_name)

                # Copy the science file in the directory
                # and rename it to something useful
                os.system('cp %s %s/%s_%s.fits' % (Files[i],
                                                   directory_name,
                                                   type_name, filename))

            # Once saved, rotate or flip if specified
            if rotate:
                # Rotate the Data by 90 degrees. Top will be left
                if file_type in ['zero', 'Bias', 'Zero', 'bias']:
                    file_name = '%s/%s_%s.fits' % (file_directory + '/bias',
                                                   type_name, filename)
                else:
                    file_name = '%s/%s_%s.fits' % (directory_name,
                                                   type_name, filename)
                fits_file = fits.open(file_name)
                fits_file[0].data = np.rot90(fits_file[0].data)  # ,k = 1)

                # Modify the bias and data sec variables
                biassec = fits_file[0].header['BIASSEC']
                bias_coma = biassec.find(',')
                one_bias = biassec[1:bias_coma]
                two_bias = biassec[bias_coma + 1:-1]
                out_bias = '[%s,%s]' % (two_bias, one_bias)

                datasec = fits_file[0].header['DATASEC']
                data_coma = datasec.find(',')
                one_data = datasec[1:data_coma]
                two_data = datasec[data_coma + 1:-1]
                out_data = '[%s,%s]' % (two_data, one_data)

                fits_file.writeto(file_name, overwrite=True)
                fits.setval(file_name,  'BIASSEC', value=out_bias)
                fits.setval(file_name,  'DATASEC', value=out_data)
                fits.setval(file_name, 'DISPAXIS', value='1')

                print('Rotated ' + file_name)

            if crop:
                # Rotate the Data by 270 degrees. Top will be left
                if file_type in ['zero', 'Bias', 'Zero', 'bias']:
                    file_name = '%s/%s_%s.fits' % (file_directory + '/bias',
                                                   type_name, filename)
                else:
                    file_name = '%s/%s_%s.fits' % (directory_name,
                                                   type_name, filename)

                # Crop the data
                fits_file = fits.open(file_name)
                fits_file[0].data = fits_file[0].data[343:617, 149:4096]
                fits_file.writeto(file_name, overwrite=True)
                print('Cropped ' + file_name)


def extract_fits_info(file_directory,
                      variable_names, data_index=0,
                      return_counts=True):
    '''
    Extract relevant information from the header of several fits files.

    Parameters
    ---------------------
    variable_names: Names of the header items
    you wish to extract, i.e. ['OBJECT', 'EXPTIME', 'RA']
    file_directory: Directory where to search
    for files in glob format, i.e. 'data/*.fits'
    data_index    : In what index of a multi-dimensional
    fits file is the data stored, i.e. 0
    return_counts : Include a column with the value of
    the hottest pixel in the image?

    Returns
    ---------------------
    Nothing, saves a file to the current directory in an Astropy Table format.
    '''

    # Import file names
    Files = sorted(glob.glob(file_directory + '/raw_data/*.fits'))

    # Make sure there are files
    if len(Files) == 0:
        print('No Files Found With %s' % file_directory)
        return

    # Variables to extract
    variables = np.array(variable_names)

    # Add the filename and value of maximum pixel if specified
    variables = np.append('FILENAME', variables)
    if return_counts:
        variables = np.append(variables, 'MAX_COUNTS')

    # Create table to put that data in
    file_list = Table(names=variables, dtype=['S'] * len(variables))

    # Get the date of the first file for the file name
    File0 = fits.open(Files[0])
    Date0 = str(File0[0].header['DATE-OBS'])

    # Extract data for each file
    for i in range(len(Files)):
        # Open File
        print(Files[i])
        File = fits.open(Files[i])

        # Get File name
        filename = Files[i][Files[i].find('ccd') + 3:Files[i].find('.fits')]

        # Get value of maximum pixel if specified
        if return_counts:
            max_counts = np.max(File[data_index].data)
            max_value = len(variables)
        else:
            max_counts = '--'
            max_value = len(variables) - 1

        # Create array with only the file name
        parameters = np.array([filename])

        # For each variable after the file name, extract it
        for j in range(1, max_value):
            try:
                if variables[j] == 'MAX_COUNTS':
                    parameters = np.append(parameters, max_counts)
                else:
                    parameters = np.append(parameters,
                                           File[0].header[variables[j]])
            # If the variable is not found, add '--'
            except KeyError:
                print('%s not found' % variables[j])
                parameters = np.append(parameters, '--')

        file_list.add_row(parameters)

    # Save output table
    output_str = file_directory + "/Nightlog_%s.txt" % Date0
    file_list.write(output_str, format='ascii.tab', overwrite=True)


parser = argparse.ArgumentParser()
parser.add_argument("input_dir", help="Directory of raw data", type=str)
args = parser.parse_args()

extract_fits_info(args.input_dir,
                  ['OBJECT', 'EXPTYPE', 'EXPTIME',
                   'RA', 'DEC', 'DATE-OBS', 'UT-TIME',
                   'APERTURE', 'FILTER'])

prepare_data(args.input_dir, variables=[''])