import numpy as np
from astropy.io import fits
import glob
from astropy.table import Table
import os

def prepare_data(file_directory = 'raw_data/*.fits', instrument = '', crop = False, rotate = False, flip = False, variables = ['DISPERSR', 'FILTER'], break_character = '-', filter_name = 'Spectroscopic2', filterkey = 'FILTER', disperser = 'DISPERSR', data_index = 0, datasec_key = 'DATASEC', biassec_key = 'BIASSEC', rotations = 3, header_index = 0, objname = 'OBJECT'):
    '''
    Copy the raw science image into directories and rename them to something useful.
    Also crop, rotate, or flip them if specified (not yet implemented)

    Parameters
    ---------------------
    file_directory : Directory where to search for files in glob format, i.e. 'data/*.fits'
    instrument     : Instrument being used, right now it only affects the crop factor.
    crop           : Crop the images?
    rotate         : Rotate the images?
    flip           : Flip the images?
    filter_name    : IMACS = Spectroscopic2
                     Something = Bessell_R2
    disperser      : Name of disperser variable
    data_index     : Where is the data?
    rotations      : 3 means top will be left

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
        file_kind = File[header_index].header[disperser]
        if filter_name != '':
            filter_kind = File[header_index].header[filterkey]
        else:
            filter_kind = ''
        full_name = File[header_index].header[objname]

        is_it_calibration = np.array([k in full_name for k in ['BIAS', 'bias', 'Bias', 'Flat', 'FLAT', 'flat', 'ZERO']])

        if ('Gri' in file_kind) or (filter_name in filter_kind) or is_it_calibration.any():
            # Get File name
            filename    = Files[i][Files[i].find('/')+1:Files[i].find('.fits')]
            full_name   = File[header_index].header[objname]

            # Break the name in two if there are two words
            name_break  = full_name.find(break_character)
            if name_break not in [-1, 0]:
                object_name = full_name[:name_break]
                type_name   = full_name[1+name_break:].replace(' ', '')
            else:
                object_name = full_name
                type_name   = full_name.replace(' ', '')

            # Get the type of file (arc, flat, object, etc.)
            try:
                try:
                    file_type = File[header_index].header['IMAGETYP']
                except:
                    file_type = File[header_index].header['EXPTYPE']
            except:
                file_type = File[header_index].header['OBSTYPE']

            # If the object is a bias frame, copy that inot the bias folder
            if file_type in ['zero', 'Bias', 'BIAS', 'Zero', 'bias', 'ZERO']:
                print_name     = 'bias'
                directory_name = 'bias'
                if variables[0] != '':
                    for k in range(len(variables)):
                        value            = File[header_index].header[variables[k]].replace(' ', '')
                        print_name      += ' \t ' + value
                        directory_name  += '_' + value

                if len(glob.glob(directory_name)) == 0:
                    os.system('mkdir ' + directory_name)

                os.system('cp %s %s/%s_%s.fits'%(Files[i], directory_name, object_name, filename))

            elif file_type.replace(' ', '') in ['object', 'Object', 'SPECTRUM', 'COMP', 'ARC', 'OBJECT', 'OBJNAME', 'arc', 'LAMPFLAT', 'FLAT', 'Flat', 'flat', 'comp', 'Comp']:
                print_name      = object_name
                directory_name  = object_name

                if variables[0] != '':
                    for k in range(len(variables)):
                        value            = File[header_index].header[variables[k]].replace(' ', '')
                        print_name      += ' \t ' + value
                        directory_name  += '_' + value

                print(print_name)
                print('\n')

                # Make directory with relevant information
                if len(glob.glob(directory_name)) == 0:
                    os.system("mkdir %s"%directory_name)

                # Copy the science file in the directory and rename it to something useful
                os.system('cp %s %s/%s_%s.fits'%(Files[i], directory_name, type_name, filename))

            if file_type.replace(' ', '') in ['zero', 'Bias', 'BIAS', 'Zero', 'bias', 'ZERO', 'object', 'Object', 'SPECTRUM', 'COMP', 'ARC', 'OBJECT', 'OBJNAME', 'arc', 'FLAT', 'LAMPFLAT', 'Flat', 'flat', 'comp', 'Comp']:
                # Once saved, rotate or flip if specified
                if rotate:
                    # Rotate the Data by 90 degrees. Top will be left
                    file_name = '%s/%s_%s.fits'%(directory_name, type_name, filename)
                    fits_file = fits.open(file_name, ignore_missing_end=True)
                    fits_file[data_index].data = np.rot90(fits_file[data_index].data, k = rotations)

                    # Modify the bias and data sec variables
                    biassec = fits_file[data_index].header[biassec_key]
                    bias_coma = biassec.find(',')
                    one_bias = biassec[1:bias_coma]
                    two_bias = biassec[bias_coma+1:-1]
                    out_bias = '[%s,%s]'%(two_bias, one_bias)

                    datasec = fits_file[data_index].header[datasec_key]
                    data_coma = datasec.find(',')
                    one_data = datasec[1:data_coma]
                    two_data = datasec[data_coma+1:-1]
                    out_data = '[%s,%s]'%(two_data, one_data)

                    fits_file.writeto(file_name, overwrite = True)
                    fits.setval(file_name,   biassec_key,  value=out_bias, ext = header_index)
                    fits.setval(file_name,   datasec_key,  value=out_data, ext = header_index)
                    fits.setval(file_name,  'DISPAXIS', value='1', ext = header_index)

                    print('Rotated ' + file_name)

                if crop:
                    # Rotate the Data by 270 degrees. Top will be left
                    if file_type in ['zero', 'Bias', 'BIAS', 'Zero', 'bias', 'ZERO']:
                        crop_name = '%s/%s_%s.fits'%(directory_name, type_name, filename)
                    else:
                        crop_name = '%s/%s_%s.fits'%(directory_name, type_name, filename)

                    # Crop the spectrum image
                    if instrument == '':
                        xmin, xmax = 0, fits_file[data_index].data.shape[0]
                        ymin, ymax = 0, fits_file[data_index].data.shape[1]
                    elif instrument == 'IMACS1':
                        xmin, xmax = 60, 4105
                        ymin, ymax = 300, 800
                    elif instrument == 'IMACS2':
                        xmin, xmax = 65, 2110
                        ymin, ymax = 150, 400
                    elif instrument == 'BlueChannel':
                        xmin, xmax = 5, 2700
                        ymin, ymax = 100, 250
                    elif instrument == 'LDDS3':
                        xmin, xmax = 650, 4095
                        ymin, ymax = 320, 620
                    elif instrument == 'LDDS3_2x2':
                        xmin, xmax = 300, 2047
                        ymin, ymax = 225, 380
                    elif instrument == 'Goodman':
                        xmin, xmax = 65, 2067
                        ymin, ymax = 100, 750
                    elif instrument == 'Kosmos':
                        xmin, xmax = 1, 4096
                        ymin, ymax = 340, 1600

                    # Crop the data
                    fits_file = fits.open(crop_name, ignore_missing_end=True)
                    fits_file[data_index].data = fits_file[data_index].data[ymin:ymax,:]
                    fits_file.writeto(crop_name, overwrite = True)
                    fits.setval(crop_name,   biassec_key,  value='[%s:%s,%s:%s]'%(xmin - xmin + 1, xmax - xmin, ymin - ymin + 1, ymax - ymin), ext = header_index)
                    fits.setval(crop_name,   datasec_key ,  value='[%s:%s,%s:%s]'%(xmin - xmin + 1, xmax - xmin, ymin - ymin + 1, ymax - ymin), ext = header_index)
                    fits.setval(crop_name,  'CCDSEC'    ,  value='[%s:%s,%s:%s]'%(xmin - xmin + 1, xmax - xmin, ymin - ymin + 1, ymax - ymin), ext = header_index)
                    fits.setval(crop_name,  'DISPAXIS', value='1', ext = header_index)
                    print('Cropped ' + crop_name)

def prepare_binospec(file_directory = 'raw_data/*.fits', crop = True, flip = False, variables = [''], disperser_name = 'x270', header_extension = 1, data_extension = 2):
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
            output_image_name = input_image_name[:input_image_name.find('sci')] + input_image_name[input_image_name.find('g_20')+10:]
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

def pre_prepare_data(file_directory = 'original/*.fits.fz'):
    # Import file names
    Files = sorted(glob.glob(file_directory))

    for i, name in enumerate(Files):
        print(i, name)
        A = fits.open(name)
        hdu = fits.PrimaryHDU(data = A[1].data, header = A[1].header)
        hdul = fits.HDUList([hdu])
        hdul.writeto('raw_data/' + name.split('/')[1])

##### IMACS #####
#extract_fits_info('raw_data/*c8.fits', ['OBJECT', 'EXPTYPE', 'EXPTIME', 'RA', 'DEC', 'DATE-OBS', 'TIME-OBS', 'FILTER', 'DISPERSR', 'BINNING', 'AIRMASS'])
#prepare_data('raw_data/*c8.fits', variables = [''], rotate = True, crop = True, break_character = ' ', instrument = 'IMACS1')

#### Binospec ####
#extract_fits_info('raw_data/*.fits', ['OBJECT', 'IMAGETYP', 'SCRN', 'EXPTIME', 'RA', 'DEC', 'DATE-OBS', 'FILTER', 'MASK', 'DISPERS1', 'DISPERS2', 'HENEAR', 'MJD', 'AIRMASS', 'EXPMODE', 'PI'], header_index = 1, data_index = 1, return_counts = False)
#prepare_binospec()

#### Binospec 600 ####
#extract_fits_info('raw_data/*.fits', ['OBJECT', 'IMAGETYP', 'SCRN', 'EXPTIME', 'RA', 'DEC', 'DATE-OBS', 'FILTER', 'MASK', 'DISPERS1', 'DISPERS2', 'HENEAR', 'MJD', 'AIRMASS', 'EXPMODE', 'PI'], header_index = 1, data_index = 1, return_counts = False)
#prepare_binospec(disperser_name = 'x600')

#### FAST ####
#extract_fits_info('raw_data/*.fits', ['OBJECT', 'EXPTIME', 'RA', 'DEC', 'DATE', 'APERTURE', 'DISPERSE'])
#prepare_data(variables = [''], rotate = False, crop = False, break_character = '', disperser = 'DISPERSE', filter_name = '')

#### WHT ####
#extract_fits_info('raw_data/*.fits', ['OBSTYPE',  'OBJECT', 'EXPTIME', 'RA', 'DEC', 'DATE-OBS', 'UT', 'AIRMASS', 'ISISLITW', 'DISPERSI', 'ISIFILTA', 'ISIFILTB', 'ISIGRAT', 'ISIARM', 'CENWAVE'])
#prepare_data(variables = ['ISIARM'], rotate = True, crop = True, break_character = ' ', disperser = 'ISIGRAT', filter_name = '', datasec_key = 'RTDATSEC', data_index = 1, rotations = 1)

#### LDSS ####
#extract_fits_info('raw_data/*c1.fits', ['OBJECT', 'EXPTYPE', 'EXPTIME', 'BINNING', 'RA', 'DEC', 'DATE-OBS', 'TIME-OBS', 'FILTER', 'GRISM'])
#prepare_data('raw_data/*c1.fits', variables = [''], rotate = True, crop = True, break_character = ' ', disperser = 'GRISM', filter_name = 'Open', rotations = 1, instrument = 'LDDS3')

#### Blue Channel ####
#extract_fits_info('raw_data/*.fits', ['IMAGETYP', 'OBJECT', 'EXPTIME', 'RA', 'DEC', 'DATE-OBS', 'UT', 'AIRMASS', 'APERTURE', 'DISPERSE', 'CENWAVE'])
#prepare_data(variables = [''], rotate = False, crop = True, break_character = ' ', disperser = 'DISPERSE', filter_name = '', instrument = 'BlueChannel')

#### SOAR Goodman ####
#extract_fits_info('original/*.fits.fz', ['OBSTYPE', 'OBJECT', 'EXPTIME', 'RA', 'DEC', 'DATE-OBS', 'UT', 'AIRMASS', 'FILTER', 'FILTER2', 'GRATING', 'SLIT'], header_index = 1, data_index = 1)
#pre_prepare_data('original/*.fits.fz')
#prepare_data('raw_data/*.fits.fz', variables = [''], rotate = False, crop = True, break_character = ' ', disperser = 'GRATING', filter_name = '', instrument = 'Goodman', header_index = 0, data_index = 0)

#### APO Kosmos ####
#extract_fits_info('raw_data/*.fits', ['OBJNAME','RA','DEC','DATE-OBS','IMAGETYP','EXPTIME','DISPERSR','QUARTZ','NEON','SLIT','FILTER2','FILTER1'])
#prepare_data(variables = ['DISPERSR'], rotate = True, crop = True, disperser = 'DISPERSR', instrument = 'Kosmos', filter_name = '', objname = 'IMAGETYP', datasec_key = 'CSEC11', biassec_key = 'BSEC11')


