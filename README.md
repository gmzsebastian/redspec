# redspec
IRAF based pipeline for single slit spectra reduction.
Tarraneh picked the name, she wants credit.

## Usage

The pipeline has 5 different scripts to be run separately. You might be able to run them together depending on whether PyRAF and Python like working with each other, they usually don't:


0. Get the data into the proper format for the pipeline
1. Reduce the data (bias, flats, etc.)
2. Do the wavelength correction interactively
3. Apply a flux calibration to a spectroscopy standard star
4. Plot the result and save the output ASCII files

## 0: Prepare Data

This is the main file you will have to edit depending on what telescope you are using, some example files are in the `prepare_files` folder. You can start by making a folder with all the raw data, i.e. `raw_data`, then this script will take those files, read in the headers, generate a logfile for the night, and separate the targets into individual folders. You want to end up with a directory structure that looks like this:

```
|-- raw_data
    |-- 1.fits
    |-- 2.fits
    |-- 3.fits
    |-- 4.fits
    ...

|-- TargetA
    |-- HeNeAr_1.fits
    |-- Flat_2.fits
    |-- TargetA_3.fits

|-- TargetB
    |-- HeNeAr_4.fits
    |-- Flat_5.fits
    |-- TargetB_6.fits

|-- bias
    |-- bias_7.fits
    |-- bias_8.fits
    |-- bias_9.fits
```
The script will append the file type at the beginning of each file, this is necessary for the rest of the pipeline to work.
If you want rotate, flip, or crop ALL the files, here is where you can specify this. Also, the `variables` parameter in the `prepare_data()` function serves to add suffixes to the target name, if for example you observed one target with two different configurations.

## 1: Reduce Data
This script will take in the data in the format produced by `0_Prepare_Data.py`: create the bias, create the flats, and reduce the data with `ccdproc` for bias and flats. If you only have one set of flats for the entire night, use the `individual_flats()` function first before running `reduce_data()`. You can also modify the default orders of the functions being fit, or change the sameple region from everything(`*`) or a spcific region (`[50:200]`).

## 2: Correct Wavelength
This script will deterime and apply the wavelength calibration to your files. You need an initial guess of the central wavelength and pixel scale of the spectrograph to be fed into `test_solution()`. You will then specify a few bright lines that are easy to identify and interactively change the central wavelength and pixel scale values until you have a good fit.

Then you can feed those answers into the `wavelength_solution()` function, this will do the proper calibration with all the lines it manages to fit, and the order you specify (2, 3, or 4).

## 3: Flux Correction
If you took data of a spectroscopic standard that is in the IRAF calibration directories, extract it with this script. `create_standard_sens()` will create the sensitivity file for the standard. Then you can correct your target with this sensitivity file.

## 4: Plot Molly
The last script will take in the reduced data, extract it into individual files, convert those to ASCII files, read those in, average them, and plot them. You can specify if you want to use the flux calibrated data or not. The output file will by default use the optimally extracted data that has been cleaned for cosmic rays, but the non-optimally extracted data will also be plotted.


