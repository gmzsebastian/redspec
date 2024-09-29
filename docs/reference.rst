.. _reference:

Before you Start
================

The way ``RedSpec`` is written is a little peculiar, mostly to prevent ``PyRAF`` from fighting with ``matplotlib`` and other packages, which 
can lead to crashes. The package is made up of 5 individual scripts that need to be ran separately. This means you should run
the first script, quit out of its Python instance, and then run the next script.

Although ``RedSpec`` is designed to work only on single slit spectroscopic data, it can work on several different instruments.
Throughout the tutorials we provide different examples on how to run the scripts for different instruments.

You will notice that most scripts have a ``def_instrument()`` function that sets the instrument specific parameters. This might
not be needed, depending on the way IRAF is setup on your machine. However, if your instance crashes due to "No Instrument Parameter Specified"
you might need to uncomment the ``def_instrument()`` function.

Scripts
-------

The five scripts that make up the package are:

1. ``0_Prepare_Data.py``: This script will read in your raw data, sort it into the expected directory structure, and modify the files to be in the format that ``RedSpec`` expects.
2. ``1_Reduce_Data.py``: This script will reduce the data. This includes doing bias and flat field corrections, background subtraction, and trace extraction.
3. ``2_Correct_Wavelength.py``: This script will obtain the wavelength solution for the data.
4. ``3_Flux_Correction.py``: This script will correct the flux of the data based on a spectroscopic standard star.
5. ``4_PlotMolly.py``: This script plots the data.

Most functions in these scripts have the same two arguments. The first one is the name of the directory where the data is stored,
and the second one is the format of the object name in the files. Sometimes objects will be named with their actual name (e.g. "2022xzc"),
while other times they might be named with a generic name such as "spec".

Instruments
-----------

``RedSpec`` expects the data to be in a specific format, which means data from different instruments needs to be modified
to fit this expected format. Currently, the package has been tested on the following instruments:

- Magellan IMACS
- MMT Binospec
- FLWO FAST
- WHT ISIS
- Magellan LDSS3c
- SOAR Goodman
- APO Kosmos
- Gemini GMOS
