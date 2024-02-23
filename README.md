# PAUS_cutouts
Python package for generating coadded cutouts of PAUS images. This package is designed to run in the PIC server (https://www.pic.es/).

## Installation
- Install the dependencies: pandas, numpy, sqlalchemy, astropy
- Clone this repository:

		git clone https://github.com/totoal/PAUS_cutouts
 
- Install the package with
  
		cd PAUS_cutouts
  		pip install -e .

## Use
The main function can be imported as

	from PAUS_cutout.coadds import generate_coadded_cutouts

See docstring for more info.
