# PAUS_cutouts
Python package for generating coadded cutouts of PAUS images. This package is designed to run in the PIC server (https://www.pic.es/).

## Installation
- Python >=3.11 recommended.
-  Install the dependencies: pandas, numpy, sqlalchemy, astropy
- Clone this repository:

		git clone https://github.com/totoal/PAUS_cutouts
 
- Install the package with
  
		cd PAUS_cutouts
  		pip install -e .

## Use
The main function can be imported as

	from PAUS_cutout.coadds import generate_coadded_cutouts

See docstring for more info.

A template of config.swarp is included in this repository. Copy it to your working directory, or use your own template.
