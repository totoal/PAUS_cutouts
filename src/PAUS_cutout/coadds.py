import os
import sys
import shutil
import glob

import pandas as pd
import sqlalchemy as sqla
import numpy as np

from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
from astropy.nddata.utils import NoOverlapError


def get_images_info(RA, DEC, square_size, NB_wavelength,
                    fully_contained_only=False):
    # Establish connection with database
    dsn = 'postgresql://readonly:PAUsc1ence@db.pau.pic.es/dm'
    engine = sqla.create_engine(dsn)

    # Search in a larger square because images are not exact
    # rectangles in RA, DEC
    ra_min = RA - square_size
    ra_max = RA + square_size
    dec_min = DEC - square_size
    dec_max = DEC + square_size

    if fully_contained_only:
        coordinate_conditions = f"""({ra_min} > i.ra_min AND {ra_max} < i.ra_max)
        AND ({dec_min} > i.dec_min AND {dec_max} < i.dec_max)"""

    else:
        coordinate_conditions = f"""NOT ({ra_min} >= i.ra_max OR {ra_max} <= i.ra_min)
        AND NOT ({dec_min} >= i.dec_max OR {dec_max} <= i.dec_min)"""

    query = f"""SELECT i.archivepath, i.filename, i.zp_nightly,
            i.ra_min, i.ra_max, i.dec_min, i.dec_max
            FROM image as i
            JOIN mosaic as m
            ON i.mosaic_id = m.id
            WHERE {coordinate_conditions}
            AND i.wavelength={NB_wavelength}
            AND m.production_id=943
            AND m.kind='RED_SCI'"""
    
    df = pd.read_sql(query, engine)

    return df

def copy_images_to_home(df, save_path):
    os.makedirs(save_path, exist_ok=True)
    
    img_list = []
    zp_list = []
    for _, row in df.iterrows():
        img_list.append(row.filename)
        zp_list.append(row.zp_nightly)
        
    with open(f'{save_path}/img_list.txt', 'w') as writer:
        [writer.write(f'{save_path}/{img_name}\n') for img_name in img_list]
        
    with open(f'{save_path}/zero_points.txt', 'w') as writer:
        [writer.write(f'{zp}\n') for zp in zp_list]

        

def crop_images(df, RA, DEC, cutout_square_size, savepath,
                suffix=''):
    excluded_images = []
    for i, (archivepath, fname) in enumerate(zip(df.archivepath, df.filename)):
        fname = '.'.join(fname.split('.')[:-1]) + f'{suffix}.fits'
        
        hdul = fits.open(f'{archivepath}/{fname}')
        img = hdul[0].data
        coords = SkyCoord(RA, DEC, unit='deg')
        wcs = WCS(hdul[0])
        
        # Check if (RA,DEC) position is contained in the wcs of the image
        if not coords.contained_by(wcs, image=img):
            print(f'Skipping {fname}: Coordinates not contained in the image.')
            excluded_images.append(i)
            continue
            
        try:
            cutout = Cutout2D(img, coords, size=cutout_square_size * u.deg,
                              wcs=wcs, mode='trim')
            cutout_hdu = hdul[0]
            cutout_hdu.data = cutout.data
            cutout_hdu.header.update(cutout.wcs.to_header())
            cutout_hdu.writeto(f'{savepath}/{fname}', overwrite=True)
        except NoOverlapError:
            print(f'Skipping {fname}: no overlap with desired area.')
            excluded_images.append(i)
            continue
        
        if not np.all(cutout.data.shape):
            print(f'Skipping {fname}: Cutout has zero dimension.')
            excluded_images.append(i)
    
    return excluded_images


def rm_tmp_dir(func):
    def inner(*args):
        tmp_files_dir = args[-1]
        # Never set exist_ok=True or very bad things could happen !!
        os.makedirs(tmp_files_dir, exist_ok=False)
        
        try:
            func(*args)
        except:
            shutil.rmtree(tmp_files_dir)
            raise
        finally:
            shutil.rmtree(tmp_files_dir)
        
    return inner

@rm_tmp_dir
def generate_coadded_cutouts(RA_Arr, DEC_Arr, ID_Arr, square_size,
                             NB_wav_Arr, tmp_files_dir='tmp_files_cutouts',
                             save_coadds_dir='out_cutouts',
                             config_template='config.swarp'):
    for NB_wav in NB_wav_Arr:
        save_coadds_to = f'{save_coadds_dir}/NB{int(NB_wav)}'
        os.makedirs(save_coadds_to, exist_ok=True)
        
        for RA, DEC, ID in zip(RA_Arr, DEC_Arr, ID_Arr):
            print(f'\n\n {RA=}, {DEC=}, {ID=}\n')
            df = get_images_info(RA, DEC, square_size, NB_wav)

            # Copy images to temporary directory
            copy_images_to_home(df, tmp_files_dir)

            # Crop images
            excluded_images = crop_images(df, RA, DEC,cutout_square_size,
                                          tmp_files_dir)
            crop_images(df, RA, DEC, cutout_square_size, tmp_files_dir,
                        suffix='.weight')
            
            # Remove excluded images from list
            if len(excluded_images) > 0:
                with open(f'{tmp_files_dir}/img_list.txt', 'r') as file:
                    img_list_lines = file.readlines()
                    
                with open(f'{tmp_files_dir}/img_list.txt', 'w') as file:
                    for i, line in enumerate(img_list_lines):
                        if i not in excluded_images:
                            file.write(line)
            ### Generate new Swarp config file from template
            with open(config_template, 'r') as file:
                filedata = file.read()
                
            out_filename = f'{save_coadds_to}/coadd_cutout_{ID}_NB{NB_wav}.fits'
            filedata = filedata.replace('coadd.fits', out_filename)
            out_filename_w = f'{save_coadds_to}/coadd_cutout_{ID}_NB{NB_wav}.weight.fits'
            filedata = filedata.replace('coadd.weight.fits', out_filename_w)
           
            with open(f'{tmp_files_dir}/config.swarp', 'w') as file:
                file.write(filedata)

            # Run swarp
            img_list_path = f'{tmp_files_dir}/img_list.txt'
            config_file_path = f'{tmp_files_dir}/config.swarp'
            swarp_path = '/data/astro/software/centos7/swarp/2.41.5/bin/swarp'

            os_out = os.system(f'{swarp_path} @{img_list_path} -c {config_file_path}')
            if os_out != 0:
                raise Exception(f'Error {os_out}')
                
            # Delete tmp files for this coadd
            tmp_files = glob.glob(f'{tmp_files_dir}/*')
            for f in tmp_files:
                os.remove(f)
                