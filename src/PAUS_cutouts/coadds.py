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


def get_images_info(RA, DEC, NB_wav_Arr, search_size=0.05,
                    fully_contained_only=False,
                    save_metadata=False):
    # Establish connection with database
    dsn = 'postgresql://readonly:PAUsc1ence@db.pau.pic.es/dm'
    engine = sqla.create_engine(dsn)

    # Search in a larger square because images are not exact
    # rectangles in RA, DEC
    ra_min = RA - search_size
    ra_max = RA + search_size
    dec_min = DEC - search_size
    dec_max = DEC + search_size

    df_list = []
    for NB_wav in np.atleast_1d(NB_wav_Arr):
        if fully_contained_only:
            coordinate_conditions = f"""({ra_min} > i.ra_min AND {ra_max} < i.ra_max)
            AND ({dec_min} > i.dec_min AND {dec_max} < i.dec_max)"""
            
        else:
            coordinate_conditions = f"""NOT ({ra_min} >= i.ra_max OR {ra_max} <= i.ra_min)
            AND NOT ({dec_min} >= i.dec_max OR {dec_max} <= i.dec_min)"""

            
        query = f"""SELECT i.archivepath, i.filename, i.zp_nightly,
                i.ra_min, i.ra_max, i.dec_min, i.dec_max, m.exp_time
                FROM image as i
                JOIN mosaic as m
                ON i.mosaic_id = m.id
                WHERE {coordinate_conditions}
                AND i.wavelength={NB_wav}
                AND m.production_id=943
                AND m.kind='RED_SCI'"""
    
        df_list.append(pd.read_sql(query, engine))
        df = pd.concat(df_list)
        if save_metadata:
            df.to_csv('cutout_metadata.csv')

    return df


def generate_image_list(df, save_path, save_exp_time=False):
    os.makedirs(save_path, exist_ok=True)
    
    img_list = []
    zp_list = []
    exp_time_list = []
    for _, row in df.iterrows():
        img_list.append(row.filename)
        zp_list.append(row.zp_nightly)
        exp_time_list.append(row.exp_time)
                
    with open(f'{save_path}/img_list.txt', 'w') as writer:
        [writer.write(f'{save_path}/{img_name}\n') for img_name in img_list]
        
    with open(f'{save_path}/zero_points.txt', 'w') as writer:
        [writer.write(f'{zp}\n') for zp in zp_list]
        
    if save_exp_time:
        with open(f'{save_path}/exp_time.txt', 'w') as writer:
            [writer.write(f'{exp_time}\n') for exp_time in exp_time_list]
        

def crop_images(df, RA, DEC, cutout_square_size, savepath,
                suffix='', only_contained=True):
    excluded_images = []
    for i, (archivepath, fname) in enumerate(zip(df.archivepath, df.filename)):
        fname = '.'.join(fname.split('.')[:-1]) + f'{suffix}.fits'
        
        hdul = fits.open(f'{archivepath}/{fname}')
        img = hdul[0].data
        coords = SkyCoord(RA, DEC, unit='deg')
        wcs = WCS(hdul[0])
        
        # Check if (RA,DEC) position is contained in the wcs of the image
        if not coords.contained_by(wcs, image=img) and only_contained:
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
    def inner(*args, **kwargs):
        if 'single_epoch_dir' in kwargs.keys():
            single_epoch_dir = kwargs['single_epoch_dir']
        elif len(args) > 5:
            single_epoch_dir = args[5]
        else:
            single_epoch_dir = 'single_epoch_cutouts'
            
        # Never set exist_ok=True or very bad things could happen !!
        os.makedirs(single_epoch_dir, exist_ok=False)
        
        try:
            func(*args, **kwargs)
            shutil.rmtree(single_epoch_dir)
        except:
            shutil.rmtree(single_epoch_dir)
            raise            
        
    return inner


#@rm_tmp_dir
def generate_cutouts(RA_Arr, DEC_Arr, ID_Arr, square_size,
                             NB_wav_Arr, single_epoch_dir='single_epoch_cutouts',
                             save_coadds_dir='coadd_cutouts',
                             config_template='config.swarp',
                             combine_type='AVERAGE', 
                             weight_type='BACKGROUND',
                             delete_single_epoch=True,
                             coadd=True,
                             only_contained=True,
                             save_exp_time=False):
    """
    Generate cutouts for given RA, DEC coordinates and ID with specific square size and narrowband wavelength. May coadd the cutouts using SWarp.

    Parameters:
    RA_Arr : Array of Right Ascension (RA) coordinates in deg.
    DEC_Arr : Array of Declination (DEC) coordinates in deg.
    ID_Arr : Array of IDs of objects used to generate file names.
    square_size : Size of the square cutout in deg. Either float or list of floats.
    NB_wav_Arr : Array of narrowband wavelength identificator (integer).
    single_epoch_dir : Directory path for single-epoch cutouts. Defaults to 'single_epoch_cutouts'.
    save_coadds_dir : Directory path to save coadds. Defaults to 'coadd_cutouts'.
    config_template : Template for Swarp config file. Defaults to 'config.swarp'.
    combine_type: COMBINE_TYPE keyword of config.swarp. Defaults to AVERAGE
    combine_type: WEIGHT_TYPE keyword of config.swarp. Defaults to BACKGROUND
    delete_single_epoch: Delete the single_epoch_dir directory after coadding. Default is False
    only_contained: Keep only cutouts that contain the exact coordinates of the object. Default is True
    coadd: Coadd the single-epoch cutouts using SWarp. Default is True
    save_exp_time: Save in a .ttxt file exposure time of every single-epoch cutout. Default is False

    NOTE: NB_wav_Arr can be:
        - Single integer: Make coadds for single NB
        - List of integers: Make coadds for all the NBs in the list
        - List of lists of integers: Make coadds for the combined NBs in each nested list.

    Returns:
    None
    """
    RA_Arr = np.atleast_1d(RA_Arr)
    DEC_Arr = np.atleast_1d(DEC_Arr)
    ID_Arr = np.atleast_1d(ID_Arr)
    try:
        len(square_size)
    except TypeError:
        square_size = [square_size] * len(ID_Arr)
    
    if isinstance(NB_wav_Arr, int):
        NB_wav_Arr = np.atleast_1d(NB_wav_Arr)
            
    for NB_wav in NB_wav_Arr:                
        for RA, DEC, ID, square_size_tmp in zip(RA_Arr, DEC_Arr, ID_Arr, square_size):
            print(f'\n\n {RA=}, {DEC=}, {ID=}\n')
            df = get_images_info(RA, DEC, NB_wav)

            # Copy images to temporary directory
            generate_image_list(df, single_epoch_dir, save_exp_time=save_exp_time)

            # Crop images. Weight maps will always be downloaded, just in case
            excluded_images = crop_images(df, RA, DEC, square_size_tmp,
                                          single_epoch_dir, only_contained=only_contained)
            crop_images(df, RA, DEC, square_size_tmp, single_epoch_dir,
                        suffix='.weight', only_contained=only_contained)

            # Read zero point list file to check if there is any nans
            with open(f'{single_epoch_dir}/zero_points.txt', 'r') as file:
                    zp_list_lines = file.readlines()
            for i, line in enumerate(zp_list_lines):
                if line[:3] == 'nan':
                    excluded_images.append(i)

            # Remove excluded images from list
            if len(excluded_images) > 0:
                with open(f'{single_epoch_dir}/img_list.txt', 'r') as file:
                    img_list_lines = file.readlines()
                with open(f'{single_epoch_dir}/img_list.txt', 'w') as file:
                    for i, line in enumerate(img_list_lines):
                        if i not in excluded_images:
                            file.write(line)
                
                with open(f'{single_epoch_dir}/zero_points.txt', 'w') as file:
                    for i, line in enumerate(zp_list_lines):
                        if i not in excluded_images:
                            file.write(line)
                            
                if save_exp_time:
                    with open(f'{single_epoch_dir}/exp_time.txt', 'r') as file:
                        exp_time_lines = file.readlines()
                    with open(f'{single_epoch_dir}/exp_time.txt', 'w') as file:
                        for i, line in enumerate(exp_time_lines):
                            if i not in excluded_images:
                                file.write(line)
                            
            ### Generate new Swarp config file from template
            if coadd:
                if len(np.atleast_1d(NB_wav)) > 1:
                    save_coadds_to = f'{save_coadds_dir}/NB' + '_'.join(NB_wav.astype(str))
                    os.makedirs(save_coadds_to, exist_ok=True)
                else:
                    save_coadds_to = f'{save_coadds_dir}/NB{int(NB_wav)}'
                    os.makedirs(save_coadds_to, exist_ok=True)
                
                with open(config_template, 'r') as file:
                    filedata = file.read()
                    out_filename = f'{save_coadds_to}/coadd_cutout_{ID}.fits'
                    filedata = filedata.replace('coadd.fits', out_filename)
                    out_filename_w = f'{save_coadds_to}/coadd_cutout_{ID}.weight.fits'
                    filedata = filedata.replace('coadd.weight.fits', out_filename_w)
                    out_filename_w = f'{save_coadds_to}/coadd_cutout_{ID}.weight.fits'
                    filedata = filedata.replace('combinetype_cutouts', combine_type)
                    filedata = filedata.replace('weighttype_cutouts', weight_type)
                    filedata = filedata.replace('tmpfiledir', single_epoch_dir)

                with open(f'{single_epoch_dir}/config.swarp', 'w') as file:
                    file.write(filedata)

                # Run swarp
                img_list_path = f'{single_epoch_dir}/img_list.txt'
                zp_list_path = f'{single_epoch_dir}/zero_points.txt'
                config_file_path = f'{single_epoch_dir}/config.swarp'
                swarp_path = '/data/astro/software/centos7/swarp/2.41.5/bin/swarp'

                os_out = os.system(f'{swarp_path} @{img_list_path} -c {config_file_path} -FSCALE_DEFAULT @{zp_list_path} -FSCALE_KEYWORD "nokeyword"')
                if os_out != 0:
                    raise Exception(f'Error {os_out}')

            # Delete tmp files for this coadd
            if delete_single_epoch:
                shutil.rmtree(single_epoch_dir)
                #tmp_files = glob.glob(f'{single_epoch_dir}/*')
                #for f in tmp_files:
                #    os.remove(f)
                    
            else:
                shutil.move(f'{single_epoch_dir}', f'{single_epoch_dir}_{ID}/{NB_wav}')

