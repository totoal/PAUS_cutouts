import pandas as pd
import sqlalchemy as sqla
import os
import shutil

from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.units as u
from astropy.nddata.utils import NoOverlapError

def get_images_info(RA, DEC, RA_size, DEC_size, NB_wavelength,
                    fully_contained_only=False):
    # Establish connection with database
    dsn = 'postgresql://readonly:PAUsc1ence@db.pau.pic.es/dm'
    engine = sqla.create_engine(dsn)

    ra_min = RA - RA_size/2
    ra_max = RA + RA_size/2
    dec_min = DEC - DEC_size/2
    dec_max = DEC + DEC_size/2

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
        
    weight_filenames = ['.'.join(img_name.split('.')[:-1]) + '.weight.fits'
                        for img_name in img_list]
    with open(f'{save_path}/weight_list.txt', 'w') as writer:
        [writer.write(f'{fn}\n') for fn in weight_filenames]
        
    with open(f'{save_path}/zero_points.txt', 'w') as writer:
        [writer.write(f'{zp}\n') for zp in zp_list]

        

def crop_images(df, RA, DEC, cutout_square_size, savepath):
    for i, (archivepath, fname) in enumerate(zip(df.archivepath, df.filename)):
        hdul = fits.open(f'{archivepath}/{fname}')
        img = hdul[0].data
        coords = SkyCoord(RA, DEC, unit='deg')
        wcs = WCS(hdul[0])
        
        try:
            cutout = Cutout2D(img, coords, size=cutout_square_size * u.deg,
                              wcs=wcs, mode='trim')
        except NoOverlapError:
            print(f'Skipping {fname}, no overlap with desired area.')

        cutout_hdu = hdul[0]
        cutout_hdu.data = cutout.data
        cutout_hdu.header.update(cutout.wcs.to_header())
        cutout_hdu.writeto(f'{savepath}/{fname}', overwrite=True)
        
