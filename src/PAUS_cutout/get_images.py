import os
import shutil
import numpy as np
import pandas as pd
import sqlalchemy as sqla

import matplotlib.pyplot as plt

from astropy.io import fits

def get_image_list(RA, DEC, RA_size, DEC_size, NB_wavelength_list,
                   fully_contained_only=False):
    '''
    Reads a list of RA,DEC of objects and generates pandas dataframes with info of the
    images.
    '''

    # Establish connection with database
    dsn = 'postgresql://readonly:PAUsc1ence@db.pau.pic.es/dm'
    engine = sqla.create_engine(dsn)

    try:
        len(NB_wavelength_list)
    except TypeError:
        wl_nbs = [wl_nbs]

    RA_list = np.atleast_1d(RA)
    DEC_list = np.atleast_1d(DEC)
    
    # Initialize a couple of arrays
    archivepath_list = []
    filename_list = []

    # Structure of dfs: for each source, stores dictionary for each filter in wl_nb
    dfs = []

    for ra_target, dec_target in zip(RA_list, DEC_list):
        ra_min = ra_target - RA_size/2
        ra_max = ra_target + RA_size/2
        dec_min = dec_target - DEC_size/2
        dec_max = dec_target + DEC_size/2

        if fully_contained_only:
            coordinate_conditions = f"""({ra_min} > i.ra_min AND {ra_max} < i.ra_max)
            AND ({dec_min} > i.dec_min AND {dec_max} < i.dec_max)"""

        else:
            coordinate_conditions = f"""NOT ({ra_min} >= i.ra_max OR {ra_max} <= i.ra_min)
            AND NOT ({dec_min} >= i.dec_max OR {dec_max} <= i.dec_min)"""


        this_df_list = []
        this_archivepath_list = []
        this_filename_list = []
        for wl_nb in wl_nbs:
            query = f"""SELECT i.archivepath, i.filename, i.zp_nightly,
                    i.ra_min, i.ra_max, i.dec_min, i.dec_max
                    FROM image as i
                    JOIN mosaic as m
                    ON i.mosaic_id = m.id
                    WHERE {coordinate_conditions}
                    AND i.wavelength={wl_nb}
                    AND m.production_id=943
                    AND m.kind='RED_SCI'"""
            this_df = pd.read_sql(query, engine)
            this_df_list.append(this_df)
        
        dfs.append(this_df_list)

    return dfs


def download_images(dfs, wl_nbs, img_folder):
    '''
    Downloads the images and stores them in separate folders by NB.
    '''
    os.makedirs(f'{img_folder}', exist_ok=True)

    for j, df_list in enumerate(dfs):
        for i, wavelength in enumerate(wl_nbs):
            save_path = f'{img_folder}/{wavelength}'
            img_list = []
            zp_list = []
            os.makedirs(save_path, exist_ok=True)

            df = df_list[i]

            for _, row in df.iterrows():
                load_dir = f'{row.archivepath}/{row.filename}'
                shutil.copyfile(load_dir, f'{save_path}/{row.filename}')
                img_list.append(row.filename)
                zp_list.append(row.zp_nightly)
                
            with open(f'{save_path}/img_list.txt', 'w') as writer:
                [writer.write(f'{img_name}\n') for img_name in img_list]
                
            with open(f'{save_path}/zero_points.txt', 'w') as writer:
                [writer.write(f'{zp}\n') for zp in zp_list]
