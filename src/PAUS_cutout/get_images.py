import os
import shutil
import numpy as np
import pandas as pd
import sqlalchemy as sqla

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
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


        dfs = []
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
            dfs.append(this_df)
        
            # Add image paths to the download list
            archivepath_list += list(this_df.archivepath)
            filename_list += list(this_df.filename)

    image_list = []
    for ap, fn in zip(archivepath_list, filename_list):
        image_list.append(f'{ap}/{fn}')

    # This is the list of paths of images to download
    unique_image_list = np.unique(image_list)

    return unique_image_list
