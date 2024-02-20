import os
import shutil
import numpy as np
import pandas as pd
import sqlalchemy as sqla


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
        NB_wavelength_list = [NB_wavelength_list]

    RA_list = np.atleast_1d(RA)
    DEC_list = np.atleast_1d(DEC)

    # Structure of dfs: for each source, stores dictionary for each filter in wl_nb
    dfs = []

    for iii, (ra_target, dec_target) in enumerate(zip(RA_list, DEC_list)):
        print(f'Retrieving image info: {iii + 1} / {len(RA_list)}',
                end='\r')

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
        for wl_nb in NB_wavelength_list:
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


def download_images(dfs, NB_wavelength_list, img_folder):
    '''
    Downloads the images and stores them in separate folders by NB.
    '''
    os.makedirs(f'{img_folder}', exist_ok=True)

    for _, df_list in enumerate(dfs):
        for i, wavelength in enumerate(NB_wavelength_list):
            save_path = f'{img_folder}/{wavelength}'
            img_list = []
            zp_list = []
            os.makedirs(save_path, exist_ok=True)

            df = df_list[i]

            for _, row in df.iterrows():
                load_dir = f'{row.archivepath}/{row.filename}'
                
                save_to = f'{save_path}/{row.filename}'
                if not os.path.exists(save_to):
                    print(f'Copy: {row.filename}')
                    shutil.copyfile(load_dir, save_to)
                else:
                    print(f'Skip: {row.filename}')
                img_list.append(row.filename)
                zp_list.append(row.zp_nightly)
                
            with open(f'{save_path}/img_list.txt', 'w') as writer:
                [writer.write(f'{img_name}\n') for img_name in img_list]
                
            with open(f'{save_path}/zero_points.txt', 'w') as writer:
                [writer.write(f'{zp}\n') for zp in zp_list]

def get_images(RA, DEC, RA_size, DEC_size, NB_wavelength_list,
               savedir='PAUS_images', ID=None):
    if len(RA) != len(DEC):
        raise Exception('RA and DEC arrays must have the same size.')
    
    # If no ID list is given, make up ids
    if ID is None:
        ID = np.arange(len(RA))

    dfs = get_image_list(RA, DEC, RA_size, DEC_size, NB_wavelength_list)

    # Download the images
    download_images(dfs, NB_wavelength_list, savedir)

if __name__ == '__main__':
    NB_wavelength_list = [455, 465, 475, 485, 495, 505,
                          515, 525, 535, 545, 555, 565,
                          575, 585, 595, 605, 615, 625,
                          635]
    LAE_selection = pd.read_csv('~/LAE_selection_200224.csv')
    RA = LAE_selection.RA
    DEC = LAE_selection.DEC

    RA_size = 10 / 3600
    DEC_size = 10 / 3600

    dfs = get_image_list(RA, DEC, RA_size, DEC_size, NB_wavelength_list) # testing
    download_images(dfs, NB_wavelength_list, '~/PAUS_images')
