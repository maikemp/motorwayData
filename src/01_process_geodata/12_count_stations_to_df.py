# -*- coding: utf-8 -*-
"""


Runtime: approx. 2 Minutes.
"""

import arcpy
import os
import pandas as pd
import numpy as np
from arcgis import GIS

import src.process_geodata.functions.f12_countst_functions as count
from project_paths import get_project_paths, get_params

arcpy.env.overwriteOutput = True


pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))
sr = arcpy.SpatialReference(param['spacial_reference'])
sr_in = arcpy.SpatialReference('WGS 1984')

dates = param['dates']

compass = {( 22.5,  67.5): 'NO', ( 67.5, 112.5): 'O',
           (112.5, 157.5): 'SO', (157.5, 202.5): 'S',
           (202.5, 247.5): 'SW', (247.5, 292.5): 'W',
           (292.5, 337.5): 'NW', (337.5, 360  ): 'N', (0, 22.5):'N'}


for year in sorted([str(x) for x in param['year']]):
    
    print('Start with year:', year)
    pp_in, pp_temp_df, pp_temp_z, pp_frame, pp_out = count.make_file_paths(
                                                                      year, pp)

    # Here: first load layer as csv and change the column names so that they are
    # still understandable after being shortened. (Ri1 und Ri2 müssen erhalten bleiben)
    df = pd.read_csv(pp_in, sep=';', decimal=",", 
                     thousands = ".", encoding="ISO-8859-1")

    # Select only columns of interest:
    keep_col = ['Abschnitt_Ast', 'Nahziel_Ri1', 'Nahziel_Ri2',
                'DZ_Nr', 'Hi_Ri1', 'Hi_Ri2', 'Anz_Fs_Q',
                'Koor_WGS84_N', 'Koor_WGS84_E', 'fer', 'Mt', 'Mn', 'pMn',  
                'DTV_Kfz_MobisSo_Ri1', 'DTV_Kfz_MobisSo_Ri2', 
                'DTV_SV_MobisSo_Ri1', 'DTV_SV_MobisSo_Ri2',
                'MSV50_Kfz_MobisSo_Ri1', 'MSV50_Kfz_MobisSo_Ri2',
                'bSV50_MobisSo_Ri1', 'bSV50_MobisSo_Ri2', 
                'DL_Ri1', 'DL_Ri2', 'bSo_Ri1', 'bSo_Ri2', 'bMo_Ri1', 'bMo_Ri2']
    
    # Subset and rename columns and save out as csv to pp_temp_df.
    df = count.select_and_rename_columns(df, keep_col)
    # as a place holder, assign -9999 for missings. Change back to nan later 
    # on. Reason: shapefile will otherwise contain missings for nans!
    df.fillna(-9999, inplace=True)
    column_list = list(df.columns)
    df.to_csv(pp_temp_df,sep=';',decimal=",",encoding="ISO-8859-1",index=False)
    del(df)
    
    # Make a feature layer from the csv file and project into UTM
    arcpy.management.MakeXYEventLayer(pp_temp_df, "KoorWGS84E", "KoorWGS84N", 
                                      r"in_memory\count_stations", sr_in, None)
    arcpy.management.Project(r"in_memory\count_stations", pp_temp_z, sr)
    # Delete layer with wrong projection to avoid confusion.
    arcpy.Delete_management(r"in_memory\count_stations")
    
    # Make a layer containing directional mean and the respecive compass value.
    # results will be saved in r"in_memory\frame".
    pp_frame_dir = r"in_memory\frame"
    count.add_dir_mean(pp_frame, compass, id_field='seg_ID', dissolve=True, 
                 make_layer=True, memory_name=pp_frame_dir, save_out=True, 
                 path_out=os.path.join(pp['data_temp'], 'frame_nosw'))
    # Drop counting stations with missing values if they are in the same 
    # sector as another station. Save to r"in_memory\counting_stations". 
    count.drop_double_stations(pp_temp_z, r"in_memory\counting_stations", 
                               os.path.join(pp['data_temp_gdb'], 'cs_temp'))

    # Make a dataframe based on pp_frame_dir that additionally contains the
    # values of a counting station in the given proximity. Drop coordinates.
    df = count.cs_to_id_df(r"in_memory\counting_stations", pp_frame_dir, 
                                  '255 Meters', column_list)
    df.drop(['KoorWGS84N', 'KoorWGS84E'], axis='columns', inplace=True)
    
    # Select only frame segments that have a counting station assigned and keep 
    # columns ending on R1 or R2 only for the direction the row corresponds to.
    df = count.select_zrows_and_direction(df, 
                                            del_cols=['compass', 'direction'])


    # Dictionary with individual directions for entire frame. Use as master df
    # to include ids and directions of all frame segments in  the output file.
    fr_dir_dict = count.add_dir_mean(pp_frame, compass, "ID", dissolve=False, 
                                       make_layer=False, save_out=False)
    df_master = pd.DataFrame.from_dict(
            fr_dir_dict, orient='index',columns=['compass', 'direction'])
    
    df_master = df_master.join(df, how='left')
    df_master = df_master.replace(-9999, np.nan).replace('-9999', np.nan)
    df_master.to_csv(pp_out, encoding="ISO-8859-1")











