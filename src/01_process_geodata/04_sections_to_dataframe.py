"""
Runtime: approx. 10 minutes
"""
import arcpy
import os

import src.process_geodata.functions.f04_sec_to_df_functions as sec_to_df
# from process_geodata.functions.general_functions import *
from project_paths import get_project_paths, get_params

arcpy.env.overwriteOutput = True


pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))

dates = param['dates']
# for date in dates....
for date in dates:
    
    print('Started creating data frame for:', date)
    # In paths:
    pp_change_points = os.path.join(
            pp['data_out_otsecpt'], 'change_points_dist'+date+'.shp'
        )
    pp_cut_points = os.path.join(
            pp['data_out_otsecpt'], 'cut_points_populated'+date+'.shp'
        )
    # List of used ids.
    id_path = os.path.join(pp['data_temp'], 'id_list'+date+'.txt')
    id_list = sec_to_df.intersect_ids(id_path, pp_cut_points)
    
    with open(
        os.path.join(pp['data_temp'], 'id_list' + date + '_df.txt'), "w"
    ) as f:
        f.write(str(id_list))
    
    # Find the maximum number of times, the same ID appears in change points.
    field_duplic_num = sec_to_df.duplic_val_num(pp_change_points, 'ID')
    
    excl_fs_cutpt = ['FID_12', 'FID', 'FID_1', 'in_seg_ID', 'official_r', 'Shape']
    cutpt_df = sec_to_df.make_cutpt_df(pp_cut_points, excl_fs_cutpt, id_list)
    
    
    excl_fs_changept = ['F_id','ref','FID', 'Shape', 'OBJECTID', 'TARGET_FID', 
                        'par_seg_ID', 'seg_ID', 'in_seg_ID', 'ORIG_FID',
                        'FID_1', 'COUNT', 'FID_12', 'BUFF_DIST']
    
    changept_df = sec_to_df.make_changept_df(
            pp_change_points, excl_fs_changept, id_list, field_duplic_num
        )
    
    
    # Concatenate dataframes
    full_data = cutpt_df.join(changept_df, how='outer')
    
    full_data = sec_to_df.remove_double_startpt(full_data)

    
    for col_name in ['bridge', 'tunnel']:
        full_data[
                    [x for x in full_data.columns if x.startswith(col_name)]
                 ].replace('missing','no', inplace=True)
    
    full_data.to_csv(os.path.join(pp['data_out_df'],'ot_df_'+date+'.csv'))
    print('Done')


