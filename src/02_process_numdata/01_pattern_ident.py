# -*- coding: utf-8 -*-
"""
Runtime: approx. 100 Minutes.
"""


import os 
import json

import src.process_numdata.functions.f01_pattern_functions as pattern 
from project_paths import get_project_paths, get_params


pp = get_project_paths()

param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))
dates = param['dates']
id_length = param['id_length']
sw_length = param['std_slipway_length']
max_link_length = param['max_link_length']

patt_param = json.load(open(os.path.join(pp['lib'], 'pattern_ident_lib.json')))


date = param['reference_date']

# First: mark links in lanes_long.
# Save out to path, so no need to rerun every time, takes little longer.
for date in dates:
    
    print('start with date:', date)
    
    pp_in_ot, pp_in_zaehl, pp_out_ot, pp_out_zaehl, pp_links_out, \
        pp_main_points_upd = pattern.get_pattern_paths(date, pp)
    # Read in data and do some cleaning.
    df_ot = pattern.read_data_select_cols_and_relabel(
        pp_in_ot, id_length, tuple(patt_param['remove_col_starts']), 
        tuple(patt_param['relabel_col_starts']), patt_param['relabel_dict']
    )
    
    # Extract lanes data in long format.
    lanes_long = pattern.get_long_ot_df(
            df_ot, 'lanes', ['Entry','Exit','Slipway_lanes','sw_type'])  
    
    ids = df_ot.index.to_list()
    last_ids = [i for i in ids if pattern.str_id(int(i)+1) not in ids]
    
    for sw_type in ['main', 'sec']:
        
        print('start with sw_type:', sw_type)

        link_long = pattern.get_and_prep_link_data(sw_type, pp, last_ids)
        
        # Update lanes_long to have identified links. 
        lanes_long, link_long = pattern.identify_links_upd_lanes(
                lanes_long, link_long, sw_type, ids, max_link_length, sw_length
            )
        # Save link_long for main because there where some minor fixes made to
        # it. Will be used for vehicle count interpolation. 
        # Otherwise keep updating directly to lanes_long and save out later.
        if sw_type == 'main':
            link_long.to_csv(pp_main_points_upd)
    
    # Save the long lanes dataframe.
    lanes_long.to_csv(pp_links_out)
    

        
        






