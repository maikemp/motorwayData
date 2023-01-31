# -*- coding: utf-8 -*-
"""

Runtime: approx. 15 minutes.
""" 


import os 
import json

import src.process_numdata.functions.f02_collapse_functions as collapse 
from src.process_numdata.functions.f04_concat_functions import load_csv_and_fix_id
from src.process_numdata.functions.f01_pattern_functions import (
        make_id_sequence_dict, read_data_select_cols_and_relabel)

from project_paths import get_project_paths, get_params

# {Variable name root: {area: length of change, new_var: name of new column}}. 
pattern_dict = {'bridge': {'area':(0 , 40 ), 'new_var':'small_bridge'}}
    
pp = get_project_paths()

param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))
dates = param['dates']
id_length = param['id_length']
in_seg_id_length = param['in_seg_id_length']
sw_length = param['std_slipway_length']
max_link_length = param['max_link_length']

patt_param = json.load(open(os.path.join(pp['lib'], 'pattern_ident_lib.json')))
len_tol = patt_param['road_length_tolerance']


date = param['reference_date']


    
# Construct new variables
# for date in dates...
# Second: Collapse lanes_long to shares. Do the same for other variables.
for date in dates:
    
    print('Start with date:', date)
    pp_lanes_links_long = os.path.join(pp['data_out_df'], 
                                       'links_on_lanes_' + date[:4] + '.csv') 
    pp_in_ot = os.path.join(pp['data_out_df'], 'ot_df_'+date+'.csv')
    pp_out_ot = os.path.join(pp['data_out_df'], 'ot_final'+date+'.csv')
    
    df_ot = read_data_select_cols_and_relabel(
        pp_in_ot, param['id_length'], tuple(patt_param['remove_col_starts']), 
        tuple(patt_param['relabel_col_starts']), patt_param['relabel_dict']
    )
    df_ot.index_name = 'id'
    # Make some lists of columns.
    rep_col, start_col, collapse_col, unique_col = collapse.get_collists(df_ot)
    # This is the df_ot that all the other values will be matched to.
    df_master = df_ot[[x for x in unique_col if x !='length']]
    
    # Create a dictionary of subsequent indices.
    id_dict = make_id_sequence_dict(df_ot, id_length, in_seg_id_length)
    
    lanes_long = load_csv_and_fix_id(pp_lanes_links_long)
    lanes_long.index.name = 'id'
    # Match check_length to the df. 
    lanes_long = lanes_long.join(df_ot['check_length'])
    
    lanes_long.set_index('num', append=True, inplace=True)
    lanes_long.set_index(lanes_long.index.rename(('id', 'num')), inplace=True)
    lanes_long.rename({'id.1': 'id', 'num.1': 'num'}, inplace=True, axis=1)
    lanes_long, sw_vars = collapse.add_differentiated_sw_type_vars(lanes_long)

    
    # Collapse the dataframe and reduce to shares of each value.
    # Concatenate to df_master
    new_collapse_col = ['lanes', 'Slipway_lanes'] + sw_vars
    for var in new_collapse_col:
        print(var)
        df_sel = lanes_long[[x for x in lanes_long.columns if x in 
                             [var, 'length', 'check_length']]]
        unique_vals = collapse.get_unique_vals(df_sel, var)
        df_sel = collapse.values_to_columns(df_sel, var, unique_vals, 
                                   delete_old_cols=True)
        df_sel = collapse.to_wide_and_sum_over_values(df_sel)
        df_master = df_master.join(df_sel)
    
    for var in sw_vars:
        df_master.drop(columns=var+':0', inplace=True)
        df_master.rename(columns={var+':1': var}, inplace=True)
    
    df_master.drop(columns=['Slipway_lanes:0', 'Slipway_lanes:1'], inplace=True)

    # Next for all other collapse variables. 
    collapse_col = [x for x in collapse_col if x not in new_collapse_col]
    for var in collapse_col:
        print(var)
        # Here, identify other patterns.
        if var in pattern_dict.keys():
            
            new_var = pattern_dict[var]['new_var']
            area = pattern_dict[var]['area']
            
            # Collapse df to its informative aspects with respect to one variable.
            df_col = collapse.select_and_collapse_columns(df_ot, var, rename=True)
            df_col = df_col.join(df_ot['check_length'])
            df_col = collapse.change_df_by_pattern(df_col, var, area, new_var, 
                                                  id_dict, change_var=True, 
                                                  return_form='long')
            
            for v in [var, new_var]:
                
                df_sel = df_col[[v, 'length', 'check_length']]
                
                unique_vals = collapse.get_unique_vals(df_sel, v)
                df_sel = collapse.values_to_columns(df_sel, v, unique_vals, 
                                           delete_old_cols=True)
                df_sel = collapse.to_wide_and_sum_over_values(df_sel)
                df_master = df_master.join(df_sel)
        else:
    
            df_sel = collapse.select_and_wide_to_long(df_ot, var, rename=True)
         
            if var == 'maxspeed_c':
                mc_dict = collapse.get_maxspeed_c_dict()
                df_sel.replace(mc_dict, inplace=True)
            
            unique_vals = collapse.get_unique_vals(df_sel, var)
            df_sel = collapse.values_to_columns(df_sel, var, unique_vals, 
                                       delete_old_cols=True)
            df_sel = collapse.to_wide_and_sum_over_values(df_sel)
            df_master = df_master.join(df_sel)
    
    df_master.to_csv(pp_out_ot)



            



        
        






