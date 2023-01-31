# -*- coding: utf-8 -*-
"""
Runtime: approx. 15 Minutes. 
"""


import os 

import src.process_numdata.functions.f03_count_functions as count

from src.process_numdata.functions.f04_concat_functions import load_csv_and_fix_id
from project_paths import get_project_paths, get_params
   
pp = get_project_paths()

param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))
dates = param['dates']
id_length = param['id_length']
in_seg_id_length = param['in_seg_id_length']


date = param['reference_date']


for date in dates: 
    
    print("Start with date:", date)
    pp_in_bisstra,pp_fs,pp_in_count,pp_out_count = count.get_count_paths(pp, date)
    
    
    # Load counting station data and add some further information on the segment. 
    df_c = load_csv_and_fix_id(pp_in_count, encoding="latin1")
    # Store index as it is and in a numerical veresion.
    df_c['num_id'] = [int(x) for x in list(df_c.index)]
    df_c['id'] = list(df_c.index)
    
    df_c = count.join_data(df_c, 
                {os.path.join(pp['data_out_df'], 'sk.csv'):    [],
                 os.path.join(pp['data_out_df'], 'frame.csv'): ['one_direct']})
    val_cols = [x for x in df_c.columns if x not in 
                    ['compass', 'direction', 'one_direct', 'AbschnittA', 'DZNr']] 
    zs_cols = [x for x in df_c.columns if x not in [
            'compass', 'direction', 'sk_kennung', 'Sector', 'one_direct_all']]
    fs_cols1 = ["DTV19", "DTVw19","DTVu19", "DTVs19", "DTVSV", "fer", "bSo",  
                "MSVI", "bSVI", "MSVII", "bSVII", "Mt", "Mn", "pt", "pn", "DZ"]
    
    # Replace sk_kennung of the frame if counting station has different identifier, 
    # if AbschnittA is not zero and the segments are clearly succeeding. 
    df_c = count.fix_sk_kennung_remove_false_stations(df_c, val_cols)
    
    val_cols_2 = [x for x in val_cols if x not in ['num_id', 'id']]
    df_c, new_val_dict = count.get_newval_dict_existing_stations(df_c, val_cols_2)
    for k, v in new_val_dict.items():
        df_c.loc[k, val_cols_2] = v
    
    df_c, df_mm = count.load_count_data_merge_and_get_minmax(df_c, pp_in_bisstra, 
                                                             pp_fs, fs_cols1)
    
    df_mm.drop(columns=['compass', 'direction', 'one_direct', 'AbschnittA', 
                        'AnzFsQ', 'pMn', 'R1', 'Nahziel', 'Hi', 'DL', 'bMo', 
                        'Sector', 'vnk', 'nnk', 'vnk_nnk', 'DTVw19', 'DTVu19', 
                        'DTVs19','pt', 'pn'], inplace=True)
    fs_cols = [x for x in [x+'_fs' if x+'_fs' in df_mm.columns else x for x in 
                           fs_cols1] if x in list(set(list(df_mm.columns))-{'DZ'})]
    ##############################################################################
    # Derive distance weights. 
    # First, make a list of all sectors and their distance from one another. 
    # Dictionaries have to be sorted by a to differentiate between directions.
    d = count.get_sec_length_dict(df_mm)
    
    weight_dict = count.get_weight_dict(df_mm, d, impute_var='TVKfzms')
    # If there is no counting station on a a-segment, just apply a general time/measure
                # scaling factor as the average over the whole motorway?
    
    # Only needed two segments to derive distances, now drop for simplicity. 
    df_mm.drop_duplicates(['sk_kennung','one_direct_all','a'], inplace=True)
    
    ###############################################################################
    # First, interpolate variables in the Fortschreibung, as these are often only 
    # single missing lines.
    df_mm = count.interpolate_values(df_mm, d, 
                                     primary_var='DTV19', 
                                     all_vars=fs_cols, 
                                     scale_vars=['DTV19', 'DTVSV', 'MSVI', 
                                                 'MSVII','Mt_fs', 'Mn_fs'], 
                                     scale_factor=0.8)
    
    
    # Before being able to impute things, need to derive the relevant relations. 
    cs_cols = ['fer','Mt','Mn','TVKfzms','TVSVms','SV5Kfzms','bSV50ms','bSo']
    
    relations = count.get_relations_dict(df_mm, impute_var='TVKfzms')
    
    var_mapping = {'fer':'fer_fs', 
                   'Mt': 'Mt_fs',
                   'Mn':'Mn_fs', 
                   'TVKfzms': 'DTV19',
                   'TVSVms': 'DTVSV',
                   'SV5Kfzms': ['MSVI','MSVII'],
                   'bSV50ms': ['bSVI','bSVII'],
                   'bSo': 'bSo_fs'}
        
    # Impute values by transfering relations observed at counting stations
    # weighted by distance measure to each segment for which it is poosible. 
    df_mm = count.impute_values_by_dist_and_relation(df_mm, weight_dict, relations)
    
    # To handle more missings, interpolate missings in the target columns, 
    # where no Fortschreibung is available. 
    df_mm = count.interpolate_values(df_mm, 
                                     d, 
                                     primary_var='TVKfzms', 
                                     all_vars=cs_cols, 
                                     scale_vars=['TVKfzms', 'TVSVms', 'SV5Kfzms', 
                                                 'Mt', 'Mn'], 
                                     scale_factor=0.8)
    
    # Last option to impute missings: if no DZ in a, simply transfer 
    # Fortschreibungs-values with some average relational factors. 
    
    df_mm = count.search_relation_and_transfer_fs(df_mm, var_mapping, relations, 
                                                  cs_cols)
    
                                
    # --> Now, only nans left should be those that are alone in a, have no counting 
        # stations AND have not Fortschreibungs values. 
    
    
    # Now, transfer info from df_mm to the whole dataframe. 
    df_c = count.transfer_info_to_df(df_c, df_mm, cs_cols)
    df_c = df_c[df_c['id'].notna()]
    # Select columns to keep.
    df_c = df_c[cs_cols+['compass', 'direction', 'id', 
                         'sk_kennung', 'Sector']].set_index('id')
    df_c.index.name = None
    # Save out data. 
    df_c.to_csv(pp_out_count)
    
    from datetime import datetime
    print(print(datetime.now().strftime("%H:%M:%S")), df_c.isna().sum())




