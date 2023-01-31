# -*- coding: utf-8 -*-
"""

Runtime: approx. 3 minutes. 
"""

import os
import pandas as pd
import numpy as np
from geojson import dump
import arcpy
arcpy.env.overwriteOutput = True


from arcgis import GIS

import src.process_numdata.functions.f04_concat_functions as concat
from project_paths import get_project_paths, get_params
from src.process_numdata.functions.f04_concat_functions import load_csv_and_fix_id


pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))

dates = sorted(param['dates'])
pp_fin_out = os.path.join(pp['fin_data'], 'df.csv')
pp_dict_out = os.path.join(pp['fin_data'], 'df_dict.json')

pp_fin_out_me = os.path.join(pp['data_temp'], 'df_me.csv')
pp_dict_out_me = os.path.join(pp['data_temp'], 'df_me_dict.json')


# Übersicht über Datensätze: 
# Name         # diff year # Description          # Required operations
#------------------------------------------------------------------------------
# BL           #  No       #                      # subset & get_dummies
# elevation    #  No       #                      # select
# zeb_inkar    #  No       # ZEB and INKAR data   # None
# sk_node_loc  #  No       # Sector,node&location # None 
# geometry     #  No       # Describes curvature  # None
# wind         #  No       # Wind over 1981-2000  # None
# -----------------------------------------------------------------------------
# weather_     #  Yes      #                      # None
# unf_df_      #  Yes      # info per accident    # aggregate
# traffic      #  Yes      #                      # select & get_dummies
# ot_final     #  Yes      # converted to shares  # select & get_dummies
# -----------------------------------------------------------------------------
type_dict = {'year_specific': ['weather', 'crash', 'traffic', 'ot'],
             'constant':      ['BL', 'elevation', 'zeb_inkar', 'sk_node',
                               'geometry','wind']}

# Make complete data frames for the individual years. 
for date in dates:
    
    year = date[:4]
    # Create dictionary of all required paths. 
    p_df = pp['data_out_df']
    path_dict = {'BL':         os.path.join(p_df, 'BL.csv'),
                 'elevation':  os.path.join(p_df, 'elevation.csv'),
                 'zeb_inkar':  os.path.join(p_df, 'zeb_inkar.csv'),
                 'sk_node_loc':os.path.join(p_df, 'sk_node_loc.csv'),
                 'geometry':   os.path.join(p_df, 'geometry.csv'),
                 'wind':       os.path.join(p_df, 'wind.csv'),
                 # Year specific: 
                 'weather':    os.path.join(p_df, 'weather_'+year+'.csv'),
                 'crash':      os.path.join(p_df, 'unf_df_'+date+'.csv'),
                 'traffic':    os.path.join(p_df, 'counts_interpolated'+year+'.csv'),
                 'ot':         os.path.join(p_df,'ot_final'+date+'.csv'),
                 # Out (year specific):
                 'out':        os.path.join(pp['fin_data'], 'df'+date+'.csv'),
                 'out_me':     os.path.join(pp['fin_data'], 'df'+date+'_me.csv')
            }
    
    # Concatenate everything to a year specific dataframe and save out. 
    df, colname_dict = concat.concat_year_data(path_dict, year, param, False)   
    df_me, colname_dict_me = concat.concat_year_data(path_dict, year, 
                                                     param, True)    

    print(date)
    print(df.shape)
    df.to_csv(path_dict['out'])
    df_me.to_csv(path_dict['out_me'])

#####
# Also make a long data frame with all years. 
df = pd.DataFrame()
df_me = pd.DataFrame()
for date in dates:
    path = os.path.join(pp['fin_data'], 'df'+date+'.csv')
    path_me = os.path.join(pp['fin_data'], 'df'+date+'_me.csv')
    
    df_n = load_csv_and_fix_id(path)
    df_n['year'] = int(date[:4])
    df = pd.concat([df, df_n])

    df_m = load_csv_and_fix_id(path_me)
    df_m['year'] = int(date[:4])
    df_me = pd.concat([df_me, df_m])

# Drop rows if these states are 1, since they have only data for 1 or 2 years!
df.update(df[['Berlin', 'NRW', 'Thüringen']].replace(np.nan,0))
df = concat.deselect_rows(df, {'':['Berlin', 'NRW', 'Thüringen']}) 

df.to_csv(os.path.join(pp['fin_data'], 'df_long.csv'))


# Now, construct a wide df and aggregate over years. 
# Not specified and thus treated as constant over the years (value from 2019):
#   BL, elevation, zeb_inkar, and wind
var_cols_sum     =colname_dict['crash']
var_cols_average =colname_dict['weather']+colname_dict['traffic']
same_cols = [x for x in colname_dict['ot'] if (
                    not x in ['maxspeedchange', 'large_sw','new_var']
                ) & (not x.startswith('lane'))]
df = concat.aggregate_years(df, year_list=[x[:4] for x in dates], 
                            same_cols=same_cols, var_cols_sum=var_cols_sum, 
                            var_cols_average =var_cols_average,
                            change_tolerance=0.2)
df_me = concat.aggregate_years(df_me, [x[:4] for x in dates], same_cols, 
                               var_cols_sum, var_cols_average, 0.2, 
                               max_extent=True)

# Drop all rows with missing values...
df.dropna(axis=0, inplace=True)    

##
# Add spatial lags 
lag_cols = ['lane_change', 'maxspeedchange', 'main_Entry', 
            'sec_Entry', 'main_Exit', 'sec_Exit', 'tunnel', 'maxspeed_none',
            'maxspeed_100', 'maxspeed_120', 'maxspeed_130', 'bridge_total', 
            'maxspeed_vari', 'down_change', 
            'up_change', 'net_change', 'max_slope', 'mean_slope', 'elev_change',
            'surf_het', 'sub_mean', 'perf_mean', 'asphalt', 
            'node_area', 'right_turn', 'left_turn', 'total_turn', 'sharp_turns',
            'n_lanes']
# How lags are imputed for starting segments. 
lag0 = {1: ['lane_change', 'maxspeedchange', 'main_Entry', 'node_area', 
            'sharp_turns', 'maxspeed_vari'],
        0: ['sec_Entry', 'main_Exit', 'sec_Exit', 'tunnel', 'maxspeed_none',
            'maxspeed_100', 'maxspeed_120', 'maxspeed_130',
            'bridge_total', 'left_turn'],
        30: ['right_turn', 'total_turn'],
        'col': []}
lag0.update({'col':[x for x in lag_cols if x not in lag0[1]+lag0[0]+lag0[30]]})

# Make a dataframe that contains all rows in the frame, even with missings. 
# df_me = concat.concat_year_data(path_dict, year, param, True)[0]
df_lag = df_me[lag_cols] 
df = concat.add_lags(df, df_lag, lag_cols, lag0)
#####
df_me = concat.add_lags(df_me, df_lag, lag_cols, lag0)

# Make a geodataframe of the max_extent data in 2019. 
pp_frame = os.path.join(pp['data_out_km'], 'frame.shp')
pp_spatial_out = os.path.join(pp['fin_data'], 'data.shp')
frame_df = pd.DataFrame.spatial.from_featureclass(pp_frame)
# Join maximum extent data to the  frame. 
frame_df.set_index('ID', drop=False, inplace=True)
frame_df = frame_df.join(df_me)

sl_cols = [x for x in frame_df.columns if ('maxspeed' in x) & 
           (x not in ['maxspeed_cond','maxspeedchange'])]
# Create one variable for maxspeed for visualization.
frame_df['maxspeed'] = pd.DataFrame(frame_df[sl_cols].idxmax(axis=1)).apply(
                                        lambda x: '_'.join(
                                                str(x[0]).split('_')[1:]),
                                        axis=1)


frame_df = frame_df[[x for x in frame_df.columns if x != 'SHAPE']+['SHAPE']]
ms_cols = [x for x in frame_df.columns if 'maxspeed' in x]
frame_df.rename(columns=dict(zip(ms_cols, 
                                 [x.replace('maxspeed','ms') for x in ms_cols])),
                inplace=True)
frame_df.spatial.to_featureclass(r"in_memory\data")
arcpy.management.CopyFeatures(r"in_memory\data", pp_spatial_out)


# Save out max extent data
df_me.to_csv(os.path.join(pp['fin_data'], 'df_max_extent.csv'))


# Save out final wide, aggregated dataframe. 
df.to_csv(pp_fin_out)

with open(pp_dict_out, 'w') as f:
    dump(colname_dict, f)



