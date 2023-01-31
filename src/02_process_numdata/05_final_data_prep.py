# -*- coding: utf-8 -*-
"""
Do all final dataset manipulations.
Runtime: <1 Minute.

"""

import pandas as pd
import numpy as np
import os 
import matplotlib.pyplot as plt
import json
import random

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from project_paths import get_project_paths, get_params
from src.process_numdata.functions.f04_concat_functions import load_csv_and_fix_id
from src.data_analysis.functions.f00_CVsplitter import BalancedStratifiedGroupKFold
import src.process_numdata.functions.f05_fin_prep_functions as fp
import src.data_analysis.functions.f08_visualize_functions as visual

cate_vars=['AADT', 'HT_share', 'elev_change', 'node_area', 'total_turn', 
           'down_change','ramps','holiday']

random.seed(845) 

pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))
predict_param = get_params(path=os.path.join(pp['lib'], 'prediction_params.json'))

# Define parameters
reps = predict_param["reps"] 
p_comps = 4

# Define paths.
in_path = os.path.join(pp['fin_data'], 'df.csv')
pp_col_dict = os.path.join(pp['fin_data'], 'df_dict.json')
out_path = os.path.join(pp['fin_data'], 'df_prep.csv')

# Load stuff.
df = load_csv_and_fix_id(in_path)
cd = json.load(open(pp_col_dict, encoding="utf8"))

# Define lists of 
crash_cols, sl_cols, x_cols, spatial_cols = fp.get_collists(df, cd, lags=False)

# Add crash rates.
df['total_rate'] = df['total']/df['AADT']*10000
df['fatal_rate'] = df['fatal']/df['AADT']*10000
df['severe_rate'] = df['severely_injured']/df['AADT']*10000
df['light_rate'] = df['lightly_injured']/df['AADT']*10000
df['fatal_severe_rate'] = (df['fatal']+df['severely_injured'])/df['AADT']*10000

rate_cols = ['total_rate','fatal_rate','severe_rate','light_rate',
             'fatal_severe_rate']

# PCA to reduce dimension of spatial columns.    
df, x_cols, cd, pca = fp.replace_by_pca(df, spatial_cols, x_cols, cd, 
                                        explained_var=p_comps)

# Set speed limit variable to the prevalent speed limit (dummy). 
df[sl_cols] = pd.get_dummies(df[sl_cols].idxmax(axis=1))

# Derive new variables. Reason: more powerful in FFS:
    # Some variables may only be useful in combination, especially when looking 
    # at changes. Therefore, it may be more relevant to have more informative 
    # features. 
#df['ms_change_area'] = df['maxspeedchange_lag'] + df['maxspeedchange']
df['maxspeed_change'] = df['maxspeedchange']
df['ramps'] = df['main_Entry'] + df['main_Exit'] + df['sec_Entry'] + df['sec_Exit'] 
df['ramps_lag'] = df['main_Entry_lag'] + df['main_Exit_lag'] + \
               df['sec_Entry_lag'] + df['sec_Exit_lag'] 

# Remove no shoulder info on slipways. 
df.loc[(df['no_shoulder'] == df['ramps'])&(df['no_shoulder']>0),'no_shoulder'] = 0

# Delete unusable data.
# All segments with variable speed limit and the respective column. 
df = df[df['maxspeed_vari']==0].drop(columns=['maxspeed_vari'])
sl_cols = [x for x in sl_cols if x != 'maxspeed_vari']
# Columns with unusable info, like various dummy variables etc. 
df = df.drop(columns=[x for x in cd['BL'] if x in df.columns] + crash_cols +
                     [x for x in df.columns if 
                      x.startswith(('lanes_','maxspeedchange', 'lane_change'))] + 
                     ['large_sw', 'sub_num', 'perf_num', 'net_change', 
                      'Sector', 'sharp_turns', 
                      'net_change_lag', 'sharp_turns_lag'] +
                     [x+'_lag' for x in sl_cols])

# Keep only information with total bridge length, not divided by small and big.
df = df.drop(columns=['bridge', 'small_bridge']).rename(columns={
                     'bridge_total':'bridge','bridge_total_lag':'bridge_lag'})
df = df.rename(columns={'random':'aux_random_spatial'})

# Add random variables. 
for i in range(reps):
    df['aux_random_'+str(i)] = np.random.uniform(0, 1, len(df))

# Add flags for holdout set 
# Start with W.hat, i.e. for each sub-dataframe. 
for limit in sl_cols[:3]:
    for i in range(reps):
        try:
            df_sub = df[(df[limit]==1)|df['maxspeed_none']==1]
            Xt, Xv, yt, yv, loc_train, loc_val = fp.get_holdout(
                df_sub, 'location', limit, ['bridge'], share=0.2, tolerance=50
            )
            var_str = 'aux_holdout_'+limit[-3:]+'_'+str(i)
            df[var_str] = np.nan
            df[var_str].loc[loc_val.index] = 1
            df[var_str].loc[loc_train.index] = 0
        except:
            df_sub = df[(df[limit]==1)|df['maxspeed_none']==1]
            X_train, X_val, y_train, y_val, loc_train, loc_val = fp.get_holdout(
                df_sub, 'location', limit, ['bridge'], share=0.2, tolerance=220
            )
            var_str = 'aux_holdout_'+limit[-3:]+'_'+str(i)
            df[var_str] = np.nan
            df[var_str].loc[loc_val.index] = 1
            df[var_str].loc[loc_train.index] = 0

# Next, continue with general hodlout set, for the outcome estimation. 
            # Balance speed limit distribution?? 
for i in range(reps):
    Xt, Xv, yt, yv, loc_train, loc_val = fp.group_train_test_split(
                                df, ['bridge'], 'total_rate', 'location', 0.2)
    var_str = 'aux_holdout_'+str(i)
    df[var_str] = np.nan
    df[var_str].loc[loc_val.index] = 1
    df[var_str].loc[loc_train.index] = 0

# Save out for further use. 
df.to_csv(out_path)


# Get cutoff values for CATEs Used for GRF:  
# Selects threshold such that both blocks contains 50% of traffic
cates_d = visual.get_cates_cutpoints_dict(df, cate_vars)

# Save out cutoff values to library.
with open(os.path.join(pp['lib'],'cates_cutoffs.json'), 'w') as f:
    json.dump(cates_d, f)


## Also make a df with the spatial variables instead of the principal components.
df_s = load_csv_and_fix_id(in_path)
df = load_csv_and_fix_id(out_path)
cols = spatial_cols + [x for x in cd['BL'] if x in df_s.columns]
df = df.join(df_s[cols]).drop(columns=[x for x in df.columns if 
                                               x.startswith('p_comp_')])

df.to_csv(os.path.join(pp['fin_data'], 'df_prep_spatial.csv'))





###############  ###################  ######################
# Make a full extent data without any missings in x_cols.
max_ext_path = os.path.join(pp['fin_data'], 'df_max_extent.csv')
df_me = load_csv_and_fix_id(max_ext_path)

df_BL = load_csv_and_fix_id(os.path.join(pp['data_out_df'], 'BL.csv'))
df_me = df_me.join(df_BL)

# Make a complete df, without missings, even if imputations are not ideal.
miss_cols = df_me.isna().sum()[df_me.isna().sum()>0].keys()
for mc in miss_cols:
    # First try sector mean.
    df_me[mc] = df_me[mc].fillna(df_me.groupby('Sector')[mc].transform('mean'))
    # Then try location mean.
    df_me[mc] = df_me[mc].fillna(df_me.groupby('location')[mc].transform('mean'))
    # Then try Kreis mean. (Because usually for tiny motorways)
    df_me[mc] = df_me[mc].fillna(df_me.groupby('GEN_KRS')[mc].transform('mean'))
    # Then simply take BL mean. (Should not even be required.)
    df_me[mc] = df_me[mc].fillna(df_me.groupby('GEN_LAN')[mc].transform('mean'))


# Replace spatial cols by pc, use pca that was fit above
df_me, x_cols, cd, pca = fp.replace_by_pca(df_me, spatial_cols, x_cols, cd, 
                                           explained_var=p_comps, pca=pca)

df_me['maxspeed_change'] = df_me['maxspeedchange']
df_me['ramps'] = df_me[['main_Entry','main_Exit','sec_Entry','sec_Exit']].sum(axis=1)
df_me['ramps_lag'] = df_me[['main_Entry_lag','main_Exit_lag',
                      'sec_Entry_lag','sec_Exit_lag']].sum(axis=1) 

# Remove no shoulder info on slipways. 
df_me.loc[(df_me['no_shoulder'] == df_me['ramps'])&(df_me['no_shoulder']>0),
          'no_shoulder'] = 0

df_me.drop(columns='bridge', inplace=True)
df_me.rename(columns={'bridge_total':'bridge','bridge_total_lag':'bridge_lag'},
             inplace=True)

df_me['ms_other'] = 1-df_me[sl_cols].sum(axis=1)
df_me[sl_cols+['ms_other']] = pd.get_dummies(df_me[sl_cols+['ms_other']].idxmax(axis=1))


df_me.to_csv(os.path.join(pp['fin_data'], 'df_max_ext_imputed.csv'))






