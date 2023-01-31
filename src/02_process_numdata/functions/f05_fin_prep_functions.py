# -*- coding: utf-8 -*-
"""
Functions for final data preparation. 

"""

import pandas as pd
import numpy as np 
import random 

from sklearn import preprocessing
from sklearn import decomposition
from sklearn.model_selection import GroupShuffleSplit




def group_train_test_split(df, rf_cols, outcome, location, test_size):
    
    splitter = GroupShuffleSplit(test_size=test_size, n_splits=1)
    split = splitter.split(df, groups=df[location])
    train_inds, test_inds = next(split)
    
    train = df.iloc[train_inds]
    test = df.iloc[test_inds]
    X_train = train[rf_cols]
    X_test = test[rf_cols]
    y_train = train[outcome]
    y_test = test[outcome]
    loc_train = train[location]
    loc_test = test[location]
    
    return(X_train, X_test, y_train, y_test, loc_train, loc_test)


def get_holdout(df_sub, group_var, limit, x_cols, share=0.2, 
                tolerance=50, sort_cols=['n_t','n'], verbose=1):
    
    k = round(100/(100*share))
    print('Actual share of hold out set:', 100/k, '%')
    
    df_aux = get_aux_df(df_sub, group_var=group_var, out_var=limit) 
    folds = draw_random_stratified_location_folds(
                df_aux, k=k, tolerance=tolerance, sort_cols=sort_cols) 
    # Set aside validation set. 
    X_train, X_val, y_train, y_val, loc_train, loc_test = execute_split(
            df_sub, x_cols, limit, folds, i=0, verbose=1)
    
    return(X_train, X_val, y_train, y_val, loc_train, loc_test)


def get_aux_df(df_sub, group_var, out_var):
    
    df_aux = pd.DataFrame(df_sub[group_var].value_counts()
                     ).sort_index().rename(columns={group_var: 'n'})
    df_aux['s_t'] = df_sub.groupby(group_var)[out_var].mean() 
    df_aux['n_t'] = df_sub.groupby(group_var)[out_var].sum() 
    df_aux['n_ut'] = df_aux['n'] - df_aux['n_t']
    
    return(df_aux)


def execute_split(df_sub, x_cols, limit, folds, i, verbose=1):
    
    
    X_train = df_sub[x_cols].loc[
            df_sub['location'].isin(folds[i]['ids'])==False]
    y_train = df_sub[limit].loc[
            df_sub['location'].isin(folds[i]['ids'])==False]
    loc_train = df_sub['location'].loc[
            df_sub['location'].isin(folds[i]['ids'])==False]    
    X_test = df_sub[x_cols].loc[
            df_sub['location'].isin(folds[i]['ids'])]
    y_test = df_sub[limit].loc[
            df_sub['location'].isin(folds[i]['ids'])]
    loc_test = df_sub['location'].loc[
            df_sub['location'].isin(folds[i]['ids'])]
    
    if verbose == 1:
        print('Testset share:',round(len(y_test)/(len(y_train)+len(y_test)),2),
              '(Ideal:{})'.format(round(1/len(folds.keys()),2)),
              'Treated ratio:', round(np.mean(y_test)/np.mean(y_train),2),
              '(Ideal:1)')
    
    return(X_train, X_test, y_train, y_test, loc_train, loc_test)


def draw_random_stratified_location_folds(df_aux, k, tolerance=0, 
                                          sort_cols=['n_t','n'], verbose=1):
    '''Returns a dictionary with k id_lists generated from df_aux. 
    It attempts to reach a roughly equal number of treated and untreated 
    observations in all folds.
    
    
    '''
    
    av_s_size = df_aux.n.sum()/k
    av_s_treated = sum(df_aux['s_t']*df_aux['n'])/df_aux.n.sum()
    av_n_t = av_s_size*av_s_treated
    av_n_ut = av_s_size-av_n_t
    
    track = dict(zip(list(range(k)), 
                     [{'ids':[], 'n':0, 'n_t':0, 'n_ut':0,
                       'open':{'t':av_n_t , 'ut':av_n_ut}} for i in range(k)]))
    
    
    # First, divide largest k*3 locations randomly and track size & n_treated.
    df_aux.sort_values(sort_cols, ascending=False, inplace=True)
    for i in list(df_aux.index):
        r = df_aux.loc[i]
        
        # if treated is 0, only untreated matters, and vice versa. 
        # ...... 
        if r.n_t == 0:
            a = [x for x in list(range(k)) if track[x]['open']['ut'] > r.n_ut]
        elif r.n_ut == 0:
            a = [x for x in list(range(k)) if track[x]['open']['t'] > r.n_t]
        else:
            a = [x for x in list(range(k)) if (track[x]['open']['t'] > r.n_t) &
                                              (track[x]['open']['ut']> r.n_ut)]
        # If it was not yet possible to assign value, increase tolerance.
        t = 1
        while (len(a) == 0) & (t < tolerance):
            a =[x for x in list(range(k)) if (track[x]['open']['t']+t>r.n_t) & 
                                             (track[x]['open']['ut']+t>r.n_ut)]
            t += 1
        if len(a) == 0:
            raise ValueError('''Not possible to assign all observations to fold.
                             Try increasing tolerance, decreasing k
                             or choosing smaller locations''')
        r_pick = random.choice(a)
    
        # Now, update the dictionary accoringly.
        track[r_pick]['ids'] += [i]
        track[r_pick]['n'] += r.n
        track[r_pick]['n_t'] += r.n_t
        track[r_pick]['n_ut'] += r.n_ut
        track[r_pick]['open']['t'] -= r.n_t
        track[r_pick]['open']['ut'] -= r.n_ut
    
    # Add share of treated. 
    for i in range(k):
        track[i]['s_t'] = track[i]['n_t']/track[i]['n']
    
    o_t = []
    o_ut = []
    for i in range(k):
        o_t.append(track[i]['open']['t'])
        o_ut.append(track[i]['open']['ut'])
        
    if verbose == 1:
        print('Standard deviation of number of treated:', np.std(o_t))
        print('Standard deviation of number of untreated:', np.std(o_ut))
    
    return(track)


def get_rf_cols(x_cols, cd):
    rf_cols = [x for x in x_cols if x not in cd['BL'] and not 
               x.startswith('lane') and not x in ['bridge', 'small_bridge', 
                           'sharp_turns', 'large_sw', 'perf_num', 'sub_num']]+\
              ['node_area_lag', 'up_change_lag', 'down_change_lag', 
               'total_turn_lag', 'n_lanes_lag']
    return(rf_cols)


def get_collists(df, cd, lags=False):
    '''If lags=True, lags are included in x_cols. Default is False.'''
    
    crash_cols = ['fatal', 'severely_injured', 'lightly_injured', 'total']
    for i in ['total_rate', 'light_rate', 'fatal_rate', 'fatal_severe_rate', 
              'severe_rate']:
        if i in df.columns:
            crash_cols = crash_cols + [i]
    sl_cols = ['maxspeed_100','maxspeed_120','maxspeed_130','maxspeed_none',
               'maxspeed_vari'] 
    
    # Instead use pca. 
    #exclude_cols = ['precipGE10mm_day', 'precipGE20mm_day', 'snowcover_days', 
    #                'ice_days', 'summer_duration']
    x_cols = [x for x in df.columns if (not x.startswith('maxspeed')) &
                  (not x in crash_cols+sl_cols+['location', 'Sector', 
                                                'random'])]
    if lags == False:
        x_cols = [x for x in x_cols if not x.endswith('_lag')] 
    # For spatial variables, take principal components rahter than all variables. 
    spatial_cols = cd['wind'] + cd['weather'] + \
        [x for x in cd['zeb_inkar'] if not x.startswith(
            ('sub','perf','surf','asphalt', 'pop_dens','random')
         )] + ['elevation', 'mean_slope']
    spatial_cols = [x for x in spatial_cols if x in x_cols]
    
    return(crash_cols, sl_cols, x_cols, spatial_cols)


def replace_by_pca(df, spatial_cols, x_cols, cd, explained_var=0.9, verbose=1,
                   pca=None):
    # Only include principal components of the spatial variables? 
    spatial_X =  df[spatial_cols]
    # Normalize columns.
    scaler_mm = preprocessing.MinMaxScaler()
    spatial_X[spatial_cols] = scaler_mm.fit_transform(spatial_X)
    # If no pca is passed, fit to data. if it is passed, use that one.
    if pca == None:
        pca = decomposition.PCA(n_components=explained_var)
        pca.fit(spatial_X)
        
    pca_cols = ['p_comp_'+str(x) for x in range(pca.n_components_)]
    x_pca = pd.DataFrame(pca.transform(df[spatial_cols]), 
                         columns=pca_cols,
                         index=spatial_X.index)
    df = df.join(x_pca).drop(columns=spatial_cols)
    
    x_cols_new = [x for x in x_cols if x not in spatial_cols] + pca_cols
    
    
    if verbose==1:
        print(len(pca_cols), 'principal components replaced', len(spatial_cols), 
              'variables. Reduction of',len(x_cols)-len(x_cols_new),'columns.',
              'Explained variance:',round(pca.explained_variance_ratio_.sum(),3))
    cd['pca'] = pca_cols
    
    return(df, x_cols_new, cd, pca)  
    
    