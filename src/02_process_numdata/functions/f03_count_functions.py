# -*- coding: utf-8 -*-
"""

"""

import numpy as np
import pandas as pd
import re
import os
from itertools import chain
from arcgis import GIS

from src.process_numdata.functions.f04_concat_functions import load_csv_and_fix_id


def transfer_info_to_df(df_c, df_mm, cs_cols):
    
    for a in df_c['a'].unique():
        for sk in df_c[df_c['a']==a]['sk_kennung'].unique():
            for r in [0,1]:
                f_c = (df_c['a']==a)&(df_c['one_direct_all']==r)&(df_c['sk_kennung']==sk)
                f_fs = (df_mm['a']==a)&(df_mm['one_direct_all']==r)&(df_mm['sk_kennung']==sk)
                if (sum(f_c)>0)&(sum(f_fs)>0):
                    df_c.loc[f_c, cs_cols] = df_mm.loc[df_mm.loc[f_fs].index.repeat(sum(f_c)),cs_cols].values
            
    return(df_c)


def search_relation_and_transfer_fs(df_mm, var_mapping, relations, cs_cols):
    df_mm['mw'] = df_mm['id'].apply(lambda x: x[:3])
        
    for var in cs_cols:
        # Get a list of as wich contain segments with missing values for this var. 
        aa = list(set([x for x in list(df_mm[df_mm[var].isna()].a)]))
        
        for a_i in aa:
            # Select all mws which are affected in this a. 
            i_s = [x for x in list(df_mm[df_mm['a']==a_i].id)]
            m = i_s[0][:3]
            # Get list of as on the same motorway.
            a_list = list(set(df_mm[df_mm['mw']==m].a))
            rel_a = dict((k, relations[k]) for k in a_list)
            # If not relations defined for this motorway, take similar motorways by
            # using only the first digit. 
            # Derive an average relations dictionary. 
            rel_new = get_new_relation(rel_a, var)
            
            if sum([len(x) for x in rel_a.values()]) > 0:
                rel_new = get_new_relation(rel_a, var)
                
            if (sum([len(x) for x in rel_a.values()]) == 0)|(pd.isna(rel_new)):
                a_list = list(set(df_mm[df_mm['id'].apply(lambda x: x[:2])==m[:2]].a))
                rel_a = dict((k, relations[k]) for k in a_list)
                rel_new = get_new_relation(rel_a, var)
                
            if (sum([len(x) for x in rel_a.values()]) == 0)|(pd.isna(rel_new)):
                a_list = list(set(df_mm[df_mm['id'].apply(lambda x: x[:1])==m[:1]].a))
                rel_a = dict((k, relations[k]) for k in a_list)
                rel_new = get_new_relation(rel_a, var)
    
            
            # Impute
            f = df_mm['id'].isin(i_s)
            if type(var_mapping[var])==str:
                df_mm.loc[f, var] = rel_new * df_mm.loc[f, var_mapping[var]]
            else:
                df_mm.loc[f, var] = rel_new * (df_mm.loc[f, var_mapping[var][0]]+df_mm.loc[f, var_mapping[var][1]])
                
    return(df_mm)


def impute_values_by_dist_and_relation(df_mm, weight_dict, relations):
    '''Impute by transfering relation observed at nearby counting stations,
        weighted with distance measures.
        May change to also pass variable mapping instead of this copying...'''
        
    for k, v in weight_dict.items():
        for s, w in v.items():
            f = (df_mm['a']==k)&(df_mm['sk_kennung']==s)
            if (len(w) == 1)&(sum(f)>0):
                r1 = relations[k][list(w.keys())[0]]
                x = df_mm[f].iloc[0]
                df_mm.loc[f,'fer']     = r1['fer']*x['fer_fs'] 
                df_mm.loc[f,'Mt']      = r1['Mt']*x['Mt_fs'] 
                df_mm.loc[f,'Mn']      = r1['Mn']*x['Mn_fs'] 
                df_mm.loc[f,'TVKfzms'] = r1['TVKfzms']*x['DTV19'] 
                df_mm.loc[f,'TVSVms']  = r1['TVSVms']*x['DTVSV'] 
                df_mm.loc[f,'SV5Kfzms']= r1['SV5Kfzms']*(x['MSVI']+x['MSVII'])
                df_mm.loc[f,'bSV50ms'] = r1['bSV50ms']*(x['bSVI']+x['bSVII'])
                df_mm.loc[f,'bSo']     = r1['bSo']*x['bSo_fs'] 
            elif (len(w) == 2)&(sum(f)>0): 
                r1 = relations[k][list(w.keys())[0]]
                w1 = list(w.values())[0]
                r2 = relations[k][list(w.keys())[1]]
                w2 = list(w.values())[1]
                x = df_mm[f].iloc[0]
    
                df_mm.loc[f,'fer']     = (w1*r1['fer'] + w2*r2['fer'])*x['fer_fs'] 
                df_mm.loc[f,'Mt']      = (w1*r1['Mt'] + w2*r2['Mt'])*x['Mt_fs'] 
                df_mm.loc[f,'Mn']      = (w1*r1['Mn'] + w2*r2['Mn'])*x['Mn_fs'] 
                df_mm.loc[f,'TVKfzms'] = (w1*r1['TVKfzms'] + w2*r2['TVKfzms'])*x['DTV19'] 
                df_mm.loc[f,'TVSVms']  = (w1*r1['TVSVms'] + w2*r2['TVSVms'])*x['DTVSV'] 
                df_mm.loc[f,'SV5Kfzms']= (w1*r1['SV5Kfzms'] + w2*r2['SV5Kfzms'])*(x['MSVI']+x['MSVII'])
                df_mm.loc[f,'bSV50ms'] = (w1*r1['bSV50ms'] + w2*r2['bSV50ms'])  *(x['bSVI']+x['bSVII'])
                df_mm.loc[f,'bSo']     = (w1*r1['bSo'] + w2*r2['bSo'])*x['bSo_fs'] 
                
    df_mm.sort_values('num_id', inplace=True)
    return(df_mm)


def get_new_relation(rel_a, var):
    rel_new = np.mean([
                       x for x in list(chain(*[[rel_a[x][i][var] for i in 
                                                list(rel_a[x].keys())]  
                                                for x in list(rel_a.keys())]))
                           if pd.isna(x) == False])
    return(rel_new)

def interpolate_values(df_mm,d,primary_var,all_vars,scale_vars,scale_factor):
    
    w = get_weight_dict(df_mm, d, impute_var=primary_var)
    # Scaling factor to scale down values that have only one weight.
    for k, v in w.items():
        for s, w in v.items():
            f = (df_mm['a']==k)&(df_mm['sk_kennung']==s)
            if (len(w) == 1)&(sum(f)>0):
                f1 = (df_mm['a']==k)&(df_mm['sk_kennung']==list(w.keys())[0])
                x1 = df_mm[f1].iloc[0]
                
                for var in all_vars:
                    if var in scale_vars:
                        scale = scale_factor
                    else:
                        scale = 1
                    
                    df_mm.loc[f,var]   = scale * x1[var]
                
            elif (len(w) == 2)&(sum(f)>0): 
                f1 = (df_mm['a']==k)&(df_mm['sk_kennung']==list(w.keys())[0])
                x1 = df_mm[f1].iloc[0]
                w1 = list(w.values())[0]
                f2 = (df_mm['a']==k)&(df_mm['sk_kennung']==list(w.keys())[1])
                x2 = df_mm[f2].iloc[0]
                w2 = list(w.values())[1]
                
                for var in all_vars:
                    df_mm.loc[f,var]   = w1 * x1[var]  + w2 * x2[var]
            
   # a=df_mm[['sk_kennung','a', 'id']+all_vars]
    
    # For the variables that still have more missings than the primary variable,
    # repeat imputation.
#    for var in list(a.isna().sum()[a.isna().sum()>a.isna().sum()[primary_var]].index):
    for var in [x for x in all_vars if x!=primary_var]:     # safer.
        
        w = get_weight_dict(df_mm, d, impute_var=var)
        for k, v in w.items():
            for s, w in v.items():
                f = (df_mm['a']==k)&(df_mm['sk_kennung']==s)
                if (len(w) == 1)&(sum(f)>0):
                    f1 = (df_mm['a']==k)&(df_mm['sk_kennung']==list(w.keys())[0])
                    x1 = df_mm[f1].iloc[0]
                    if var in scale_vars:
                        scale = scale_factor
                    else:
                        scale = 1
                    
                    df_mm.loc[f,var]   = scale * x1[var]
    
                    
                elif (len(w) == 2)&(sum(f)>0): 
                    f1 = (df_mm['a']==k)&(df_mm['sk_kennung']==list(w.keys())[0])
                    x1 = df_mm[f1].iloc[0]
                    w1 = list(w.values())[0]
                    f2 = (df_mm['a']==k)&(df_mm['sk_kennung']==list(w.keys())[1])
                    x2 = df_mm[f2].iloc[0]
                    w2 = list(w.values())[1]
                    
                    df_mm.loc[f,var]   = w1 * x1[var]  + w2 * x2[var]
    return(df_mm)


def get_relations_dict(df_mm, impute_var='TVKfzms'):
    relations = {}
    
    # Get all relations for sections with existing counting station.
    for i in df_mm['a'].unique():
        dt = df_mm[(df_mm[impute_var].isna()==False)&(df_mm['a']==i)]
        rel_dict = {}
        
        for cs in list(dt.sk_kennung):
            x = dt[dt['sk_kennung']==cs].iloc[0]
            rel_dict.update({cs: {'fer':     x['fer']/x['fer_fs'], 
                                  'Mt':      x['Mt']/x['Mt_fs'],
                                  'Mn':      x['Mn']/x['Mn_fs'],
                                  'TVKfzms': x['TVKfzms']/x['DTV19'],
                                  'TVSVms':  x['TVSVms']/x['DTVSV'],
                                  'SV5Kfzms':x['SV5Kfzms']/(x['MSVI']+x['MSVII']),
                                  'bSV50ms': x['bSV50ms']/(x['bSVI']+x['bSVII']),
                                  'bSo':      x['bSo']/x['bSo_fs']}})
        
        relations.update({i: rel_dict.copy()})
    return(relations)


def get_weight_dict(df_mm, d, impute_var='TVKfzms'):
    
    weight_dict = {}
    no_cs = []
    # Go through all a values. 
    for i in df_mm['a'].unique():
        temp = {}
        # Take only the sk_kennung values in i. 
        df_temp_i = df_mm[(df_mm['a']==i)].drop_duplicates('sk_kennung')
        # Merge the information about number of segments in sector to it. 
        e = pd.merge(df_temp_i,
                     pd.DataFrame.from_dict(d[i],
                                            orient="index").reset_index().rename(
                                                   columns={'index':'sk_kennung'}),
                     on='sk_kennung', how='left')
        # Take a dictionary of all counting sectors with counting stations. 
        cs_i = dict(zip(list(e[e[impute_var].isna()==False].sk_kennung),
                        list(e[e[impute_var].isna()==False].num_id)))
        # For all sectors without counting station, find at most two closest 
        # counting stations. 
        for j in [x for x in list(e['num_id']) if x not in list(cs_i.values())]:
            # If there is only one counting station just take it. 
            # No need to do distance weighting, as there is only this one station. 
            # Same holds if segment does not lie between two stations.
            if len(cs_i) == 0:
                no_cs.append(j)
            elif (len(cs_i) == 1)|(min(cs_i.values())>j)|(max(cs_i.values())<j):
                k = list(e[e['num_id']==j].sk_kennung)[0] # To take sector id as key.
                # k = j # to take segment id as key. 
                temp.update({k:
                    {list(cs_i.keys())[0]:1}})
            else: 
                 lower = max([x for x in cs_i.values() if x < j])
                 upper = min([x for x in cs_i.values() if x > j])
                 # Easily switch to sectors instead of segmetns, by taking len(...) 
                 # instead of ....sum() for distance measure.
                 # Take average of distance in sections and in segments. 
                 dist_seg_1 = e[(e['num_id']>lower)&(e['num_id']<j)][0].sum()+1
                 dist_seg_2 = e[(e['num_id']<upper)&(e['num_id']>j)][0].sum()+1
                 dist_sec_1 = len(e[(e['num_id']>lower)&(e['num_id']<j)][0])+1
                 dist_sec_2 = len(e[(e['num_id']<upper)&(e['num_id']>j)][0])+1
                 dist1 = (dist_seg_1 + dist_sec_1)/2
                 dist2 = (dist_seg_2 + dist_sec_2)/2             
                 
                 w1 = 1/(dist1) / (1/(dist1) + 1/(dist2))
                 w2 = 1/(dist2) / (1/(dist1) + 1/(dist2))
    
    
                 k = list(e[e['num_id']==j].sk_kennung)[0] # to take sector id as key.
                 # k = j # to take segment id as key.
                 temp.update({k: {list(e[e['num_id']==lower].sk_kennung)[0]: w1,
                                  list(e[e['num_id']==upper].sk_kennung)[0]: w2}})
                
        weight_dict.update({i:temp})        
    return(weight_dict)


def get_sec_length_dict(df_mm):
    
    d = {}
    for i in df_mm['a'].unique():
        sks_i = list(df_mm[df_mm['a']==i].sk_kennung.unique())
        d_t = {}
        for s in sks_i: 
            df_temp_si = df_mm[(df_mm['a']==i)&(df_mm['sk_kennung']==s)]
            # Count number os segments in roads a with sk_kennung s. 
            nsi = df_temp_si.num_id.max() - df_temp_si.num_id.min() + 1 
            d_t.update({s:nsi})
        
        d.update({i:d_t})
        
    return(d)


def load_count_data_merge_and_get_minmax(df_c, pp_in_bisstra, pp_fs, fs_cols1):
    
    df_bisstra = pd.DataFrame.spatial.from_featureclass(pp_in_bisstra)
    # Join vnk & nnk from df_bisstra to the data frame. 
    df_c = df_c.merge(df_bisstra[['sk_kennung','vnk', 'nnk']].drop_duplicates(), 
                      on="sk_kennung", how="left")
    
    df_c = df_c[df_c['id'].notna()]
    df_c['a'] = df_c['id'].apply(lambda x: x[:6]) 
    df_c['a'] = df_c.groupby(['a', 'one_direct_all']).ngroup()
    
    # Make df with only smallest and largest id for each segment by direction. 
    df_mm = pd.concat([df_c.loc[df_c.groupby(['sk_kennung',
                                              'one_direct_all']).num_id.idxmin()],
                       df_c.loc[df_c.groupby(['sk_kennung',
                                              'one_direct_all']).num_id.idxmax()]])
        
    df_mm.sort_values('num_id', inplace=True)
    df_mm[['vnk', 'nnk']] = df_mm[['vnk', 'nnk']].astype(np.float).astype("Int32")
    # Next, match information from Fortschreibung to it.
    df_fs = pd.read_excel(pp_fs, sheet_name="Zeilenformat")
    df_fs = df_fs.rename(columns={"VonNK":"vnk", "NachNK":"nnk"})
    
    df_mm['vnk_nnk'] = df_mm['vnk'].apply(str) + '_'+df_mm['nnk'].apply(str)
    df_fs['vnk_nnk'] = df_fs['vnk'].apply(str) + '_'+df_fs['nnk'].apply(str)
    # Drop nans from df_fs
    df_fs = df_fs[df_fs['DTV19'].isna()==False]
    
    df_mm = df_mm.merge(df_fs[fs_cols1+['vnk_nnk']], on="vnk_nnk", how="left", 
                        suffixes=('','_fs'))
    
    return(df_c, df_mm)


def get_count_paths(pp, date):
    pp_in_bisstra = os.path.join(pp['data_temp'], 'BISSTra_lage.shp')
    pp_fs = os.path.join(pp["data_in_count"], "Fortschreibung_BAB.xlsx")
    pp_in_count = os.path.join(pp['data_out_df'],'zahlst_dir_'+date[:4]+'.csv')
    pp_out_count = os.path.join(pp['data_out_df'], 'counts_interpolated'
                                                   +date[:4]+'.csv')
    return(pp_in_bisstra, pp_fs, pp_in_count, pp_out_count)


def get_newval_dict_existing_stations(df, val_cols):
    '''Get a dictionary with the ids and new values for the ids for which 
    a permanent counting station already exists.'''
    # Make a dict with the values of each counting station. 
    ct_val_dict = {}
    for ct in list(df['DZNr'].unique()[1:]) + [6831,6730]: # Part of manual correct. 
        ct_val_dict.update({ct: {'0': np.nan, '1': np.nan}})
        for d in ['0', '1']:
            if len(df.loc[(df['DZNr']==ct)&(df['one_direct']==int(d))]) > 0:
                ind = df.loc[(df['DZNr']==ct)&(df['one_direct']==int(d))].index[0]
                ct_val_dict[ct].update({d: list(df.loc[ind][val_cols].values)})
    
    
    # Insert counting station number to segments. 
    sk_countstation_dict = {}
    for val in df['sk_kennung'].unique():
        dznr = [int(x) for x in 
                df.loc[(df['AbschnittA']==val)&
                              ((df['TVKfzms'].notnull())|
                               (df['TVSVms'].notnull()))].DZNr.unique()]
        try:
            sk_countstation_dict.update({val: dznr[0]})
        except:
            sk_countstation_dict.update({val: np.nan})
    df.loc[df['DZNr'].isna(), 'DZNr'] = df[df['DZNr'].isna()]['sk_kennung'].replace(sk_countstation_dict)
    # Some manual corrections... 
    manual_ct_ids = {6831: ['6714490020','6714500000'],
                     6730: ['6714490014','6714490015','6714490016','6714490017',
                            '6714490018','6714490019',
                            '6714500001','6714500002','6714500003','6714500004',
                            '6714500005','6714500006','6714500007']}
    for key in manual_ct_ids.keys():
        df.loc[df.index.isin(manual_ct_ids[key]),'DZNr'] = key
    
    
    ct_ids_dict = {}
    for ct in df['DZNr'].unique():
        ids = [i for i in df.loc[df['DZNr']==ct].index]
        ct_ids_dict.update({ct: sorted(ids)})
    ct_ids_dict.update(manual_ct_ids)
    
    new_val_dict = {}
    unclear_ids = []
    for ct, id_list in ct_ids_dict.items():
        temp_dict = {'0': [], '1':[]}
        for i in id_list:
            a = 0
            if type(ct_val_dict[ct]['0']) != float: # float means its nan.
                if df.loc[i].one_direct_all == 0:
                    new_vals = ct_val_dict[ct]['0']
                    new_val_dict.update({i: new_vals})
                    temp_dict.update({'0': temp_dict['0']+[i]})
                    a += 1
            if type(ct_val_dict[ct]['1']) != float:
                if df.loc[i].one_direct_all == 1:
                    new_vals = ct_val_dict[ct]['1']
                    new_val_dict.update({i: new_vals})
                    temp_dict.update({'1': temp_dict['1']+[i]})
                    a += 1
            if a == 0:
                unclear_ids.append(i)
                
    return(df, new_val_dict)
    

def fix_sk_kennung_remove_false_stations(df, val_cols):
    ''' Replace sk_kennung of the frame if counting station has different 
    identifier, if AbschnittA is not zero and the segments are clearly 
    succeeding and remove falsely matched counting stations.'''

    # One manual correction due to error in data (maybe identify generically?).
    df['AbschnittA'].loc[df['DZNr']==6210] = '4926020O4926030O'
    df['aux1'] = df['sk_kennung'].apply(lambda x: re.split(r'\D+', x)[0])
    df['aux2'] = df['sk_kennung'].apply(lambda x: re.split(r'\D+', x)[1])
    df['aux11'] = df[['aux1', 'AbschnittA']
                         ].fillna('-9').apply(lambda x: x[0] in x[1], axis=1)
    df['aux12'] = df[['aux2', 'AbschnittA']
                         ].fillna('-9').apply(lambda x: x[0] in x[1], axis=1)
    df['sk_kennung'].loc[(df['AbschnittA'].notnull()) & 
                                (df['AbschnittA']!=df['sk_kennung']) &
                                ((df['aux11'])|(df['aux12']))
        ] = df['AbschnittA'][(df['AbschnittA'].notnull()) &
                                (df['AbschnittA']!=df['sk_kennung']) &
                                ((df['aux11'])|(df['aux12']))]
    # Remove counting station that are matched falsely. (from neighbor streets) 
    df.loc[((df['AbschnittA'].notnull()) & 
                                (df['AbschnittA']!=df['sk_kennung']) &
                                ((df['aux11'].apply(lambda x: not x))&
                                 (df['aux12'].apply(lambda x: not x)))),
           val_cols+['one_direct','AbschnittA','DZNr']
        ] = np.nan
    df.drop(columns=[x for x in df.columns if x.startswith('aux')], 
            inplace=True)
    
    return(df)
    


def join_data(df, data_dict): 
    '''Add data to df, path as in keys and columns in values of data_dict.'''
    for path in data_dict.keys():
        if path.endswith('frame.csv'):
            new_df=load_csv_and_fix_id(path)
        else:
            new_df = load_csv_and_fix_id(path)
        if data_dict[path] == []:
            cols = list(new_df.columns)
        else:
            cols = data_dict[path]
        df = df.join(new_df[cols], rsuffix='_all')
    return(df)




