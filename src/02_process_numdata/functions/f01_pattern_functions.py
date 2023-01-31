# -*- coding: utf-8 -*-
"""

"""
import pandas as pd
import numpy as np
import os

from src.process_numdata.functions.f04_concat_functions import load_csv_and_fix_id
from src.process_numdata.functions.f02_collapse_functions import select_and_collapse_columns


def identify_links_upd_lanes(lanes_long, link_long, sw_type, ids, max_link_length, 
                   sw_length):
    
    for i in range(len(link_long)):
        # Make a subdataframe with the respective id and if existing, the 
        # previous and the succeeding segments.
        df, r = get_link_row_and_ot_subdf(lanes_long, link_long, i, ids)
        
        new_ids = []
        if (len(df)!=0) & (r.id in df.id.to_list()):
            position, need_split = find_position(df, r)
            # Add split if no lane change near the Entry/Exit point. 
            if need_split == True: 
                df, new_ids = add_missing_split(df, r, position,new_ids)
                position, need_split = find_position(df, r)
            # Select positions of thicker segments and respective lengths.
            pos_list, l_length, len_list_cum, a_lanes, sl_lanes = \
                find_thicker_seg_positions_and_lengths(df, position, r.Type)
            
            # Some labels are wrong because a link segment goes from one 
            # motorway directly into the other. Fix them. 
            if a_lanes > sl_lanes:
                (r, position, need_split, pos_list, l_length, len_list_cum, 
                    a_lanes, sl_lanes, link_long) = reverse_swtype(df, r, 
                                                                   link_long)
            # If the thicker segment part is too long, select only a 
            # segment of length of sw_length and denote it at the slipway.
            if l_length > 400 :
                split_pos, pos_list, need_split = find_split_on_long_sw(
                                     df, pos_list, len_list_cum, sw_length)
                # Split up one segment if no change at around 200m.
                if need_split:
                    df, pos_list, new_ids = split_segment(
                            df, split_pos, pos_list, len_list_cum, 
                            sw_length, r.Type, new_ids)
            # Execute changes in the df.
            df = edit_subdf(df, r, pos_list, sl_lanes, 
                                    a_lanes, sw_type)
            for new_id in new_ids:
                lanes_long = lanes_long.append(df.loc[new_id])
            lanes_long.sort_index(inplace=True)
            
            lanes_long.update(df)
    
    return(lanes_long, link_long)


def find_split_on_long_sw(df, pos_list, len_list_cum, sw_length):
    
    val = min(len_list_cum, key=lambda x:abs(x-sw_length))
    pos = len_list_cum.index(val)
    
    if (sw_length - 15) < val < (sw_length + 15):
        need_split = False
        # Find position of 
        split_pos = None
        pos_list = pos_list[:pos+1]
    # If the closest value is too small, split the next segment
    elif val < sw_length: 
        need_split = True 
        split_pos = pos_list[pos+1]
        pos_list = pos_list[:pos+2]
    # if the closest value is too large, split thsi segment 
    else: 
        need_split = True
        split_pos = pos_list[pos]
        pos_list = pos_list[:pos+1]

    return(split_pos, pos_list, need_split)


def edit_subdf(df, r, pos_list, sl_lanes, a_lanes, sw_type):
    '''Execute the required changes on the df.
    
    '''
    
    # Mark Entries/Exits
    df.loc[pos_list, r.Type] = 1
    df.loc[pos_list, 'sw_type'] = sw_type
    # HERE: Change the number of lanes in these segments
    df.loc[pos_list, 'Slipway_lanes'] = int(max(sl_lanes - a_lanes,1))
    df.loc[pos_list, 'lanes'] = a_lanes    
    
    return(df)


def split_segment(df, split_pos, pos_list, len_list_cum, 
                  sw_length, sw_type, new_ids=[]):
    
    # Determine the length of the slipway before the split position. 
    prev_len = len_list_cum[pos_list.index(split_pos)-1]
    prev_len = prev_len if (pos_list.index(split_pos)-1) >= 0 else 0
    
    # Get the right position of the split within the split-segment. 
    if sw_type == 'Entry':
        cut_len = sw_length - prev_len
    elif sw_type == 'Exit':
        cut_len = df['length'].loc[split_pos] - (sw_length - prev_len)
    else:
        raise TypeError('''sw_type '{}' not recognized! 
                        Must be either 'Entry' or 'Exit'!'''.format(
                                                              sw_type)) 
    # Copy the row that is being split and set its length. 
    # new_row will come first, then teh row that is already in the dataset.
    new_row = df.loc[split_pos]
    new_row.at['length'] = cut_len
    # Adjust length, num, and dist, wherever neccessary. 
    df.loc[(df['id']==split_pos[0]) & (df['num']>split_pos[1]), 'num'] += 1
    df.at[split_pos, 'length'] = df.loc[split_pos, 'length'] - cut_len
    df.at[split_pos, 'dist'] = df.loc[(df['id']==split_pos[0]) & 
                                      (df['num']<split_pos[1]),'length'].sum() \
                               + cut_len
    df.at[split_pos, 'num'] += 1
    
    old_ids = df.index.to_list()
    # Add the new_row, change the index accordingly and sort.
    df = df.append(new_row)
    if sw_type == 'Exit': 
        df.loc[(df['id']==split_pos[0])&(df['num']==split_pos[1]), 'Exit'] = 0
    else:
        df.loc[(df['id']==split_pos[0])&(df['num']==split_pos[1]+1), 'Entry'] = 0
    
    df['num'] = df['num'].astype("int")
    df.set_index(['id', 'num'], drop=False, inplace=True)    
    df.sort_index(inplace=True)
    
    all_ids = df.index.to_list()
    # correct pos_list to only include the segments where the slipway will be.
    if sw_type == 'Entry':
        pos_list_new = [x for x in pos_list if x <= split_pos]
    elif sw_type == 'Exit':
        pos_list = pos_list + [(split_pos[0], split_pos[1]+1)]
        pos_list_new = [x for x in pos_list if x > split_pos]
    
    new_id = [i for i in all_ids if i not in old_ids][0]
    new_ids.append(new_id)
    
    return(df, pos_list_new, new_ids)        


def add_missing_split(df, r, position, new_ids):
    
    old_ids = df.index.to_list()
    new_row = df.iloc[position]
    new_row.at['length'] = r.dist - df['dist'].iloc[position]
    df.at[df.index[position], 'dist'] = r.dist
    df.at[df.index[position], 'length'] = df.iloc[position].loc['length'] - new_row.length
    
    df.loc[(df['id']==new_row.id) & (df['num']>=new_row.num),'num'] +=1
    
    df = df.append(new_row)
    df['num'] = df['num'].astype("int")
    df.set_index(['id', 'num'], drop=False, inplace=True)    
    df.sort_index(inplace=True)
    all_ids = df.index.to_list()
    new_id = [i for i in all_ids if i not in old_ids][0]
    
    new_ids.append(new_id)
    
    return(df, new_ids)            


def reverse_swtype(df, r, link_long):
    
    r.Type = 'Entry' if r.Type == 'Exit' else 'Exit'
    #Repeat the above retrieval of information in this case.
    position, need_split = find_position(df, r)    
    pos_list, l_length, len_list_cum, a_lanes, sl_lanes = \
        find_thicker_seg_positions_and_lengths(df, position, r.Type)
    # Fix the link_long df as well.
    link_long.at[r.name, 'Type'] = r.Type
    
    return(r, position, need_split, pos_list, l_length,  
           len_list_cum, a_lanes, sl_lanes, link_long)

       
def get_long_ot_df(df_ot, var, new_vars):
    
    lanes_df = df_ot[sorted([x for x in df_ot.columns if 
                             x.startswith((var, 'length'))])]
    lanes_df = select_and_collapse_columns(lanes_df, var, rename=True)
    
    # Derive distance from length and add empty Exit&Entry variables where needed. 
    lanes_df = derive_dist_and_add_empty_vars_wide(lanes_df, new_vars)
    
    lanes_long = prepare_df_and_wide_to_long(lanes_df)
    lanes_long['id'] = lanes_long.index.get_level_values('id')
    lanes_long['num'] = lanes_long.index.get_level_values('num')
    lanes_long = insert_from_neighboring_segments_long(lanes_long, var, 
                                                       'missing')
    return(lanes_long)


def insert_from_neighboring_segments_long(long_df, var, replace_val):
    '''Take a dataframe in the long format, replace all values for which the 
    variable defined in var has the value defined in replace_val and replace 
    this value by the next value of the preceding or succeeding segment that is 
    not euqal to the replace_val. 
    
    '''
    
    a = 0
    locs = long_df.loc[long_df[var]==replace_val].index.to_list()
    ids = list(long_df['id'])
    
    while len(locs) >0 and a < 10:
        for i in locs:
            r = long_df.loc[i]
            subids = [str_id(int(r.id)-1) if 
                      str_id(int(r.id)-1) in ids else np.nan]  \
                    + [r.id] \
                    + [str_id(int(r.id)+1) if 
                       str_id(int(r.id)+1) in ids else np.nan]
            df = long_df[long_df['id'].isin(subids)]
            position = np.flatnonzero((df['num']==r.num) & (df['id']==r.id))[0]
            if len([x for x in subids if str(x) != 'nan']) > 1:
                
                if df.iloc[position-1].loc[var] != replace_val:
                    df.at[i, var] = df.iloc[position-1].loc[var]
                elif df.iloc[position+1].loc[var] != replace_val:
                    df.at[i, var] = df.iloc[position+1].loc[var]
                
                long_df.update(df)
            else:
                # Drop if no way to impute lanes.
                long_df.drop(r.id, inplace=True)
        locs = long_df.loc[long_df[var]==replace_val].index.to_list()
        a+=1
                
    return(long_df)


def find_thicker_seg_positions_and_lengths(df, position, sl_type):
    '''
    '''
    
    # Number of lanes on the Autobahn (after the exit).
    if sl_type == 'Exit':
        a_lanes = df.iloc[position+1].lanes
    elif sl_type == 'Entry':
        a_lanes = df.iloc[position-1].lanes
    else:
        raise TypeError('''sl_type '{}' not recognized! 
                        Must be either 'Entry' or 'Exit'!'''.format(
                                                              sl_type))                
    # Find all indices for which the lanes are more than after the Exit
    # and track the length of this segment at the same time.
    sl_lanes = df.iloc[position].lanes
    pos_list = []
    a = position
    l_length = 0
    len_list_cum = []
    # If there is no lane change after the Exit/before the Entry, compare track
    # segments that have the same number of lanes as the a_lanes or more...
    b = a_lanes
    if df.iloc[position].lanes == a_lanes:
        b = a_lanes - 1
    while a >= 0 and a < len(df) and df.iloc[a].lanes > b:
        l_length += df.iloc[a].length
        len_list_cum.append(l_length)
        pos_list.extend(df.iloc[[a]].index.tolist())
        if sl_type == 'Exit':
            a -= 1
        else:
            a += 1
    
    # If the last exit segment is at the same time identified as an Entry point,
    # Take only the preceding segments, that are identified as Entry points.
    # Reason: due to the performed change in the number of lanes due to the 
    # Entry, the slipway cannot be identified properly anymore.   
    if a_lanes <= sl_lanes:     # Only if Exit/Entry was labeled right. 
        if (sl_type=='Exit') & (df.loc[pos_list[0], 'Entry'] == 1):
            pos_new = []
            while len(pos_new)!=len(pos_list):
                pos_new = [x for x in pos_list if (df.loc[x, 'Entry']==1)]
                pos_list = [x for x in pos_new if (pos_list.index(x)==0) |
                                 (pos_list[pos_list.index(x)-1] in pos_new)]
            len_list_cum = len_list_cum[:len(pos_list)]
            l_length = max(len_list_cum)
        
    return(pos_list, l_length, len_list_cum, a_lanes, sl_lanes)
  
    
def find_position(df, r, proximity=50):
    '''Find the position in the dataframe in which the Exit or Entry
    point lies. 
    
    Note that the Exit segment position denotes the last slipway-segment so
    that r.dist corresponds to df
    while the Entry segment denotes the first segment with a slipway.
    
    Check also, whether the matched position is close to the given exit/entry
    position and if it isn't, set need_split to True and return the position 
    of the segment that needs splitting instead. 
    '''
    
    need_split = False
    try:
        position = np.flatnonzero((df['dist']==r.dist) & (df['id']==r.id))[0]
        if r.Type == 'Exit':
            # Deduct one to get the last segment containing the slipway. 
            position -= 1
    except:
        # If couldn't find a perfect match, change dist before finding next
        # best position to take into account breaks between segments.
        # Do not select the break points because nothing ever changes here.
        t_df = df[df['num']!=0]
        t_df.loc[t_df['id']<r.id, 'dist'] -= 500
        t_df.loc[t_df['id']>r.id, 'dist'] += 500
        
        if len(t_df) == 0:
            position = np.flatnonzero(df['id']==r.id)[0]
            need_split = True
        else:            
            t = t_df[t_df['dist']==min(t_df['dist'], 
                     key=lambda x:abs(x-r.dist))].index[0]
            dist_list = [abs(x - r.dist) for x in t_df['dist']]
            position = df.index.get_loc(t)
            
            if r.Type == 'Exit':
                position -=1
            
            if min(dist_list) >= proximity:
                need_split = True
                position = max(np.flatnonzero((df['dist']<r.dist) & 
                                              (df['id']==r.id)))
        
    return(position, need_split)


def get_link_row_and_ot_subdf(lanes_long, link_long, i, ids): 
    
    if i in range(0,len(link_long),1000):
        print(i)
    r = link_long.iloc[i]
    subids = [str_id(int(r.id)-1) if 
              str_id(int(r.id)-1) in ids else np.nan] \
           + [r.id] \
           + [str_id(int(r.id)+1) if 
              str_id(int(r.id)+1) in ids else np.nan]    
    df = lanes_long[lanes_long['id'].isin(subids)]  
    
    return(df, r)


def str_id(index):        
    return('0'*(10-len(str(index))) + str(index))


def get_pattern_paths(date, pp):

    pp_in = os.path.join(pp['data_out_df'], 'ot_df_' +date+'.csv')
    pp_zaehl = os.path.join(pp['data_out_df'], 'zahlst_dir_'+date[:4]+'.csv')
    pp_out_ot = os.path.join(pp['data_out_df'], 'ot_final'+date+'.csv')
    pp_out_zaehl = os.path.join(pp['data_out_df'], 
                                'zaehlst_interpolated'+date[:4]+'.csv')
    pp_links = os.path.join(pp['data_out_df'], 
                            'links_on_lanes_'+date[:4]+'.csv')
    pp_main_points_upd = os.path.join(pp['data_out_df'], 
                                  'main_link_points_fixed_'+date[:4]+'.csv')
    
    return(pp_in, pp_zaehl, pp_out_ot, pp_out_zaehl, 
           pp_links, pp_main_points_upd)


def read_data_select_cols_and_relabel(pp_in, id_length, remove_col_starts, 
                                      relabel_col_starts, relabel_dict):
    
    df = pd.read_csv(pp_in, index_col=0)
    df['new_index'] = ['0'*(id_length-len(str(x)))+str(x) for x in df.index]
    df.set_index('new_index', inplace=True)
    df.index.name = None
    cols = list(df.columns)
    
    remove_cols = [x for x in cols if x.startswith(remove_col_starts)]
    df = df.drop(remove_cols, axis='columns')
    
    # Relabel missings to 0 and yes to 1 for defined variables
    relable = [x for x in cols if x.startswith(relabel_col_starts)]
    df.update(df[relable].replace(relabel_dict))
    
    return(df)


def get_and_prep_link_data(sw_type, pp, last_ids):
        
    pp_links = os.path.join(pp['data_out_df'], sw_type+'_links.csv')
    link_long = prepare_df_and_wide_to_long(pp_links)
    link_long['dist'] = round(link_long['dist']) 
    link_long['id'] = link_long.index.get_level_values('id')
    # Remove Entries at the beginning and the end of a subsequent road segment.
    link_long.drop(link_long[(link_long['id'].apply(lambda x: x[-4:])=='0000')&
                             (link_long['dist']==0)].index, inplace=True)
    link_long.drop(link_long[(link_long['id'].apply(lambda x: x in last_ids))&
                             (link_long['dist']>=480)].index, inplace=True)    
    return(link_long)


def derive_dist_and_add_empty_vars_wide(df, empty_vars):
    '''Derive the distance from the length in a wide dataframe and add fields
    only where there are (not nan) entries in for the other variables.
    '''
    
    len_cols = [x for x in df.columns if x.startswith('length')]
    width = len(len_cols)
    
    for i in range(width):
        if i == 0:
            df['dist_'+str(i)] = 0
        else:
            df['dist_'+str(i)] = round(
                        df[[x for x in len_cols if int(x[-1])<i]].sum(axis=1))
        df['dist_'+str(i)][df['length_'+str(i)].isnull()]=np.nan
        
        for var in empty_vars:
            df[var+'_'+str(i)] = np.nan
            df[var+'_'+str(i)][df['length_'+str(i)].notnull()] = 0

    return(df)
    

def prepare_df_and_wide_to_long(df):
    '''df can be a path to a dataframe or a dataframe.
    
    '''
    
    if type(df) == str:
        df = load_csv_and_fix_id(df)
    cols = list(set(['_'.join(x.split('_')[:-1])+'_' for x in df.columns 
                     if x[-1].isnumeric()]))
    df['id'] = df.index
    df_long = pd.wide_to_long(df, cols, i='id', j='num')
    
    df_long.columns=[x[:-1] if x[-1]=='_' else x for x in df_long.columns]
    df_long.sort_index(inplace=True)
    df_long.dropna(axis=0, how='any', inplace=True)
    
    return(df_long)


def make_id_sequence_dict(df, id_length, in_seg_id_length):
    
    id_dict = {}
    sep = id_length - in_seg_id_length
    for ind in df.index:
        if ind[:sep] in id_dict.keys():
            id_dict[ind[:sep]].append(ind)
        else:
            id_dict.update({ind[:sep]:[ind]})
    
    for key in list(id_dict.keys()):
        i = 0
        old_pos=0
        for pos in range(1, len(id_dict[key])):
            if int(id_dict[key][pos-1]) != int(id_dict[key][pos])-1:
                id_dict.update({key+'_'+str(i): id_dict[key][old_pos:pos]})
                old_pos = pos
                i+=1
        if i>0:
            id_dict.update({key+'_'+str(i): id_dict[key][old_pos:pos+1]})
            id_dict.pop(key)
    
    return(id_dict)
  







