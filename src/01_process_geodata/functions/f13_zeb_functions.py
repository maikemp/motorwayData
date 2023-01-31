# -*- coding: utf-8 -*-
"""

"""

import os
import arcpy
import pandas as pd
import numpy as np
import math
from numbers import Number
from src.process_geodata.functions.f10_ramps_functions import locate_points_and_to_df
from src.process_numdata.functions.f04_concat_functions import load_csv_and_fix_id


arcpy.env.overwriteOutput = True

def aggregate(df, lanes, newvar_dict):
    
    newvars = []
     
    for var, stats in newvar_dict.items():
        if var == 'bauw':
            for stat in stats:
                df['bauw_'+stat] = df[cols(df,('bauw','length'))].apply(
                        bw_stat, args=[stat], axis=1)
                newvars.append('bauw_'+stat)
        else:
            for stat in stats:
                # If lanes=1 instead of 4, only primary lane woulb be considered. 
                df[var+'_'+stat]=df[cols(df,(var,'length'))].apply(
                                             w_stat, args=[stat,lanes], axis=1)
                newvars.append(var+'_'+stat)
    return(df, newvars)



def w_stat(r, stat, lanes):
    ''' stats: must be in ['mean', 'var', 'max', 'min', 'num'] for Mean, 
        Variance, Minimmum, Maximum or number of considered values. 
        lanes: gives number of lanes (maximally considered in computations). 
               E.g. if set to 1, only primary lane will be considered. 
    '''
    
    weights = []
    vals = []
    out = np.nan
    # Find varname. 
    var = '_'.join([x for x in r.keys() if not 
                    x.startswith('length_')][0].split('_')[:2])
    for n, l in enumerate(r[[x for x in r.keys() if x.startswith('length')]]):
        for i in range(1,lanes+1):
            if isinstance(r[var+'_'+str(i)+'_'+str(n)], Number):
                if r[var+'_'+str(i)+'_'+str(n)] > 0:
                    weights.append(r['length_'+str(n)])
                    vals.append(r[var+'_'+str(i)+'_'+str(n)])
    if len(weights)>0:
        if stat == 'mean':
            out = np.average(vals, weights=weights)
    
        elif stat == 'var':
            av = np.average(vals, weights=weights)
            out = np.average((vals-av)**2, weights=weights)
        
        elif stat == 'max':
            out = max(vals)
        
        elif stat == 'min':
            out = min(vals)
            
        elif stat == 'num':
            out = len(vals)
              
    return(out)

def cols(df, var_starts):
    return([x for x in df.columns if x.startswith(var_starts)])

def bw_stat(r, stat):
    A_length = 0
    B_length = 0
    
    out_stat = np.nan
    
    for n, l in enumerate(r[[x for x in r.keys() if x.startswith('length')]]):
        for i in range(1,5):
            if r['bauw_3_'+str(i)+'_'+str(n)] == 'A':
                if l>0:
                    A_length += l
            elif r['bauw_3_'+str(i)+'_'+str(n)] == 'B':
                if l>0:
                    B_length += l
    if stat == 'mode':
        if A_length > B_length:
            out_stat = 'A'
        elif A_length < B_length:
            out_stat = 'B'
        
        # If they have equal length, take bauw of primary lane. 
        if (A_length == B_length)&(A_length != 0):
            A_length = 0
            B_length = 0
    
            for n, l in enumerate(r[[x for x in r.keys() if 
                                     x.startswith('length')]]):
                if r['bauw_3_1_'+str(n)] == 'A':
                    if l>0:
                        A_length += l
                elif r['bauw_3_1_'+str(n)] == 'B':
                    if l>0:
                        B_length += l
            if A_length > B_length:
                out_stat = 'A'
            elif A_length < B_length:
                out_stat = 'B'
    elif stat == 'het':
        if (A_length > 0)&(B_length > 0):
            out_stat = 1
        elif (A_length > 0)|(B_length > 0): 
            out_stat = 0
    
    return(out_stat)


def get_length(r, ncols):
    
    length = []
    r = list(r)
    r = [x for x in r if x==x]
    r.append(500)
    
    for n, i in enumerate(r): 
        if i != 500:
            l = r[n+1]-i
            if l > 120:
                l=100
            length.append(l)
            
    length = length + (ncols-len(length))*[np.nan]
            
    return(pd.Series(length))


def add_startpoints(df, id_fields, set_to_start):
    
    # Set 'FMEAS' to 0 when it is close to it according to set_to_start.
    df['FMEAS'][(df['FMEAS']<set_to_start)&(df['FMEAS']>0)]=0
    df_0 = df[['ID', 'FMEAS']+id_fields]
    # Keep only the point closest to the start.   
    df_0 = df_0.loc[df_0.groupby('ID')['FMEAS'].idxmin()]
    # Drop if 'FMEAS' close to 0 since those startpoints are already part of df. 
    df_0 = df_0[df_0['FMEAS']>=set_to_start]

    # Separate by direction
    df_0_l = df_0[df_0['l_r']=='L']
    df_0_r = df_0[df_0['l_r']=='R'] 
    
    # For right df, find zeb info that matches 'vnk', 'nnk', 'l_r' and has the 
    #    next smaller 'sk_dist'.
    df_0_r_m  = df_0_r[df_0_r['sk_dist']>=100]
    df_0_r_m['sk_dist'] -= 100 
    # Merge df_0_r_m with mutated sk_dist with the columns from df. These must
    # be what happens at the start of each segment. 
    df_0_r_f = df_0_r_m.merge(df[[x for x in df.columns if x not in 
                                  ['ID', 'FMEAS']]], how='left', on=id_fields)
    # Take those with successfull merge.
    df_0_r_f = df_0_r_f[df_0_r_f['bauw_3_1'].notna()]
    # Set FMEAS to 0 for Start of frame segment. 
    df_0_r_f['FMEAS'] = 0
    # Drop duplicates to be ready to be merged back to df
    df_0_r_f.drop_duplicates(subset=id_fields, inplace=True)
    # Now take df_0_r that has a smaller sk_dist. 
    df_0_r_m  = df_0_r[df_0_r['sk_dist']<100]
       
    # Take only maximum distance in bisstra section and in frame 
    a = df[df['l_r']=='R'].loc[df[df['l_r']=='R'].groupby('ID')['FMEAS'].idxmax()]
    # a = a.loc[a.groupby('sk_kennung')['sk_dist'].idxmax()]
    # Derive preceeding frame segment.
    df_0_r_m['ID_s'] = df_0_r_m['ID'].apply(int).apply(lambda x: x-1).apply(
            str).apply(lambda x: (10-len(x))*'0'+x)
    # Merge columns from a if their nnk is vnk in frame segment of interest
    # and if it is from the preceeding frame segment. 
    df_0_r_m = df_0_r_m.merge(
            a[[x for x in a.columns if x not in ['ID_s','FMEAS']]], how='left', 
           left_on=['ID_s','l_r'], right_on=['ID','l_r'])
    #prev. left_on=['vnk','ID_s','l_r'], right_on=['nnk','ID','l_r'])
    # Clean up. Take info from bisstra from previous segment and ID from frame.
    df_0_r_m.drop(columns=['vnk_x', 'nnk_x', 'sk_dist_x', 'ID_y', 'ID_s'], 
                  inplace=True)
    df_0_r_m.rename(columns={'vnk_y':'vnk', 'nnk_y':'nnk', 
                             'sk_dist_y':'sk_dist','ID_x':'ID'}, inplace=True)
    # Mark as startpoint. 
    df_0_r_m['FMEAS'] = 0
    df_0_r_m = df_0_r_m[df_0_r_m['bauw_3_1'].notna()]

    
    # Append the df with this
    df = df.append([df_0_r_m, df_0_r_f])
    df.drop_duplicates(subset=['ID', 'FMEAS'], inplace=True)
        
    # Continue with Left.
    # For left df, find zeb info that matches 'vnk', 'nnk', 'l_r' and has the 
    #    the next larger 'sk_dist'. 
    df_0_l['sk_dist'] += 100 
    # Columns where successful merge could be made with previous segment.
    df_0_l = df_0_l.merge(df[[x for x in df.columns if x not in 
                                  ['ID', 'FMEAS']]], how='left', on=id_fields)
    # Take those with successfull merge.
    df_0_l_f = df_0_l[df_0_l['bauw_3_1'].notna()]
    # Set FMEAS to 0 for startpoint. 
    df_0_l_f['FMEAS'] = 0
    # Drop duplicates to be ready to be merged back to df.
    df_0_l_f.drop_duplicates(subset=id_fields, inplace=True)
    # Now, take the rows of df_0_l where nothing was merged. 
    df_0_l = df_0_l[df_0_l['bauw_3_1'].isna()][[x for x in df_0.columns]]
    # Reverse sk_dist back to what it was before. 
    df_0_l['sk_dist'] -= 100 
    
    # Take only minimuum distance in bisstra section and maximum in frame 
    a = df[df['l_r']=='L'].loc[df[df['l_r']=='L'].groupby('ID')['FMEAS'].idxmax()]
    # a = a.loc[a.groupby('sk_kennung')['sk_dist'].idxmin()]
    # Derive preceeding segment
    df_0_l['ID_s'] = df_0_l['ID'].apply(int).apply(lambda x: x-1).apply(
            str).apply(lambda x: (10-len(x))*'0'+x)
    # Merge columns from a if their vnk is nnk in frame segment of interest
    # and if it is from the preceeding frame segment.     
    df_0_l = df_0_l.merge(
            a[[x for x in a.columns if x not in ['ID_s','FMEAS']]], how='left', 
           left_on=['ID_s','l_r'], right_on=['ID','l_r'])
    #prev. left_on=['nnk','ID_s','l_r'], right_on=['vnk','ID','l_r'])

    # Clean up. Take info from bisstra from previous segment and ID from frame.
    df_0_l.drop(columns=['vnk_x', 'nnk_x', 'sk_dist_x', 'ID_y', 'ID_s'], 
                inplace=True)
    df_0_l.rename(columns={'vnk_y':'vnk', 'nnk_y':'nnk', 
                             'sk_dist_y':'sk_dist','ID_x':'ID'}, inplace=True)
    # Make it a startpoint.
    df_0_l['FMEAS'] = 0
    df_0_l = df_0_l[df_0_l['bauw_3_1'].notna()]
    
    # Append the df with this
    df = df.append([df_0_l, df_0_l_f])
    df.drop_duplicates(subset=['ID', 'FMEAS'], inplace=True)
    
    return(df)


def merge_zeb_to_df(zeb_df, df_dir, col_map, value_cols):
    
    for fs in zeb_df['fs'].unique():
        
        zeb_temp = zeb_df[zeb_df['fs']==fs]
        rename_dict = {col: col+'_'+str(fs) for col in value_cols}
        zeb_temp.rename(columns=rename_dict, inplace=True)
        
        df_dir = df_dir.merge(zeb_temp, how='left', 
                              left_on=list(col_map.keys()),
                              right_on=list(col_map.values()))     
    return(df_dir)


def get_and_clean_dfs(zeb_df, df_r, df_l, col_map_r, col_map_l, ident_cols):
    
    zeb_df['vnk'] = zeb_df['vnk'].apply(str).apply(
                        cut_letters).apply(shorten_nk)
    zeb_df['nnk'] = zeb_df['nnk'].apply(str).apply(
                        cut_letters).apply(shorten_nk)
    df_r['vnk'] = df_r['vnk'].apply(str)
    df_r['nnk'] = df_r['nnk'].apply(str)
    df_l['vnk'] = df_l['vnk'].apply(str)
    df_l['nnk'] = df_l['nnk'].apply(str)
    
    # Subset zeb_df by 'fs' and 
    # then merge zeb over vnk, nnk, lage and vst (for lage=R) or bst (for lage=L) 
    # Before transform zeb:bst,vst and df_l/df_r:fmeas to int or so... 
    # (larger tolerance for endpoint of lage=L may be needed) (fuzzy matching or so??)
    zeb_df_r = zeb_df[zeb_df['lage']=="R"]
    zeb_df_l = zeb_df[zeb_df['lage']=="L"]
    # Where possible, correct bst by fmeas to enable merge over this. 
    zeb_df_l = adjust_L_bst(zeb_df_l, df_l)

    zeb_df_r.drop_duplicates(subset=[x for x in ident_cols if x!='bst'], 
                             inplace=True)
    zeb_df_l.drop_duplicates(subset=[x for x in ident_cols if x!='vst'], 
                             inplace=True)
    
    return(zeb_df_r, zeb_df_l, df_r, df_l)


def adjust_L_bst(zeb_df_l, df_l):
    
    # Select only sections with 'fs'==1 and odd endpoint in zeb data.
    zeb_bst = zeb_df_l[['vst', 'vnk', 'nnk', 'bst', 'gw_15']][
            (zeb_df_l['bst'].apply(str).apply(lambda x: x[-2:])!='00')&
            (zeb_df_l['fs']==1)]
    # Search only odd startpoints in the point data.
    dfl_fmeas = df_l[['l_r', 'vnk', 'nnk', 'fmeas']][
            df_l['fmeas'].apply(str).apply(lambda x: x[-2:])!='00']
    # Merge them together on vnk and nnk (enough since only one row per sector.)
    zeb_bst = dfl_fmeas.merge(zeb_bst, how='outer', on=['vnk', 'nnk'])
    # Keep only those that need fixing due to different end point distance, 
    # otherwise matchable
    zeb_bst = zeb_bst[zeb_bst['fmeas']!=zeb_bst['bst']].dropna()
    zeb_bst.drop_duplicates(inplace=True)
    zeb_bst = zeb_bst[zeb_bst['vst']<zeb_bst['fmeas']]
    
    zeb_bst = zeb_df_l.merge(zeb_bst[['fmeas','bst','vnk','nnk']], 
                       on=['bst','vnk','nnk'], how='left')

    zeb_bst['new'] = zeb_bst[['bst', 'fmeas']].apply(
        lambda x: x.fmeas if not math.isnan(x.fmeas) else x.bst, axis=1)
    zeb_bst['bst']=zeb_bst['new']
    zeb_bst.drop(columns=['fmeas','new'], inplace=True)
    
    return(zeb_bst)


def cut_letters(a):
    if not a[-1].isnumeric(): 
        a = a[:-1]
    return(a)


def shorten_nk(a):
   
    if len(a)==9:
       if '0000' in a:
           a = a.replace('0000', '00')
       elif '000' in a:
           a = a.replace('000', '0')
       elif a.count('00') == 2:
           a = a.replace('00', '0')
       elif '00' in a:
           a = a.replace('00', '')

    return(a)


def choose_coordinates_and_subset_df(sdf):
    
    sdf['Compass_dist_0'] = sdf[['compass_0','compass_a']].apply(
            compass_distance, axis=1
    )
    sdf['Compass_dist_1'] = sdf[['compass_1','compass_a']].apply(
            compass_distance, axis=1
    )
    
    sdf['Near_dir'] = sdf[['NEAR_DIST_0', 'NEAR_DIST_1', 
                           'Compass_dist_0', 'Compass_dist_1']].apply(
            choose_nearpoint, axis=1
    )
    
    sdf['Near_X'] = sdf[['NEAR_X_0', 'NEAR_X_1', 'Near_dir']].apply(
            get_nearcoord, axis=1
    )
    sdf['Near_Y'] = sdf[['NEAR_Y_0', 'NEAR_Y_1', 'Near_dir']].apply(
            get_nearcoord, axis=1
    )
        
    return(sdf[['sk_kennung','str_kennun', 'sk_laenge_', 'vnk', 
                'nnk', 'l_r', 'fmeas', 'Near_X', 'Near_Y']])

def get_nearcoord(r):
    
    pt = np.nan
    
    for key in r.keys():
        if pd.notnull(r.Near_dir):
            if str(int(r.Near_dir)) == key[-1]:
                pt = r[key]
    
    return(pt)


def choose_nearpoint(r):
    
    dist_0 = r[0]
    dist_1 = r[1]
    c_0 = r[2]
    c_1 = r[3]
    
    nearpt = np.nan
    
    if min(dist_0, dist_1) == dist_0:
        # If 0 is closer, check if compass distance matches okay.
        if c_0 < 90:
            nearpt = 0
        # If not, check if compass distance of point 1 matches well.
        else: 
            if c_1 < 45:
                nearpt = 1
    elif min(dist_0, dist_1) == dist_1:
        # If 1 is closer, check if compass distance matches okay.
        if c_1 < 90:
            nearpt = 1
        # If not, check if compass distance of point 0 matches well.
        else: 
            if c_0 < 45:
                nearpt = 0
        
    return(nearpt)


def compass_distance(a, abs_val=True):
    ''' Derive distance between point x and y on a compass.
    '''
    x = a[0]
    y = a[1]
    
    if max(x,y) > 360:
        raise ValueError(str(max(x,y)),"is not a compass value. Must lie between 0 and 360.")
    if min(x,y) < 0:
        raise ValueError(str(min(x,y)),"is not a compass value. Must lie between 0 and 360.")   
    
    dist = x-y
    if abs(dist) > 180:
        dist = np.sign(y-x)*(360 - abs(x-y)) 

    if x == y:
        dist = 0
        
    if abs_val == True:
        dist = abs(dist)
        
    return(dist)


def find_nearpoints_to_sdf_by_A(in_pts, in_lines_dict, tolerance, sdf_in):

    sdf_sr = pd.DataFrame.spatial.from_featureclass(in_pts)
    sdf_sr['str_kennun'] = sdf_sr['str_kennun'].apply(lambda a: a[0]+' '*(1-a.count(' '))+(4-len(a))*'0'+a[1:])
    sdf_sr.spatial.to_featureclass(in_pts)
    
    
    df = pd.DataFrame()
    for A in sdf_in['ref'].unique():
           
        arcpy.management.MakeFeatureLayer(in_pts, 
                                           r"in_memory\pts_A",
                                           "str_kennun = '{}'".format(A))
        arcpy.management.CopyFeatures(r"in_memory\pts_A", r"in_memory\A")
        for key in in_lines_dict.keys():
            arcpy.management.MakeFeatureLayer(r"in_memory\frame_"+key, 
                                               r"in_memory\frame_"+key+"_A",
                                               "ref = '{}'".format(A))
        sdf_A = sdf_in[sdf_in['ref']==A]
        
        sdf_np_A = find_nearpoints_and_to_sdf(r"in_memory\A", in_lines_dict,
                                                 tolerance, sdf_A)
        df = df.append(sdf_np_A, ignore_index=True)
    
    return(df)


def find_nearpoints_and_to_sdf(in_pts, in_lines_dict, tolerance, sdf_in):
    
    for key in in_lines_dict.keys():
        # Find closest point on frame_0. Fields are added to input directly.
        arcpy.analysis.Near(in_pts, in_lines_dict[key], 
                            tolerance, "LOCATION", "NO_ANGLE", "PLANAR")
        for field in ['NEAR_FID', 'NEAR_DIST', 'NEAR_X', 'NEAR_Y']:
            arcpy.management.AlterField(in_pts, field, field+'_'+key)        
    
    sdf_new = pd.DataFrame.spatial.from_featureclass(in_pts)

    join_columns = ['FID', 'one_direct', 'compass', 'ID']
    sdf_new = sdf_new.merge(sdf_in[join_columns], left_on='NEAR_FID_0', 
                            right_on='FID', how='left'
                            ).rename(columns={k:k+'_0' for k in join_columns}
                ).merge(sdf_in[join_columns], left_on='NEAR_FID_1', 
                right_on='FID', how='left').rename(
                        columns={k:k+'_1' for k in join_columns})
    return(sdf_new)


def load_frame_add_dir_by_direction(pp_frame, ldm_source):
    
    
    sdf = pd.DataFrame.spatial.from_featureclass(pp_frame).set_index('ID')
    zs_df = load_csv_and_fix_id(ldm_source, encoding='latin1')[['compass']]
    sdf = sdf.join(zs_df)
    sdf['ID'] = sdf.index
    
    arcpy.management.MakeFeatureLayer(pp_frame, r"in_memory\frame_0", 
                                      "one_direct = '0'")
    arcpy.management.MakeFeatureLayer(pp_frame, r"in_memory\frame_1", 
                                      "one_direct = '1'")
    print(r"Frame direction 0 saved in r'in_memory\frame_0'","\n"
          r"Frame direction 1 saved in r'in_memory\frame_1'")
    
    return(sdf)


def join_to_startpts(split_lines, split_pts, out_pts):
    
    arcpy.management.FeatureVerticesToPoints(split_lines, 
                                             r"in_memory\start_pts", "START")
    
    arcpy.analysis.SpatialJoin(r"in_memory\start_pts", split_pts, 
                               out_pts, "JOIN_ONE_TO_MANY", 
                               "KEEP_ALL", match_option="WITHIN_A_DISTANCE", 
                               search_radius="0.5 Meters")
    
    sdf = pd.DataFrame.spatial.from_featureclass(out_pts)
    # Only correctly matched 
    sdf = sdf[(sdf['sk_kennung_1']==sdf['sk_kennung'])&(sdf['Join_Count']>=1)]
    sdf = sdf[[x for x in sdf.columns if not x.endswith('1')]]
    
    sdf.spatial.to_featureclass(out_pts)
    
    print("Point Layer of start points of", "r'"+str(split_lines)+"'", 
          "with information on Distance along segment from", 
          "r'"+str(split_pts)+"'","in 'fmeas'",
          "and information on directional mean from", 
          "r'"+str(split_lines)+"'", "in 'compass_a'",       
          "saved out to:", "r'"+str(out_pts)+"'")
    
    return()


def split_line_at_points_and_add_dir_mean(in_lines, cut_points, out_lines):
    
    # Cut lines at points.
    arcpy.management.SplitLineAtPoint(in_lines, cut_points, out_lines,
                                      "0.1 Meters")
    print("Finished splitting", "r'"+str(in_lines)+"'", 
          "on", "r'"+str(cut_points)+"'")
    # Add a unique identifier.
    arcpy.AddField_management(out_lines, "id_temp", "Long")
    c = arcpy.da.UpdateCursor(out_lines, ["id_temp"])
    for i, r in enumerate(c):
        r[0] = i
        c.updateRow(r)    
    arcpy.stats.DirectionalMean(out_lines, r"in_memory\dir_mean", 
                                "DIRECTION", "id_temp")
    
    # Merge this information back to the point feature layers.
    sdf_s = pd.DataFrame.spatial.from_featureclass(out_lines)
    sdf_dm = pd.DataFrame.spatial.from_featureclass(r"in_memory\dir_mean")    
    
    # Prepare merge to start/end points.
    sdf_s.set_index('id_temp', inplace=True, drop=False)
    sdf_dm.set_index('id_temp', inplace=True, drop=False)
    
    # Add the variable of interest to start/end points. 
    sdf_s = sdf_s.join(sdf_dm['CompassA'])
    
    # Back to feature class.
    sdf_s.spatial.to_featureclass(out_lines)  
    
    print("Layer of lines from", "r'"+str(in_lines)+"'", "that where split at", 
          "r'"+str(cut_points)+"'", 
          "with information on directional mean of segment in 'CompassA'",
          "saved as:", "r'"+str(out_lines)+"'")
    return()


def generate_and_locate_points(in_lines, out_points, spacing="100 Meters"):
    
    # Create point every 100 meters along the line. Includes start and end point!
    arcpy.management.GeneratePointsAlongLines(in_lines, out_points, "DISTANCE", 
                                              spacing, None, "END_POINTS")

    # Make sure there is a unique identifier. 
    arcpy.AddField_management(out_points, "id_temp", "Long")
    c = arcpy.da.UpdateCursor(out_points, ["id_temp"])
    for i, r in enumerate(c):
        r[0] = i
        c.updateRow(r)

    # Get distance from start along lines.
    df_pts_loc = locate_points_and_to_df(
            in_lines, "0,1 Meters",  out_points, 
            ['id_temp', 'sk_kennung', 'sk_laenge_', 'vnk', 'nnk', 'l_r'], 
            "nr_", False, "ALL"
    )
    
    # delete rows if nr_!=nr_2?
    df_pts_loc = df_pts_loc[df_pts_loc['nr_']==df_pts_loc['nr_2']]
    
    # Merge this information back to the point feature layers.
    sdf_pts = pd.DataFrame.spatial.from_featureclass(out_points)
    # Prepare merge to start/end points.
    sdf_pts.set_index('id_temp', inplace=True)
    df_pts_loc.set_index('id_temp', inplace=True)
    # Add the variable of interest to start/end points. 
    sdf_pts = sdf_pts.join(df_pts_loc['FMEAS'])
    
    sdf_pts.spatial.to_featureclass(out_points)
    
    print("Point Layer of points along line every", spacing, 
          "with information on Distance along segment in 'FMEAS'",
          "saved out to:", "r'"+str(out_points)+"'")
    
    return()


def find_side_and_join_to_bisstra(sdf_startpts, sdf_endpts, bisstra_in):
    '''
    '''
    # Join the dfs on OBJECTID (from bisstra line) and give suffix 'v' or 'n'. 
    sdf_startpts = sdf_startpts.set_index('OBJECTID').add_suffix('_v')
    sdf_endpts = sdf_endpts.set_index('OBJECTID').add_suffix('_n')
    vn_df = pd.concat([sdf_startpts, sdf_endpts], axis=1, join="inner")
    # Obsolete to have twice: SK_Kennung.
    vn_df = vn_df.rename(columns={'Sk_Kennung_v':'Sk_Kennung'}
                                         ).drop(columns='Sk_Kennung_n')
    
    # Next, decide, whether to join left or right zeb values. 
    vn_df['vnk'] = vn_df['Sk_Kennung'].apply(lambda x: x[0: 7])
    vn_df['nnk'] = vn_df['Sk_Kennung'].apply(lambda x: x[8: 15])
    vn_df['L_R'] = 0
    
    # Decide Whether the segment is left or right of the stationary direction. 
    vn_df['L_R'][(vn_df['vnk']==vn_df['NK_Kennung_v']) & 
                 (vn_df['nnk']==vn_df['NK_Kennung_n'])] = "R"
    vn_df['L_R'][(vn_df['vnk']==vn_df['NK_Kennung_n']) & 
                 (vn_df['nnk']==vn_df['NK_Kennung_v'])] = "L"
    vn_df['L_R'][((vn_df['vnk']==vn_df['NK_Kennung_v']) |
                  (vn_df['nnk']==vn_df['NK_Kennung_n'])) 
                 & (vn_df['L_R']==0)] = "R"
    vn_df['L_R'][((vn_df['vnk']==vn_df['NK_Kennung_n']) | 
                  (vn_df['nnk']==vn_df['NK_Kennung_v'])) 
                 & (vn_df['L_R']==0)] = "L"
    
    # Check if length of Kennungen fits: NK_Kennungen: 7, Sk_Kennung: 16
    #vn_df['NK_Kennung_n'][vn_df['NK_Kennung_n'].apply(
    #              lambda x: x is not None)].apply(str).apply(len).min()
    # Futher checks...
    #len(vn_df[vn_df['L_R']==0]) # should be 8
    #len(vn_df[vn_df['L_R']=='R']) # should be 3075
    #len(vn_df[vn_df['L_R']=='L']) # should be 2953
    
    # Add this information to the bisstra network. 
    sdf_bisstra = pd.DataFrame.spatial.from_featureclass(bisstra_in)
    sdf_bisstra.set_index('OBJECTID', inplace=True)
    sdf_bisstra = pd.concat([sdf_bisstra, vn_df[['vnk', 'nnk', 'L_R']]], 
                            axis=1, join='inner')
    
    return(sdf_bisstra)



def take_line_start_end_and_join_pts(in_line, join_pts, search_radius):
    '''Take the start and end vertice of in_line, in separate layers, join 
    information from closest point in join_pts within search_radius to them,
    save them into a spatially enabled arcpy dataframe and select required 
    columns (currently hard coded...).
    
    '''
    # Take end and start points from bisstra sections in separate layers. 
    arcpy.management.FeatureVerticesToPoints(
            in_line, r"in_memory\start_pts", "START")
    arcpy.management.FeatureVerticesToPoints(
            in_line, r"in_memory\end_pts", "END")
    
    # Join info from closest NK_Point to start and end points.
    arcpy.analysis.SpatialJoin(
            r"in_memory\start_pts", join_pts, 
            r"in_memory\NK_start_Join", "JOIN_ONE_TO_ONE", "KEEP_ALL", 
            match_option="CLOSEST", search_radius="71 Meters")
    arcpy.analysis.SpatialJoin(
            r"in_memory\end_pts", join_pts, 
            r"in_memory\NK_end_Join", "JOIN_ONE_TO_ONE", "KEEP_ALL", 
            match_option="CLOSEST", search_radius="71 Meters")
    
    # Make pandas df from joined features, keeping only relevant columns.
    sdf_startpts = pd.DataFrame.spatial.from_featureclass(r"in_memory\NK_start_Join")
    sdf_startpts = sdf_startpts[['Join_Count', 'OBJECTID', 'Sk_Kennung', 'NK_Kennung']]
    sdf_endpts = pd.DataFrame.spatial.from_featureclass(r"in_memory\NK_end_Join")
    sdf_endpts = sdf_endpts[['Join_Count', 'OBJECTID', 'Sk_Kennung', 'NK_Kennung']]
    
    return(sdf_startpts, sdf_endpts)



def get_paths(pp):
    
    # Only one measurement point per segment. Thus, we can ignore years.
    zeb_files = [os.path.join(pp['data_in_zeb'],x) for x in 
                     os.listdir(pp['data_in_zeb']) if x.endswith('.xlsx')]
    pp_frame = os.path.join(pp['data_out_km'], 'frame.shp')
    pp_bisstra = os.path.join(pp['data_out_bast'], 'BFStr_Netz_SKPolyline.shp')
    pp_bisstra_nk = os.path.join(pp['data_out_bast'], 'BFStr_Netz_NKPoint.shp')
    pp_temp = os.path.join(pp['data_temp'], 'BISSTra_lage.shp')

    pp_out_pts = os.path.join(pp['data_out_zeb'], "zeb_pts.shp")
    pp_out_csv = os.path.join(pp['data_out_df'], "zeb.csv")
    
    return(zeb_files, pp_frame, pp_bisstra, pp_bisstra_nk, pp_temp, pp_out_pts, 
        pp_out_csv)


