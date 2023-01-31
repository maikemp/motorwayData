# -*- coding: utf-8 -*-
"""

"""


import os, glob
from datetime import datetime
import pandas as pd
import numpy as np
import arcpy








def select_zrows_and_direction(df, del_cols):
    '''Select only frame segments that have a zaehlstelle assigned and keep 
    columns ending on R1 or R2 only for the direction the row corresponds 
    to. Also delete the columns defined in del_cols
    '''
    df_z = df.loc[df['DZNr'].notnull()]
    
    columns = list(df_z.columns)
    columns_r1 = [x for x in columns if x.endswith('R1')]
    columns_r2 = [x for x in columns if x.endswith('R2')]
    columns_join = [x[:-2] for x in columns_r1 if x[:-2]+'R2' in columns_r2]
    
    # Make a new clumn 'R1' that is one for all segments corresponding to R1
    # and zero otherwise. 
    df_z['R1'] = df_z.apply(assign_dir, axis=1)

    for col in columns_join:
        df_z[col] = df_z.apply(take_colums, axis=1, col=col)
        df_z = df_z.drop([col+'R1', col+'R2'], axis='columns')
    df_z = df_z.drop(del_cols, axis='columns')

    return(df_z)
    

def drop_double_stations(pp_temp_z, result_string, temp_path):
    # Drop counting stations with missing values if they are in the same 
    # Sector as another counting station.
    cs_df = pd.DataFrame.spatial.from_featureclass(pp_temp_z)
    d = list(cs_df.AbschnittA.value_counts()[
                                  cs_df.AbschnittA.value_counts()>1].index)
    a=cs_df[cs_df['AbschnittA'].isin(d)]
    droplist = []
    for j in d:
        b = len(a.loc[a['AbschnittA']==j].loc[a['TVKfzmsR1']!=-9999])
        if b == 1:
            droplist.append(a['DZNr'][a['AbschnittA']==j][
                                            a['TVKfzmsR1']==-9999].values[0])
        else:
            # If both nans or not-nans, just select one.
            droplist.append(a['DZNr'][a['AbschnittA']==j].values[0])
    # Some manual corrections due to sector that should be two sectors.
    droplist = [x for x in droplist if x not in [6730, 6831]]
    
    cs_df = cs_df[cs_df["DZNr"].isin(droplist) == False]
    # Neccessary operation to avoid spatial.to_featureclass kill the colnames.
    rename_dict = dict(zip([x for x in cs_df.columns if x!="SHAPE"],
                           [x.lower() for x in cs_df.columns if x!="SHAPE"]))
    cs_df = cs_df.rename(columns=rename_dict)
    cs_df.spatial.to_featureclass(result_string)
    
    # Renaming works only in a gdb...
    # Due to a bug, changing only cases of letters not working, therefore,
    # double renaming with another column name required in between.
    arcpy.CopyFeatures_management(result_string, temp_path)
    rename_dict_rand = dict(zip([x for x in cs_df.columns],
                           ['a'+x+'b' for x in cs_df.columns]))
    for fn in rename_dict_rand.keys():
        arcpy.AlterField_management(temp_path, fn, rename_dict_rand[fn])    
    rename_dict_rev = {'a'+v+'b': k for k, v in rename_dict.items() if v in 
                   [f.name for f in arcpy.ListFields(result_string)]}
    for fn in rename_dict_rev.keys():
        arcpy.AlterField_management(temp_path, fn, rename_dict_rev[fn], "")
    
    arcpy.CopyFeatures_management(temp_path, result_string)
    
    return()

def cs_to_id_df(pp_zaehlst, pp_frame_dir, match_radius, column_list):
    '''
    '''
    
    # Next, project zaehlstellen to frame split by one_direct
    arcpy.CopyFeatures_management(pp_zaehlst, r"in_memory\zaehlst")
    for direct in [0,1]:
        print(direct)
        arcpy.MakeFeatureLayer_management(pp_frame_dir, 
                r"in_memory\one_dir", "\"one_direct\" = '{}'".format(direct))
                
        arcpy.analysis.SpatialJoin(
                 r"in_memory\one_dir",  r"in_memory\zaehlst", 
                 r"in_memory\frame_z_t", "JOIN_ONE_TO_ONE", "KEEP_ALL", 
                 match_option="CLOSEST", search_radius=match_radius) 
                
        data_dict = {}
        with arcpy.da.SearchCursor(r"in_memory\frame_z_t", 
                    ['ID', 'compass', 'direction', 'one_direct']+column_list) as c:
            for r in c:
                data_dict.update({r[0]:r[1:]})
        
        frame_df = pd.DataFrame.from_dict(data_dict, orient='index', 
                        columns=['compass', 'direction','one_direct']+column_list)
        if direct == 0:
            frame_add = frame_df
        if direct == 1:
            df = frame_add.append(frame_df)
    return(df)
    

def add_dir_mean(
    path_in, compass, id_field, dissolve=True, make_layer=True,
    memory_name=None, save_out=True, path_out=None
):
    '''
    '''
    
    # Get the ids and the compass angle mean direction per ID in a dictionary. 
    # after dissolving by seg_ID to get predominant route direction.
    if dissolve == True:
        arcpy.Dissolve_management(path_in, r"in_memory\diss_frame", 
                                  dissolve_field=id_field)
    else:
        arcpy.CopyFeatures_management(path_in, r"in_memory\diss_frame")
    # Create the directional mean layer.         
    arcpy.stats.DirectionalMean(r"in_memory\diss_frame", r"in_memory\dir",
                                    "DIRECTION", id_field)
    arcpy.Delete_management(r"in_memory\diss_frame")
    # Make a dictionary 
    dir_dict = {}
    with arcpy.da.SearchCursor(r"in_memory\dir", [id_field, 'CompassA']) as c:
        for r in c:
            for key in compass.keys():
                if key[0] <= c[1] < key[1]:
                    hr = compass[key]
            dir_dict.update({c[0]:[c[1], hr]})
    # Delete layer with only mean direction to avoid confusion.
    arcpy.Delete_management(r"in_memory\dir")
    
    if make_layer==False:
        return(dir_dict)
    else:
        
        # make a new layer from the frame, and add compass infos to it.
        arcpy.CopyFeatures_management(path_in, memory_name)
        # Add new fields to it.
        arcpy.AddField_management(memory_name, 'compass', 'FLOAT')
        arcpy.AddField_management(memory_name, 
                                  'direction', 'TEXT', field_length=2)
        cursor = arcpy.da.UpdateCursor(memory_name, 
                                       [id_field, 'compass', 'direction'])
        for r in cursor:
            r[1] = dir_dict[r[0]][0]
            r[2] = dir_dict[r[0]][1]
            cursor.updateRow(r)
        
        if save_out==True:
            # Save out new frame file for inspection
            arcpy.CopyFeatures_management(memory_name, path_out)
            

def select_and_rename_columns(df, keep_col):
    '''Select only the columns in keep_col and shorten the column names as 
    as defiend here. return the subsetted dataframe.
    
    '''
    
    df = df[keep_col]
    newcollist = []
    for col in df.columns:
           newcol = col.replace('MSV50','SV5').replace('MSV30','SV3')
           newcol = newcol.replace('DTV', 'TV').replace('MobisSo','ms')
           newcol = newcol.replace('W_MobisFr','mf').replace('DiMiDo','dd')
           newcol = newcol.replace('_','').replace('Ri', 'R').replace('NZB','NB')
           newcol = newcol[:10]
           newcollist.append(newcol)
    df.columns = newcollist
    for x in [len(x) for x in newcollist]:
        if x>10:
            print('''WARNING: Not all column names have been shortened enough 
                  to fit into a shapefile. Will be shortened automatically
                  when transformed to shapefile!''')
    return(df)


def make_file_paths(year, pp):
    '''
    '''
    
    pp_in = os.path.join(pp['data_in_count'], 'Jawe'+year+'.csv')

    stamp = datetime.now()
    pp_temp_df = os.path.join(pp['data_temp'], 
                              'zdf'+year+'_'+str(stamp.timestamp())+'.csv')
    old_paths = glob.glob(os.path.join(pp['data_temp'],'zdf'+year+'*'))
    for filepath in old_paths:
        if filepath != pp_temp_df:
            os.remove(filepath)
    
    pp_temp_z = os.path.join(pp['data_temp'], 'zaehlstellen' + year + '.shp')
    
    pp_frame = os.path.join(pp['data_out_km'], 'frame.shp')
    pp_out = os.path.join(pp['data_out_df'], 'zahlst_dir_' + year + '.csv')
    
    return(pp_in, pp_temp_df, pp_temp_z, pp_frame, pp_out)


def assign_dir(row):
    out = -6
    if row.HiR1 in row.direction:
        out=1
    elif row.HiR2 in row.direction:
        out=0
    else:
        if row.direction=='S':
            if row.compass < 180:
                new_dir = 'O'
            else:
                new_dir = 'W'
        elif row.direction=='N':
            if row.compass < 45:
                new_dir = 'O'
            else:
                new_dir = 'W'
        elif row.direction=='O':
            if row.compass < 90:
                new_dir = 'N'
            else:
                new_dir = 'S'
        elif row.direction=='W':
            if row.compass < 270:
                new_dir = 'S'
            else:
                new_dir = 'N'
        if row.HiR1 in new_dir:
            out=1
        elif row.HiR2 in new_dir:
            out=0
    return(out)


def take_colums(row, col):
    if row.R1 == 1:
        out = row[col+'R1']
    else:
        out = row[col+'R2']
    return(out)


