"""
"""

import arcpy
import pandas as pd
import numpy as np
import ast

import src.process_geodata.functions.f00_general_functions as general
from src.process_numdata.functions.f02_collapse_functions import get_collists


def remove_double_startpt(full_data):

    full_data.loc[:,'check_length'] = full_data[
                       [x for x in full_data.columns if x.startswith('length')]
                    ].sum(axis=1)
    
    a = full_data.loc[(full_data['length']==full_data['length_0']) & 
                      (full_data['check_length']>505)]
    rep_col, start_col, collapse_col, unique_col = get_collists(a)
    
    a.drop([x+'_0' for x in rep_col], axis=1, inplace=True)
    a.loc[:,'check_length'] = a[
                  [x for x in a.columns if x.startswith('length')]].sum(axis=1)
    
    new_cols = {}
    nums = []
    for col in a.columns:
        split_col = col.split('_')
        if split_col[-1].isnumeric():
            nums.append(split_col[-1])
            new_num = (len(split_col[-1])-len(str(int(split_col[-1])-1)))*'0' \
                       + str(int(split_col[-1])-1)
            new_col = '_'.join(split_col[:-1]) + '_' + new_num
            new_cols.update({col: new_col})
            
    a.rename(columns=new_cols, inplace=True)
    new_num = '0'*(max([len(x) for x in nums]) \
                   - len(str(max([int(x) for x in nums]))))\
                     + str(max([int(x) for x in nums]))
    a = pd.concat([a,pd.DataFrame(columns=[x+'_'+new_num for x in rep_col])])

    # Then update full data, including nans.
    full_data.loc[a.index.values,:] = a
    
    print('Number of double startpoints fixed: ', str(len(a)))

    return(full_data)


def make_changept_df(pp_change_points, excl_fields_list, id_list, 
                     field_duplic_num):
    '''
    '''
    # Save all field names and the respective type in a list and a dictionary.
    changept_fields = general.make_fieldlist(
        pp_change_points, ['ID'] + excl_fields_list, 
        last_fields = ['dist'], make_type_dict = False
    )
    
    changept_cursor = arcpy.da.SearchCursor(
        pp_change_points, ['ID'] + changept_fields
    )
    
    changept_dict = {}
    for row in changept_cursor:
        if row[0] in id_list:
            if row[0] not in list(changept_dict.keys()):
                changept_dict.update({row[0]: {row[-1]:list(row[1:])}})
            else:
                changept_dict[row[0]].update({row[-1]: list(row[1:])})
    
    # Make dictionary of form {id: wide format values.}
    changept_dict_ord = {}
    for ind in list(changept_dict.keys()):
        #ind = '6434170000'
        a = []
        b = list(changept_dict[ind].keys())
        b.sort()
        for key in b:
            a.extend(changept_dict[ind][key])
        changept_dict_ord.update({ind: a})
    
    # Fill up dictionary with nans for missing values.
    for key in list(changept_dict_ord.keys()):
        changept_dict_ord[key] = changept_dict_ord[key] + list(
                np.repeat(
                    np.nan, 
                    len(changept_fields) * field_duplic_num \
                        - len(changept_dict_ord[key]))
            )
    # Create column names for changepoints...
    changept_cols = []
    for i in range(field_duplic_num):
        for field_name in changept_fields:
            #changept_cols.append('0'*(2-len(str(i))) + str(i) +'_'+ field_name)
            changept_cols.append(
                    field_name + '_' + '0'*(
                            len(str(field_duplic_num-1))-len(str(i))
                        ) + str(i))

    # ordered changept dict to dataframe.
    changept_df = pd.DataFrame.from_dict(
        changept_dict_ord, orient='index', columns = changept_cols
    )
    
    return(changept_df)

def make_cutpt_df(pp_cut_points, excl_fields_list, id_list):
    '''
    '''
    
    # Select fields for cutpts
    cutpt_fields = general.make_fieldlist(
        pp_cut_points, ['ID'] + excl_fields_list, make_type_dict = False
    )
    
    # Fill up the dictionary with values. 
    cutpt_cursor = arcpy.da.SearchCursor(pp_cut_points, ['ID'] + cutpt_fields)
    
    cutpt_dict = {}
    for row in cutpt_cursor:
        if row[0] in id_list:
            cutpt_dict.update({row[0]: list(row[1:])})
    # Covnert to a dataframe.
    cutpt_df = pd.DataFrame.from_dict(
        cutpt_dict, orient='index', columns=cutpt_fields
    )
    
    return(cutpt_df)



def duplic_val_num(path, field):
    '''Count how often each value appears in field in the feature class in path
    and return the maximum number of it.
    '''
    
    cursor = arcpy.da.SearchCursor(path, [field])
    mylist = []
    for row in cursor:
        mylist.append(row[0])
    return(mylist.count(max(set(mylist), key = mylist.count)))



def intersect_ids(id_path, feature_path):
    '''Take ids from frame in frame_path and intersect them with ids from 
    (text) file in id_path. Return the intersection and print a report.
    '''
    
    with open(id_path, "r") as f:
        id_list_buff = ast.literal_eval(f.read())
    
    # Collect all ids of cut points.
    id_list_add = []
    with arcpy.da.SearchCursor(feature_path, "ID") as cursor:
        for row in cursor:
            id_list_add.append(row[0])
    # Only keep the ids, that are in buffer and in cut points.    
    id_list = list(set(id_list_buff).intersection(set(id_list_add)))   
     
    print("Number of selected elements from id-file:",len(id_list_buff))
    print("Number of selected elements from feature path:",len(id_list_add))
    print("Number of elements in both:",len(id_list))

    return(id_list)



    
    