# -*- coding: utf-8 -*-
"""

"""
import pandas as pd
import numpy as np
import re



def add_differentiated_sw_type_vars(lanes_long):    
    '''
    '''
    sw_vars = []
    for direc in ['Entry', 'Exit']:
        for typ in ['main', 'sec']:
            lanes_long[typ+'_'+direc] = pd.get_dummies(
                                            lanes_long['sw_type']
                                        )[typ]*lanes_long[direc]
            sw_vars.append(typ+'_'+direc)
    print('More than one type of slipway assigned:',
          sum(lanes_long[sw_vars].sum(axis=1)>1))
    print('Exactly one type of slipway assigned:',
          sum(lanes_long[sw_vars].sum(axis=1)==1))
    
    return(lanes_long, sw_vars)
    

def to_wide_and_sum_over_values(df):

    cols = list(df.columns)
    df = long_to_wide(df, cols)
    # Sum up all columns of the same value
    for col in cols:
        sum_cols = [x for x in df.columns if x.startswith(col+'_')]
        df[col] = df[sum_cols].sum(axis=1)
        df.drop(sum_cols, axis='columns', inplace=True)
        
    return(df)


def values_to_columns(df, var, values, delete_old_cols=True):
    
    if ('length' not in df.columns or 'check_length' not in df.columns):   
        print(''''length' and 'check_length' must be in data frame columns.''')
    for val in values:
        # Make a new column for each unique value.
        df[var + ':' + str(val)] = np.repeat(0,len(df))
        df.loc[df[var]==val,var + ':' + str(val)] = \
            df.loc[df[var]==val, 'length']/df.loc[df[var]==val, 'check_length']
            
    if delete_old_cols == True:
        df.drop([var, 'length', 'check_length'], axis='columns', inplace=True)

    return(df)


def get_unique_vals(df, var):

    unique_vals = df[var].unique()
    a = []
    for i in unique_vals:
        if np.issubdtype(type(i), np.number):
            if int(i) == i:
                a.append(int(i))
            else:
                a.append(i)
        else:
            a.append(i)
            
    return(a)


def get_collists(df):
    
    rep_col = list(set(['_'.join(y.split('_')[:-1]) for y in [
                x for x in list(df.columns) if x[-1].isdigit()]]))
    start_col = [x for x in list(df.columns) if not x[-1].isdigit()]
    collapse_col = [x for x in rep_col if (x in start_col and x!='length')]
    unique_col = [x for x in start_col if x not in collapse_col]
    
    return(rep_col, start_col, collapse_col, unique_col)


def change_df_by_pattern(df, var, area, new_var, id_dict, change_var=True,
                         return_form='wide'):
    '''Take a variable name var, search for changes in this varriable that are
    in place only for a length between the two values defined in area, 
    and then return to the previous value. For these segments, remove the 
    temporary change in var and create a new binary variable new_var, that is 1
    for the temporary change and 0 otherwise. 
    Use for example to identify small bridges instead of real bridges, or 
    junctions instead of additional lanes for a short while.
    Return the changed dataframe.
    
    '''

    # Convert df to a long df, and drop complete nan columns.
    df['id'] = df.index
    df_long = pd.wide_to_long(df, [var+'_', 'length_'], i='id', j='num')
    df_long.columns=[x[:-1] if x[-1]=='_' else x for x in df_long.columns]
    df_long.dropna(axis=0, how='any', inplace=True)
    # Sort df 
    df_long.sort_index(inplace=True)
    
    # Initialize new variable. 
    df_long[new_var] = np.repeat(0, len(df_long))
    
    
    for key in id_dict.keys():
        
        # Select by ids. 
        df_sub = df_long.loc[id_dict[key],:]
        a = df_sub.iloc[0]['length']
        b = df_sub.iloc[0][var]
        ids = [df_sub.index[0]]
        
        for i in range(1,len(df_sub)):
            if df_sub.iloc[i][var] == b:
                a += df_sub.iloc[i]['length']
                ids.append(df_sub.index[i])
            else:
                # evalutate here.... 
                if (area[0]  <= a <= area[1]):
                    if b == 'missing':
                        pass
                    elif df_sub.iloc[i][var] == b - 1:
                        df_sub.loc[ids, new_var] = 1
                        if change_var == True:
                            df_sub.loc[ids, var] = df_sub.loc[ids, var] - 1
                    
                a = df_sub.iloc[i]['length']
                b = df_sub.iloc[i][var]
                ids = [df_sub.index[i]]
                
        df_long.loc[id_dict[key],:] = df_sub
    
    if return_form == 'wide':
        df_wide = long_to_wide(df_long, [var, new_var, 'length'])
        return(df_wide)
    elif return_form == 'long':
        return(df_long)
    else:
        print('''form must either be 'long' or 'wide' to get return value.''')


def select_and_wide_to_long(df, var, rename=True):
    '''
    '''
        
    df_sel = select_and_collapse_columns(df, var, rename=True)
    df_sel = df_sel.join(df['check_length'])
    df_sel['id'] = df_sel.index
    
    if type(var) == list:
        varlist = [x+'_' for x in var]
    else:
        varlist =  [var+'_']
    
    df_sel = pd.wide_to_long(df_sel, varlist+['length_'], i='id', j='num')
    df_sel.columns=[x[:-1] if x[-1]=='_' else x for x in df_sel.columns]
    df_sel.dropna(axis=0, how='any', inplace=True)
    
    # Sort by ID.
    df_sel.sort_index(inplace=True)

    return(df_sel)


def long_to_wide(df, stubnames):    
    
    # Group by num and change into columns. 
    df['idx'] = df.index.get_level_values('num')
    df_wide = df.pivot_table(index=['id'], columns='idx', 
                             values=stubnames, aggfunc='first')
    df_wide.columns = [f'{x}_{y}' for x,y in df_wide.columns]
    # Sort columns like they were before.

    return(df_wide)


def select_and_collapse_columns(df, var_name_root, rename=False):  
    '''Inputs:
        * df: pandas data frame in wide format. 
        * var_name_root: can be a list or a string of variable name roots.
            'length' should not be part of it but will be added by default.
        * rename: bool, if True, the variables will be renamed so that the 
            first variable has the index 0, and then count upwards.
    '''
    
    if type(var_name_root) != list:
        var_name_root = [var_name_root]
        
    col_slice=[]
    for var in var_name_root:
        col_slice.extend(get_colslice(df, var))
        
    
    col_lengthes = get_colslice(df, 'length')
    
    # collapse df for the individual column.
    df_col = df[col_slice+col_lengthes]
    df_col = df_col.apply(collapse_row, axis=1, d_type = 'numeric')
    # Remove columns that contain only  nans.
    df_col.dropna(axis=1, how='all', inplace=True)
    
    # Rename the columns to keep the right ordering.
    
    if rename==True:
        new_colnames = {'length':'length_0'}
        
        for var in var_name_root:
            new_colnames.update({var:var +'_0'})
            
        for col in df_col.columns:
            if col in col_slice[1:]+col_lengthes[1:]:
                new = '_'.join(col.split('_')[0:-1]) + '_' + \
                       str(int(col.split('_')[-1]) + 1)
                new_colnames.update({col: new})

        df_col.rename(columns=new_colnames, inplace=True)
    
    return(df_col)


def get_colslice(df, var_name_root):
    
    cols = list(df.columns)
    # Make the pattern the string must match
    pattern=re.compile("^{}_\d+$".format(var_name_root))
    
    if var_name_root in df.columns:
        col_slice = [var_name_root] + \
                    [x for x in cols if (x.startswith(var_name_root) 
                                     and pattern.match(x) != None)]
    else:
        col_slice = [x for x in cols if (x.startswith(var_name_root) 
                                     and pattern.match(x) != None)]
        
    return(col_slice)


def collapse_row(row, d_type):
    '''Take a row that consists of the columns from one variable and the 
    repective length values in the back (all ordered appropriately).
    Collapse the lengths if there are two times the same values in the variable
    columns and assign nan to all variable columns of the collapsed position.
    '''
    
    # Handle false 
    if d_type == 'numeric':
        for j, i in enumerate(row):
            if type(i) is str:
                if i.isnumeric():
                    row[j] = int(row[j])
                
    a = int(len(row)/2) 
    for i in range(a-1,0,-1):
        if pd.notna(row[i]):
            if row[i]==row[i-1]:
                row[i] = np.nan
                row[i+a-1] = row[i+a-1] + row[i+a]
                row[i+a] = np.nan
    repl_list = [x for x in row[0:a] if pd.notna(x)] 
    row[0:a] = repl_list + list(np.repeat(np.nan, (a-len(repl_list))))
    repl_list = [x for x in row[a:] if pd.notna(x)]
    row[a:] = repl_list + list(np.repeat(np.nan, (a-len(repl_list))))

    return(row)

   
def get_maxspeed_c_dict():
    mc_dict = {
             '100 @ (06:00-19:00)':'100 @ 6-9',
             '100 @ (06:00-19:00;Su,PH off)':'100 @ 6-9',
             '100 @ (2014 Sep 29 - 2014 Nov 14 00:00-24:00)':0,
             '100 @ (22:00-06:00)': '100 @ 22-6',
             '100 @ (Mo-Fr 06:00-10:00)': '100 @ Mo-Fr 6-10',
             '100 @ (Mo-Fr 07:00-09:00)': '100 @ Mo-Fr 7-9',
             '100 @ (wet)':'100 @ wet',
             '100 @ (signal)': 'yes',
             '100 @ signal': 'yes',
             '120 @ (06:00-20:00)': '120 @ 6-20',
             '120 @ (06:00-21:00)': '120 @ 6-21',
             '120 @ (21:00-05:00)': '120 @ 21-5',
             '120 @ (22:00-06:00)': '120 @ 22-6',
             '120 @ (Mo-Fr 07:00-09:00)': '120 @ Mo-Fr 7_9',
             '120 @ (Mo-Su 06:00-20:00)': '120 @ 6-20',
             '130 @ (06:00-19:00)': '130 @ 6-19',
             '130 @ (06:00-19:00);80 @ wet': 'yes',
             '130 @ (06:00-20:00)': '130 @ 6-20',
             '130@(wet)':'130 @ wet',
             '60 @ wet': 'yes',
             '60 @ hgv': 'yes',
             '60 @ truck safety control by BAG': 'yes',
             '80 @ (2014 May 01-2015 Jun 22 00:00-24:00)':0,
             '80 @ (2014 May 01-2015 June 22 00:00-24:00)':0,
             '80 @ (22:00-06:00)': '80 @ 22-6',
             '80 @ (May 09-Jul 23)':0,
             '80 @ weg': '80 @ wet',
             '80@(wet)': '80 @ wet',
             '80@wet': '80 @ wet',
             '40 @ wet': 'yes',
             'missing':0,
             'none':0,
             'none @ (20:00-06:00)': 'none @ 20-6',
             '120 @ (06:00-09:00)': '120 @ 6-9',
             '100@wet': '100 @ wet',
             '120@wet': '120 @ wet',
             '800@wet': '80 @ wet'
             }
    return(mc_dict)





