# -*- coding: utf-8 -*-
"""

"""

import pandas as pd
import numpy as np





def add_lags(df, df_lag, lag_cols, lag0): 
    '''Add spatial lags from df_lag to df.'''
    lag_dict = get_lag_dict(df)
    
    print("No direct lags available for ", 
          str(len([x for x in list(lag_dict.values()) if len(str(x))!=10])),
          "of the ", df.shape[0] , " segments")
    
    lag_id_df = pd.DataFrame.from_dict(lag_dict, orient='index',columns=['lag_id'])
    
    df = df.join(lag_id_df)
    df['id'] = df.index
    df.set_index('lag_id', inplace=True, drop=False)
    df.index.name = None
    
    df = df.join(df_lag, how='left', rsuffix='_lag')
    df.set_index('id', inplace=True)
    df.index.name = None
    
    # Impute missing lags. 
    # First where 'lag_id' is not available
    for col in lag0['col']:
        df[col+'_lag'].loc[df['lag_id'].isna()] = df[df['lag_id'].isna()][col]
    for val in [0, 1, 30]:
        for col in lag0[val]:
            df[col+'_lag'].loc[df['lag_id'].isna()] = val
    
    # When lags are just missing values, impute same value as in segment. 
    for col in lag_cols:
        df[col+'_lag'].loc[df[col+'_lag'].isna()] = df[df[col+'_lag'].isna()][col]
    
    df.drop(columns=['lag_id'], inplace=True)
    
    return(df)


def get_lag_dict(df):
    '''Get a dictionary with the lags of each id in df, if it exists.'''
    
    lag_dict = {}
    ids = list(df.index)
    for i in ids:
        if i[-4:] == '0000':
            lag = np.nan
        else:
            t = str(int(i)-1)
            lag = (10-len(t))*'0' + t  
        lag_dict.update({i:lag})
    
    return(lag_dict)


def aggregate_years(df, year_list, same_cols, var_cols_sum, var_cols_average,
                    change_tolerance, max_extent=False):
    '''Take a long dataframe and aggregate over the years specified in the 
    'year' columns. 
    
    
    Arguments:
    df -- a long-format dataframe with the column 'year' indicating the year 
          the index being set to the segment index.
    year_list -- a list of numerical years to aggregate over. 
    same_cols -- a list of columns that are allowed to have changed at most by 
                 change_tolerance values over the aggregation years. 
                 Otherwise rows are dropped. 
    var_cols_sum -- columns that vary over time that will be summed up.
    var_cols_average -- columns that vary over time that will be averaged. 
    round_val -- determines how close the values of same_cols have to be to 
                 to each other to keep the column. 
    
    All columns not in same_cols + var_cols_sum + var_cols_average are treated
    as time-constant and the value from the last year is taken without further
    checks.
    
    '''
    
    year_list.sort()
    
    df_year_dict = {}
    for year in year_list: 
        df_year_dict.update({year: df[df['year'].apply(int)==int(year)]})
    
    cols = df.columns
    var_cols = var_cols_sum + var_cols_average
    invar_cols = [x for x in df.columns if x not in var_cols]
    
    df_w = df_year_dict[year_list[0]].join(df_year_dict[year_list[1]],
                                           lsuffix='_'+str(year_list[0])[-2:],
                                           rsuffix='_'+str(year_list[1])[-2:],
                                           how='inner')
    for i in range(2, len(year_list)):
        df_w = df_w.join(df_year_dict[year_list[i]], how='inner')
        df_w.rename(columns=dict(zip(cols, [x+'_'+str(year_list[i])[-2:] 
                                                for x in cols])), inplace=True)
    # List of all possible year pairs, required for comparison.
    if max_extent==False:
        comb_list = []
        for i, j in enumerate(year_list):
            for x in year_list[i+1:]:
                comb_list.append([j, x])
        
        comb_list = [x for x in comb_list if not (x[0]==min(year_list)) &
                                                 (x[1]==max(year_list))]
        
        for col in same_cols:
            comb = comb_list[0]
            df_w[col+'_aux'] = abs(df_w[col+'_'+str(comb[0])[-2:]] - \
                                    df_w[col+'_'+str(comb[1])[-2:]])
            
            for comb in comb_list[1:]:
                df_w[col+'_aux'] = df_w[col+'_aux'] + \
                        abs(df_w[col+'_'+str(comb[0])[-2:]] - \
                                 df_w[col+'_'+str(comb[1])[-2:]]) 
                
            df_w[col+'_aux'] = (df_w[col+'_aux'] < change_tolerance) * 1
        
        # Only select rows for which same_cols havn't changed significantly
        df_w['_aux'] = df_w[[x for x in df_w.columns if 
                                x.endswith('_aux')]].sum(axis=1)
        df_w = df_w[df_w['_aux']==len(same_cols)]
        # Drop auxiliary columns again.
        df_w.drop(columns=[x for x in df_w.columns if x.endswith('_aux')], 
                           inplace=True)
    # Invar cols are all columns not in var_cols_sum or var_cols_average.
    for col in invar_cols:
        df_w[col] = df_w[col+'_'+str(max(year_list))[-2:]]
    
    for col in var_cols_sum:
        df_w[col] = df_w[[x for x in df_w.columns if x.startswith(col) & 
                          (len(x)==len(col)+3)]].sum(axis=1)
    
    for col in var_cols_average:
        df_w[col] = df_w[[x for x in df_w.columns if x.startswith(col) & 
                          (len(x)==len(col)+3)]].mean(axis=1)
    
    df_w = df_w[var_cols + invar_cols]
    
    # Drop all columns for which there is no variation in the variable.
    a = df_w.describe()
    drop_cols = a.loc[:,a.loc['std',:]==0].columns
    df_w.drop(columns=drop_cols, inplace=True)
    print('The following columns were dropped due to lack of variation:',
          drop_cols.values)
    
    return(df_w)


def concat_year_data(path_dict, year, param, max_extent=False): 
    
    pp_BL         = path_dict['BL']
    pp_elevation  = path_dict['elevation']
    pp_zeb_inkar  = path_dict['zeb_inkar']
    pp_sk_node_loc= path_dict['sk_node_loc'] 
    pp_geometry   = path_dict['geometry']
    pp_wind       = path_dict['wind']
    pp_weather    = path_dict['weather']
    pp_crash      = path_dict['crash']
    pp_traffic    = path_dict['traffic']
    pp_ot         = path_dict['ot']
    
    colname_dict = {}
    ###########################################################################
    # Start with BL data.
    df = process_from_path(pp_BL, 
                           rename_dict={'GEN_LAN':'BL', 'GEN_KRS':'Landkreis'},
                           get_dummies=['BL'])
    BL_list = [x.replace('ue', 'ü') for x in param['BL_'+year]]
    subset_dict={'exact': {'BL': BL_list}}
    if not max_extent:
        df = subset_rows(df, subset_dict, delete_cols=False)
    df = df[BL_list]
    df.rename(columns={'Nordrhein-Westfalen': 'NRW', 'Baden-Württemberg': 'BW',
                       'Rheinland-Pfalz': 'RP',      'Schleswig-Holstein':'SH',
                       'Sachsen-Anhalt':'SA'}, 
              inplace=True)
    colname_dict.update({'BL':list(df.columns)})
    
    ###########################################################################
    ### Elevation data ###
    df_ele = process_from_path(pp_elevation)
    # Recode down change for ease of interpretation
    df_ele['down_change'] = abs(df_ele['down_change'])
    df_ele['elev_change'] = abs(df_ele['up_change']) + abs(df_ele['down_change'])
    colname_dict.update({'elevation':list(df_ele.columns)})
    df = df.join(df_ele)
    
    ###########################################################################
    ### ZEB (road quality) data and sociodemographic data. ###
    df_zi = process_from_path(pp_zeb_inkar, deselect=['sk_kennung'])
    df_zi.rename(columns={'asphalt_share': 'asphalt'},inplace=True)
    df_zi.drop(columns=['bauw_mode'], inplace=True)
    # Substanzwert = substantial index, Gebrauchswert = performance index
    cols = [x for x in df_zi.columns if 'bauw' in x or 'geb' in x or 'bip' in x]
    newc = [x.replace('bauw','surf').replace('geb','perf').replace('bip','gdp') 
            for x in cols]
    df_zi.rename(columns=dict(zip(cols, newc)),inplace=True)
    colname_dict.update({'zeb_inkar':list(df_zi.columns)})
    df = df.join(df_zi)
    
    ###########################################################################
    ### Sector and node influence area. ###
    df_sk_node = process_from_path(pp_sk_node_loc, deselect=['sk_kennung'])
    colname_dict.update({'sk_node_loc':list(df_sk_node.columns)})
    df = df.join(df_sk_node)   
    
    ###########################################################################
    ### Geometry ###
    df_geom = process_from_path(pp_geometry)
    colname_dict.update({'geometry':list(df_geom.columns)})
    df = df.join(df_geom)
    
    ###########################################################################
    ### Wind ###  
    df_wind = process_from_path(pp_wind)
    colname_dict.update({'wind':list(df_wind.columns)})
    df = df.join(df_wind)
    
    ###########################################################################
    ### Weather ###
    df_weather = process_from_path(pp_weather)
    colname_dict.update({'weather':list(df_weather.columns)})
    df = df.join(df_weather)
    
    ###########################################################################
    ### Crash data. ###
    aggr = {'count_category': 
               {'crash': {'mapping': {'Getötete':1, 'Schwerverletzte':2, 
                                      'Leichtverletzte':3},
                          'columns': {'startswith': ['UKATEGORIE']}}},
            'sum': {'Insgesamt': {'exact': ['Getötete', 'Schwerverletzte', 
                                            'Leichtverletzte']}}} 
    df_crash = process_from_path(pp_crash, aggregate_dict=aggr, 
                                 select=['Getötete', 'Schwerverletzte', 
                                         'Leichtverletzte', 'Insgesamt'],
                                 fillnans={'Getötete':0, 'Schwerverletzte':0,
                                           'Leichtverletzte':0,'Insgesamt':0})
    df_crash.rename(columns={'Insgesamt': 'total',
                             'Leichtverletzte':'lightly_injured',
                             'Schwerverletzte':'severely_injured',
                             'Getötete':'fatal',}, 
                    inplace=True)
    colname_dict.update({'crash':list(df_crash.columns)})
    df = df.join(df_crash)
    
    ###########################################################################
    ### Traffic count data. ###    
    df_traffic = process_from_path(pp_traffic, 
                                   rename_dict={'TVKfzms':'AADT',
                                                'TVSVms':'AADT_HT',
                                                'fer':'holiday',
                                                'bSo':'sunday',
                                                'SV5Kfzms': 'MSV50' ,
                                                'bSV50ms': 'HT_share_MSV50',
                                                'Mt': 'AADT_day' ,
                                                'Mn': 'AADT_night'},
                                   select=['holiday', 'AADT_day', 'AADT_night', 
                                           'AADT', 'AADT_HT', 'MSV50', 
                                           'HT_share_MSV50', 'sunday'])
    # Heavy traffic share in data not for each direction, thus calculate here. 
    df_traffic['HT_share'] = df_traffic['AADT_HT']/df_traffic['AADT']
    # All these variables are on the cross section. 
    df_traffic['night_share'] = df_traffic['AADT_night']/(
                            df_traffic['AADT_day']+df_traffic['AADT_night'])
    colname_dict.update({'traffic':list(df_traffic.columns)})
    df = df.join(df_traffic)
    
    ###########################################################################
    ### Open street maps data. ###
    deselect_cols = {'startswith': ['width', 'maxspeed_r', 'note:'],
                     'exact': ['tunnel:0', 'bridge:0', 'small_bridge:0', 
                               'junction:0', 'length', 'maxspeed_w:missing', 
                               'Slipway_lanes:0']}
    df_ot = process_from_path(pp_ot, deselect=deselect_cols)
    # Aggregate over Slipway_lanes:2 and Slipway_lanes_3 
    df_ot['large_sw'] = df_ot['Slipway_lanes:2']+df_ot['Slipway_lanes:3']
    
    # Create a variable containing various conditional speed limits. 
    df_ot['ms_vari_aux'] = df_ot[['ms_vari:yes','ms_vari:peak_traffic', 
                                  'ms_vari:obstruction']].sum(axis=1)
    if 'maxspeed:variable' in df_ot.columns:
        df_ot['maxspeed_aux'] = df_ot[['maxspeed:variable', 'maxspeed:signals']].sum(axis=1)
        df_ot['maxspeed_vari'] =  _get_var_union(df_ot, 'ms_vari_aux','maxspeed_aux')
    else:
        df_ot['maxspeed_vari'] =  _get_var_union(df_ot, 'ms_vari_aux','maxspeed:signals')
    drop_cols = [x for x in df_ot.columns if x.startswith('ms_vari')]+[x for x in [
            'maxspeed:signals', 'maxspeed:variable', 'maxspeed_aux'] if x in df_ot.columns]
    
    # Create a variable for missing shoulder.
    df_ot['shoulder_aux'] = df_ot['shoulder:no']+df_ot['shoulder:none']
    var='ms_reason:Fehlender Seitenstreifen (120)'
    if var in df_ot.columns:
        df_ot['no_shoulder'] = _get_var_union(df_ot, 'shoulder_aux', var)
    else:
        df_ot['no_shoulder'] = df_ot['shoulder_aux']
    drop_cols.extend([x for x in df_ot.columns if x.startswith('shoulder')])
    
    # Create a variable for conditional speed limits. 
    cond_vars = [x for x in ['ms_cond:missing', 'ms_cond:none', 'ms_cond:100 @ damage', 
                             'ms_cond:100'] if x in df_ot.columns]
    df_ot['ms_cond_a1'] = 1 - df_ot[cond_vars+[x for x in df_ot.columns if 
                             (x.startswith('ms_cond:'))&('2014' in x)]
            ].sum(axis=1)
    df_ot['ms_cond_a2'] = df_ot[['ms_wet:80','ms_wet:100', 'ms_wet:120']].sum(axis=1)
    df_ot['maxspeed_cond'] = _get_var_union(df_ot, 'ms_cond_a1','ms_cond_a2')
    drop_cols.extend([x for x in df_ot.columns if x.startswith(('ms_cond','ms_wet'))])
    
    df_ot['ot_tra_aux'] = df_ot[[x for x in df_ot.columns if x.startswith('ot_tra:no')]].sum(axis=1)
    df_ot['ot_hgv_aux'] = df_ot[[x for x in df_ot.columns if x.startswith('ot_hgv:no')]].sum(axis=1)
    df_ot['overtaking_ht'] =  _get_var_union(df_ot, 'ot_tra_aux','ot_hgv_aux')
    drop_cols.extend([x for x in df_ot.columns if x.startswith('ot')])
    
    # Create a variable that inform about construction sites. 
    df_ot['constr_aux1'] = df_ot[['constructi:minor', 'constructi:motorway', 
                                 'constructi:yes', 'ms_reason:construction']].sum(axis=1)
    df_ot['constr_aux2'] = df_ot[['note_maxsp:Baustelle', 'note_maxsp:Baustelle auf Brücke']].sum(axis=1)
    df_ot['construction'] = _get_var_union(df_ot, 'constr_aux1','constr_aux2')
    
    if max_extent==False:
        # Drop segments with construction site. 
        df_ot = df_ot[df_ot['construction']==0].drop(columns=['construction'])
    drop_cols.extend([x for x in df_ot.columns if x.startswith(
            ('constructi:', 'constr_', 'ms_reason:','under_cons', 'lit', 'hgv', 
             'ms_note', 'note', 'Slipway_lanes', 'ms_hgv', 'temp'))])
    
    # Some additional columns. 
    df_ot['straight'] = df_ot['starttoend']/df_ot['tot_length']
    df_ot['bridge_total'] = df_ot['bridge:1'] + df_ot['small_bridge:1']

    df_ot.drop(columns=drop_cols,inplace=True)
    
    # Average number of lanes on segment. 
    df_ot['n_lanes'] = 0
    for col in [x for x in df_ot.columns if x.startswith('lanes')]:
        df_ot['n_lanes'] += int(col[-1])* df_ot[col]
        
    # Remove inconsistent segmetns. 
    if not max_extent:
        df_ot = df_ot[abs(df_ot['tot_length']-df_ot['check_length'])<0.01]
        # Remove segments with wrong length
        df_ot = df_ot[(490<df_ot['check_length'])&(df_ot['check_length']<510)]
        # Remove segments with inplausible starttoend value. 
        df_ot = df_ot[df_ot['starttoend']>400]
        # drop if one of the columns is > 0.1: and drop respective columns.
        deselrow_dic = {'maxspeed:':   ['missing', '10', '20', '30', '40',
                                         '50', '60', '70', '80', '90', '140',
                                         'fixme:höchster üblicher Wert', '999']}
        drop_cols = ['starttoend', 'tot_length', 'check_length']
        df_ot = deselect_rows(df_ot, deselrow_dic, drop_cols)
    if max_extent:
        # Only for lagged variables: take absolute values since there are some
        # errors in the ot data (fo 3-4 segments not in the main df).
        for col in df_ot.columns:
            df_ot[col] = abs(df_ot[col])
    # Rename some columns
    rename_dict = dict(zip(list(df_ot.columns), 
                           [x.replace(':','_') for x in df_ot.columns]))
    rename_dict.update({'bridge:1':'bridge', 'small_bridge:1':'small_bridge',
                          'tunnel:1':'tunnel'})
    df_ot.rename(columns=rename_dict, inplace=True)
    # Create some new variables.
    # Dummy for containing lane change. 
    df_ot['lane_change'] = 0
    for x in [x for x in df_ot.columns if x.startswith('lanes')]:
        df_ot.loc[(round(df_ot[x],2)!=1)&
                  (round(df_ot[x],2)!=0),'lane_change'] = 1
    # Dummy for containing a change in speed limit.
    df_ot['maxspeedchange'] = 0
    for x in [x for x in df_ot.columns if x.startswith('maxspeed_') & (x!='maxspeed_cond')]:
        df_ot.loc[(round(df_ot[x],2)!=1)&
                  (round(df_ot[x],2)!=0),'maxspeedchange'] = 1

    # Ot and master must be joined with method 'inner' 
    colname_dict.update({'ot':list(df_ot.columns)})
    if not max_extent:
        df = df_ot.join(df, how='inner')
    else:
        df = df.join(df_ot)
    
    return(df, colname_dict)

def _get_var_union(df_ot, var1, var2):
    df_ot['new_var'] = df_ot[var1]*(
                df_ot[var1].round(2)==df_ot[var2].round(2)
            ) + df_ot[[var1,var2]].apply(max,axis=1)*(
                df_ot[var1].round(2)!=df_ot[var2].round(2))  
    return(df_ot['new_var'])

def deselect_rows(df, deselect_dict, drop_cols=[], tolerance=0.1):
    
    '''Drop row if one of the columns is > tolerance:
        and drop respective columns'''

    for key in deselect_dict.keys():
        for val in deselect_dict[key]:
            if key+val in df.columns:
                df = df.loc[df[key+val] <= tolerance]
                drop_cols.append(key+val)
                
    drop_cols = [x for x in drop_cols if x in df.columns]
    df.drop(columns=drop_cols, inplace=True)
    
    return(df)


def subset_rows(df, subset_dict, delete_cols=False):
    '''Inputs:
        * df: the data frame to perform the operations on and return.
        * subset_dict: dictionary with information on the subsetting operations
            to perform. Possible keywords:
            * query: contains a subdict with variables as keys and the queries
                as values. 
                Only keep rows fow which the query is True for this variable.  
            * exact: contains a subdict with variables as keys and the values
                that will be allowed for these variables as a list.
                Only keep rows that have one of the values in the lists for 
                this variable.
        * delete_cols: If True, the variables that are used for subsetting will
            be deleted in the end. By default they will not be deleted.
            
    '''
    varlist = []
    
    if subset_dict != None:
        for key in subset_dict.keys():
            
            if key == 'exact':
                for var in subset_dict[key].keys():
                    
                    df = df[df[var].isin(subset_dict[key][var])]   
                    varlist.append(var)
                        
            elif key == 'query':
                for var in subset_dict[key].keys():
                    if var in df.columns:
                        df = robust_query(df, var, subset_dict[key][var])
                        varlist.append(var)
                
    if delete_cols == True:
        df.drop(list(set(varlist)), axis='columns', inplace=True)
    
    return(df)


def clean_string(input_string, reverse=False):
    
    if reverse == False:
        ouput_string = input_string.replace(':','____').replace('-','__minus_'
                                  ).replace('@','__at_').replace('(','__ob_'
                                  ).replace(')','__cb_').replace(' ','__s_'
                                  ).replace('.','__d_').replace(';','__sem_')
    if reverse == True:
        ouput_string = input_string.replace('____',':').replace('__minus_','-'
                                  ).replace('__at_','@').replace('__ob_','('
                                  ).replace('__cb_',')').replace('__s_',' '
                                  ).replace('__d_','.').replace('__sem_',';')
        
    return(ouput_string)


def robust_query(df, variable, condition):
    
    # Rename columns so that queries will be possible...
    df.columns = [clean_string(x) for x in list(df.columns)]
    variable = clean_string(variable)
    
    query = variable + ' ' + condition
    df = df.query(query)
    # Reverse column rename to get original names.
    df.columns = [clean_string(x, reverse=True) for x in list(df.columns)]
    
    return(df)


def get_colnames(df, col_dict):
    
    cols = []
    
    for key in col_dict.keys():
        if key == 'startswith':
            for va in col_dict[key]:
                cols.extend([x for x in df.columns if 
                             x.startswith(va)])
        elif key == 'endswith':
            for va in col_dict[key]:
                cols.extend([x for x in df.columns if
                             x.endswith(va)])
        elif key == 'exact': 
            cols.extend(col_dict[key])
        
        else:
            raise TypeError("Unknown column selection keyword: '{}'".format(key))
        
    return(cols)



def process_from_path(path, relabel_dict=None, rename_dict=None, 
                      aggregate_dict=None, get_dummies=None, 
                      select=None, deselect=None, fillnans=None):
    '''Inputs:
        * path: path to csv dataset. (required)
        * aggregate_dict: ...
        * rename: rename dictionary.
        * select:   list or dictionary of columns to select.   
                    (Alternative to deselect)
        * deselect: list or dictionary of columns to deselect. 
                    (Alternative to select)
        * subset_dict: dictionary for which variable which values will be 
                       accepted. All other rows will be dropped.
                -> Only argument that selects row... maybe move to other function...
        * get_dummies: list of columns to get dummies for. 
        
    '''
    
    df = load_csv_and_fix_id(path)

    if relabel_dict != None:
        for var in relabel_dict:
            pass
        
    if rename_dict != None:
        df.rename(columns=rename_dict, inplace=True)
        
    if aggregate_dict != None:
        for key in aggregate_dict.keys():
            if key == 'count_category':
                for label in aggregate_dict[key].keys():
                    
                    sub_dict = aggregate_dict[key][label]
                    # Get desired columns first
                    cols = get_colnames(df, sub_dict['columns'])
                    
                    for var in sub_dict['mapping']:
                        
                        # Count occasions of categories. 
                        df[var] = df.apply(count_kat, axis=1,
                                           args=(cols, 
                                                 sub_dict['mapping'][var]))
            elif key == 'sum':
                for var in aggregate_dict[key].keys():
                    
                    cols = get_colnames(df, aggregate_dict[key][var])
                    df[var] = df[cols].sum(axis=1)
            else:
                raise TypeError("Unknown aggregation keyword: '{}'".format(key))
    
    if get_dummies != None:
        for var in get_dummies:
            df = df.join(pd.get_dummies(df[var]))
        # df = df.drop(get_dummies, axis='columns')
    
    if select != None:
        if type(select) == list:
            select = {'exact': select} 
        select = get_colnames(df, select)
        select = [x for x in select if x in df.columns]

        df = df[select]
        
    if deselect != None:
        if type(deselect) == list:
    	    deselect = {'exact': deselect}
        deselect_ = get_colnames(df, deselect)    	
        select_ = [x for x in df.columns if x not in deselect_]
        # select is stronger than deselect: keep all in select by force.
        if select != None:
            select_ = list(set(select + select_)) 
        select_ = [x for x in select_ if x in df.columns]
        df = df[select_]

    if fillnans != None:
        for var in fillnans.keys():
           df.loc[:,[var]] = df[[var]].fillna(0)

    return(df)


def count_kat(row, columns, category):
    
    a=0
    for col in columns:
        if row[col] == category:
            a += 1
    return(a)


def load_csv_and_fix_id(path, index_col=0, encoding=None):
        
    df = pd.read_csv(path, index_col=index_col, encoding=encoding)
    df['new_index'] = ['0'*(10-len(str(int(x))))+str(int(x)) for x in df.index]
    df.set_index('new_index', inplace=True)
    df.index.name = None
    
    return(df)



