# -*- coding: utf-8 -*-
"""

"""

import pandas as pd
from scipy.stats import t
from scipy.stats import ttest_ind


def add_t_test(a, df_temp, speed_c, crash_c):
    gr_none = df_temp[df_temp['None']==1]
    for sl in speed_c[:-1]:
        gr_sl = df_temp[df_temp[sl]==1]
        for cr in crash_c:
            t_test = ttest_ind(gr_none[cr], gr_sl[cr], equal_var=False)
            print(sl, cr, ":  ", round(t_test.statistic, 3))
            a.loc[cr,sl] = str(a.loc[cr,sl]) + "  " +\
                               "\textit{ \scriptsize{("+direction(t_test.statistic) + \
                                   stars(t_test.statistic, 
                                               len(gr_none)+len(gr_sl))+")}} "
    return(a)

def direction(t_stat):
    if t_stat < 0:
        s = ">"
    elif t_stat > 0: 
        s = "<"
    else:
        s = "="
    return(s)

def stars(t_stat, N):
    t_stat = abs(t_stat)
    if t_stat >= t.ppf(0.995, N):
        s = "$^{***}$"
    elif t_stat >= t.ppf(0.975, N):
        s = "$^{**}$"
    elif t_stat >= t.ppf(0.95, N):
        s = "$^{*}$"
    else:
        s = ""
    
    return(s)

def aggregate_years(df, year_list, same_cols, var_cols_sum, var_cols_average,
                    round_val):
    '''Take a long dataframe and aggregate over the years specified in the 
    'year' columns. 
    
    
    Arguments:
    df -- a long-format dataframe with the column 'year' indicating the year 
          the index being set to the segment index.
    year_list -- a list of numerical years to aggregate over. 
    same_cols -- a list of columns that have to have the same values over the
                 aggregation years. Otherwise rows with this ID are dropped. 
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
        df_year_dict.update({year: df[df['year']==year]})
    
    cols = df.columns
    var_cols = var_cols_sum + var_cols_average
    invar_cols = [x for x in df.columns if x not in var_cols_sum
                                                  + var_cols_average]
    
    df_w = df_year_dict[year_list[0]].join(df_year_dict[year_list[1]],
                                           lsuffix='_'+str(year_list[0])[-2:],
                                           rsuffix='_'+str(year_list[1])[-2:],
                                           how='inner')
    for i in range(2, len(year_list)):
        df_w = df_w.join(df_year_dict[year_list[i]], how='inner')
        df_w.rename(columns=dict(zip(cols, [x+'_'+str(year_list[i])[-2:] 
                                                for x in cols])), inplace=True)
    comb_list = []
    for i, j in enumerate(year_list):
        for x in year_list[i+1:]:
            comb_list.append([j, x])
    
    for col in same_cols:
        comb = comb_list[0]
        df_w[col+'_aux'] = pd.DataFrame(
                round(df_w[col+'_'+str(comb[0])[-2:]],round_val) == 
                round(df_w[col+'_'+str(comb[1])[-2:]],round_val)
        )
        for comb in comb_list[1:]:
            df_w[col+'_aux'] = pd.DataFrame(
                    (round(df_w[col+'_'+str(comb[0])[-2:]],round_val) == 
                     round(df_w[col+'_'+str(comb[1])[-2:]],round_val)) &
                    (df_w[col+'_aux'])
            )
        df_w[col+'_aux'] = df_w[col+'_aux'].astype(int)
    
    # Only select rows for which same_cols havn't changed significantly
    df_w['_aux'] = df_w[[x for x in df_w.columns if 
                            x.endswith('_aux')]].sum(axis=1)
    df_w = df_w[df_w['_aux']==len(same_cols)]
    # Drop auxiliary columns again.
    df_w.drop(columns=[x for x in df_w.columns if x.endswith('_aux')], 
                       inplace=True)
    
    for col in invar_cols:
        df_w[col] = df_w[col+'_'+str(max(year_list))[-2:]]
    
    for col in var_cols_sum:
        df_w[col] = df_w[[x for x in df_w.columns if x.startswith(col)]].sum(axis=1)
    
    for col in var_cols_average:
        df_w[col] = df_w[[x for x in df_w.columns if x.startswith(col)]].mean(axis=1)
    
    df_w = df_w[var_cols + invar_cols]
    
    # Drop all columns for which there is no variation in the variable.
    a = df_w.describe()
    drop_cols = a.loc[:,a.loc['std',:]==0].columns
    df_w.drop(columns=drop_cols, inplace=True)
    print('The following columns were dropped due to lack of variation:',
          drop_cols.values)
    
    return(df_w)


def average_x_per_y_df(df, x_cols, y_cols, ges_col_x, col_name_y, round_to=4):
    '''Derive the average of x_cols for each col in y_cols. 
    ges_col is the total of x_cols if available. Otherwise it has to be
    derived first. 
    
    col_name_y is the overall name for the cols in y_cols.
    
    '''
    
    # Make a column containing the speed limit as a string.
    df[col_name_y] = df[y_cols][
        df[y_cols]==1
    ].stack().reset_index().set_index('level_0').drop(columns=0)
    
    round(pd.crosstab(df[ges_col_x], 
                      df[col_name_y]).apply(lambda r: r/r.sum(), axis=0)*100,0)

    
    b = pd.DataFrame(pd.crosstab(df[ges_col_x], df[col_name_y]).sum()).T
    b.index = ['anzahl']
    
    for col in x_cols:
        a = pd.crosstab(df[col], df[col_name_y])
        a[ges_col_x] = a.index
        a=pd.DataFrame(a[y_cols].multiply(a[ges_col_x],axis='index').sum()).T
        a.index=[col]
        b = b.append(a)
    
    # Derive average number of accidents on segments with respective speed limti. 
    c = pd.DataFrame(columns=y_cols)
    for col in x_cols:
        a = pd.DataFrame(round(b.loc[col]/b.loc['anzahl'],round_to)).T
        a.index=[col]
        c=c.append(a)
    return(c)






