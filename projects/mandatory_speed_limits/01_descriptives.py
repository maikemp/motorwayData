'''Finish cleaning and reshaping data.

'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
import os 
import re
from src.process_numdata.functions.concat_functions import load_csv_and_fix_id
import projects.mandatory_speed_limits.functions.descriptives_functions as descr

from project_paths import get_project_paths, get_params


pp = get_project_paths(laptop=False)
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))

pp_col_dict = os.path.join(pp['fin_data'], 'df_dict.json')
cd = json.load(open(pp_col_dict, encoding="utf8"))



in_path = os.path.join(pp['fin_data'], 'df_prep.csv')
max_ext_path = os.path.join(pp['fin_data'], 'df_max_extent.csv')


pp_out = pp['out_tables']

df = load_csv_and_fix_id(in_path)
df_me = load_csv_and_fix_id(max_ext_path)
crash_types = ['lightly_injured','severely_injured','fatal', 'total']
crash_c = ['Light', 'Severe', 'Fatal', 'Total']

df = df.join(df_me[crash_types])
speed_cols = ['maxspeed_100', 'maxspeed_120', 'maxspeed_130', 'maxspeed_none']
speed_c    = ['100', '120', '130', 'None']

# Add crash rates.
df_me['Total rate'] = df_me['total']/df_me['AADT']*10000
df_me['Fatal rate'] = df_me['fatal']/df_me['AADT']*10000
df_me['Severe rate'] = df_me['severely_injured']/df_me['AADT']*10000
df_me['Light rate'] = df_me['lightly_injured']/df_me['AADT']*10000

df.rename(columns={'total_rate':'Total rate','fatal_rate':'Fatal rate',
                   'severe_rate':'Severe rate', 'light_rate':'Light rate'}, 
            inplace=True)
crashrate_c = ['Light rate', 'Severe rate', 'Fatal rate', 'Total rate']

df_me = df_me.loc[(df_me['NRW']!=1) & (df_me['Berlin']!=1)&
                  (df_me['Thüringen']!=1) & (df_me[cd['BL']].sum(axis=1)==1)]
d_dict = {'df':df,'me':df_me}
me_ms_cols = [x for x in df_me.columns if x.startswith('maxspeed_') if (x not in
              ['maxspeed_cond']+speed_cols)&(not x.endswith('_lag'))] + speed_c


# Make descriptives tables of crashes per speed limit.
rv = 3
for k, v in d_dict.items():
    
    df_temp = v.copy()
    df_temp.rename(columns=dict(zip(speed_cols, speed_c)), inplace=True)
    df_temp.rename(columns=dict(zip(crash_types, crash_c)), inplace=True)
    if k == 'df':
        sec_share = [str(round(x/sum(list(df_temp[speed_c].sum()))*100,2))+" %"  
                     for x in list(df_temp[speed_c].sum())]
        p = ''
    elif k == 'me':
        df_temp.loc[df_temp[speed_c].sum(axis=1)==1,[x for x in me_ms_cols if x not in speed_c]]=0
        sec_share = [str(round(x/sum(list(df_temp[me_ms_cols].sum()))*100,2))+" %"  
                     for x in list(df_temp[speed_c].sum())]
        p = '_max_extent'

    # Give average number of crashes per speed limit. 
    a = descr.average_x_per_y_df(df_temp, crash_c, speed_c, 
                                    'Total', 'maxspeed', round_to=rv)
    a.columns = [(x, "  ") for x in a.columns]
    for i in ["100", "120", "130"]:
        a[(i, "   ")] = [x.replace("+-", "-") for x in 
                           ["(+"+str(x).replace(".0", " %")+")" for x in 
                            list(round(100*(a[(i,"  ")] - 
                                            a[('None',"  ")])/
                                       a[('None',"  ")]))]]
    
    
    a.loc['Speed limit share'] = sec_share + [" "]*3
    
    a['Overall'] = round(df_temp[crash_c].mean(), rv)
    a.replace(np.nan, " ", inplace=True)
    # Add N
    n = [int(x) for x in round(df_temp[speed_c].sum())] 
    a.loc["N"] = n + [" "]*3 + [str(int(sum(n)))]
    # Add total number of crashes.
    crash_num = [int(x) for x in df_temp[df_temp[speed_c].sum(axis=1)>0.5][crash_c].sum()] + ["",""]
    a['Count'] = crash_num
    print(df_temp[df_temp[speed_c].sum(axis=1)>0.5][crash_c].shape)
    # t-test. 
    # a = descr.add_t_test(a, df_temp, speed_c, crash_c)
    a = a.reindex(['Total', 'Fatal','Severe', 'Light', 'Speed limit share', 'N'])
    a.columns= pd.MultiIndex.from_tuples([x if type(x)==tuple else (x, " ") 
                                          for x in a.columns], names=["Speed Limit", " "])
    a = a[sorted(a.columns)[:6]+[a.columns[3]]+list(a.columns[7:])]
    
    # Give average crash rates per speed limit. 
    b = descr.average_x_per_y_df(df_temp, crashrate_c, speed_c, 
                                    'Total', 'maxspeed', round_to=rv)
    #b.loc['Speed limit share'] = sec_share
    #b['All'] = round(df_temp[crashrate_c].mean(), rv)
    #b.replace(np.nan, " ", inplace=True)
    b.columns = [(x, "  ") for x in b.columns]
    for i in ["100", "120", "130"]:
        b[(i, "   ")] = [x.replace("+-", "-") for x in 
                           ["(+"+str(x).replace(".0", " %")+")" for x in 
                            list(round(100*(b[(i,"  ")] - 
                                            b[('None',"  ")])/
                                       b[('None',"  ")]))]]
    
    
    b.loc['Speed limit share'] = sec_share + [" "]*3
    b['Overall'] = round(df_temp[crashrate_c].mean(), rv)
    b.replace(np.nan, " ", inplace=True)
    # Add N
    n = [int(x) for x in round(df_temp[speed_c].sum())] 
    b.loc["N"] = n + [" "]*3 + [str(int(sum(n)))]
    # Add total number of crashes.
    crash_num = [int(x) for x in df_temp[df_temp[speed_c].sum(axis=1)>0.5][crashrate_c].sum()] + ["",""]
    b['Count'] = crash_num
    print(df_temp[df_temp[speed_c].sum(axis=1)>0.5][crashrate_c].shape)
    
    b = b.reindex(['Total rate', 'Fatal rate','Severe rate', 'Light rate', 'Speed limit share', 'N'])
    b.columns= pd.MultiIndex.from_tuples([x if type(x)==tuple else (x, " ") 
                                          for x in b.columns], names=["Speed Limit", " "])
    


    print('------------------------------------------------------')
    print("For", k)
    print(a)
    print(b)
    a = a.to_latex().replace('llllllllll', 'llllllll|ll'
                              ).replace('\nSpeed limit share',
                            '\n\hline \nSpeed limit share')
    b = b.to_latex().replace('llllll', 'lllll|l'
                              ).replace('\nSpeed limit share',
                            '\n\hline \nSpeed limit share').replace("\\$",
                            "$").replace(r"\textasciicircum ", 
                            "^").replace(r"\{", "{").replace(r"\}", 
                            "}").replace(r"\textbackslash ", "\\")
    # Add thousands separator.
    a = re.sub("[0-9](?=(?:[0-9]{3})+(?![0-9]))", "\g<0>,", a)
    b = re.sub("[0-9](?=(?:[0-9]{3})+(?![0-9]))", "\g<0>,", b)
    
    
    
    
    a = a.replace('{ll', '{lll').replace("Speed", "  & Speed").replace("Total", 
              "\multirow{4}{*}{\\rotatebox[origin=c]{90}{Crash severity}} & Total").replace(
              "Severe"," & Severe").replace("Fatal"," & Fatal").replace(
                  "Light"," & Light").replace("\\\\\nN", "\\\\\n & N").replace("{l}","{c}")
    a = a.replace("{sl}", " Speed limit ")
    
    b = b.replace('{ll', '{lll').replace("Speed", "  & Speed").replace("Total", 
              "\multirow{4}{*}{\\rotatebox[origin=c]{90}{Crash severity}} & Total").replace(
              "Severe"," & Severe").replace("Fatal"," & Fatal").replace(
                  "Light"," & Light").replace("\\\\\nN", "\\\\\n & N").replace("{l}","{c}")
    b = b.replace("{sl}", " Speed limit ")  
    
    a = re.sub(r'(?<=\.)(\d+)(?=\s[^\\])',lambda m:m.group(1)+'0'*(rv-len(m.group(1))),a)
    b = re.sub(r'(?<=\.)(\d+)(?=\s[^\\])',lambda m:m.group(1)+'0'*(rv-len(m.group(1))),b)
    
    with open(os.path.join(pp_out, 'crash_freq_by_sl'+p+'.tex'), 'w') as tf:
     tf.write(a)
    with open(os.path.join(pp_out, 'crash_rate_by_sl'+p+'.tex'), 'w') as tf:
     tf.write(b)


#__________________________________________________________________________####
#   Make descriptives tables.                                              ####
rename_dict ={"air_temperature_mean": "air_temp", 
              "precipGE10mm_days":    "precip10mm",
              "precipGE20mm_days":    "precip20mm",
              "precipGE30mm_days":    "precip30mm",
              "sunshine_duration":    "sunshine_dur",
              "maxspeed_cond":        "ms_cond",
              "maxspeed_change":      "ms_change",
              "snowcover_days":       "snowcov_days",
              'HT_share_MSV50':       'HT_share_MSV',
              "p_comp_0": "pc_0",     "p_comp_1": "pc_1",
              "p_comp_2": "pc_2",     "p_comp_3": "pc_3"}

var_descr = json.load(open(os.path.join(pp['lib'], 'variable_descriptions.json')))
# Clean variable names to only contain actual name.
var_list = [x.replace("$","").replace("^","").replace("{","").replace("}",
            "").replace("*","").replace("+","").replace("(","").replace(")",
            "").replace(" ","") 
            for x in var_descr.keys()]

get_cols = []
rename_dict_rev = dict(zip(list(rename_dict.values()),
                           list(rename_dict.keys())))
for x in var_list:
    if x in rename_dict_rev.keys():
        if rename_dict_rev[x] not in df.columns:
            get_cols.append(rename_dict_rev[x])
    elif x not in df.columns:
        get_cols.append(x)

df = df.join(df_me[get_cols])
    
## Throw out later when columns where actually renamed!
df = df.rename(columns=rename_dict)
# Sort for better auto-correlation derivation.
df = df.sort_index()

m = []
s = []
ac = []
for col in var_list:
    
    if col == 'location':
        m.append(" ")
        s.append(" ")
        ac.append(" ")
    else:
        if np.mean(df[col]) >10:
            m.append(int(np.mean(df[col])))
            s.append(int(np.std(df[col])))
        else:
            m.append(round(np.mean(df[col]),2))
            s.append(round(np.std(df[col]),2))
        ac.append(round(df[col].autocorr(),2))
    
  
# Load variables for descriptives table and add metrics derived before.     
descr = pd.DataFrame({"Variable name":list(var_descr.keys()),
                      "Description":  list(var_descr.values()),
                      "Mean":         m, 
                      "Std.":         s,
                      "AC":           ac})    
    
 
with pd.option_context("max_colwidth", 1000):
    a = descr.to_latex(index=False)

a = a.replace(r"\textbackslash ","\\").replace("\\$","$").replace(r"textdegree", 
             " \degree ").replace("\\degreeC", " \degree C").replace("{\\}", 
             "").replace(r"\textasciicircum ", "^").replace("textit", 
             "\textit").replace(r"\{","{").replace(r"\}","}").replace("\\C",
             "C").replace("cdot", "\cdot") 
a = a.replace("\\\textit", "\\textit")
a = re.sub("[0-9](?=(?:[0-9]{3})+(?![0-9]))", "\g<0>,", a)

# Introduce horizontal lines before variables defined in this list.
h_vars = ['maxspeed_100', 'AADT', 'asphalt', 'down_change', 'air_temp', 
          'pop_dens', 'location']
for v in h_vars:
    v = v.replace("_", "\\\\_")
    vm = v.replace("\\\\_","\\_")
    a = re.sub('\\n\s+{}\s'.format(v), '\n slash_hline \n {} '.format(vm), a)
    a = re.sub('\\n\s+{}\$'.format(v), '\n slash_hline \n {}$'.format(vm), a)
    a = a.replace('slash_hline', '\hline')
# Save out descriptives table. 
with open(os.path.join(pp_out, 'variable_description.tex'), 'w') as tf:
     tf.write(a)


colors = ['black', 'dimgray', 'darkgray', 'lightgray']
for var in [x for x in df.columns if not 
            (x.startswith(('aux','maxspeed','location')) | x.endswith('_lag'))]:
    x = [df[df['maxspeed_none']==1][var], df[df['maxspeed_130']==1][var],
         df[df['maxspeed_120']==1][var],  df[df['maxspeed_100']==1][var]]

    plt.hist(x, 15, density=True, histtype='bar', color=colors,
             label=["None", "130", "120", "100"])
    plt.xlabel(var, fontsize=14)
    plt.ylabel("Density", fontsize=14)
    plt.legend(prop={'size': 12})
    plt.title('Density of {} by speed limit.'.format(var.replace("pc_", "p_comp_")), fontsize=16)
    plt.savefig(os.path.join(pp['out_graphics'], var + "_by_sl"), dpi=400)
    plt.show()
    

    
    
    
    
    
    
    




