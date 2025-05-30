"""
Pocess oob predictions and all results into nice tables and a final 
data frame for the grf estimation. 

Runtime: approx. 2 minutes
"""

import pandas as pd
import json
import os 
from src.process_numdata.functions.concat_functions import load_csv_and_fix_id
import re 

from project_paths import get_project_paths, get_params




pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))

pp_col_dict = os.path.join(pp['fin_data'], 'df_dict.json')
cd = json.load(open(pp_col_dict, encoding="utf8"))

in_path = os.path.join(pp['fin_data'], 'df_prep.csv')
in_path_oob = pp['data_out_a_oob']
out_path = os.path.join(pp['fin_data'], 'analysis_df.csv')

short_dict = {"maxspeed_100":"100", "maxspeed_120":"120", "maxspeed_130":"130", 
              "fatal_rate": "FR", "severe_rate":"SR",
              "light_rate":"LR", "total_rate":"T"}

df = load_csv_and_fix_id(in_path)
df.drop(columns=[x for x in df.columns if x.startswith("aux")], inplace=True)
print(df.shape)


# Concatenate required data into the analysis_data
for s in short_dict.values():
    if s.isnumeric():  
        path = os.path.join(in_path_oob, 'oob_predictions_'+s+'.csv')
        df_temp = load_csv_and_fix_id(path).drop(columns=["id"])
        df = pd.concat([df, df_temp], axis=1, join='outer')
        
        # path = os.path.join(in_path_oob, 'oob_predictions_'+s+'_ffs.csv')
        # df_temp = load_csv_and_fix_id(path).drop(columns=["id"])
        # df = pd.concat([df, df_temp], axis=1, join='outer')
        print(s, df.shape)
    else:
        path = os.path.join(in_path_oob, 'oob_predictions_'+s.lower()+'.csv')
        df_temp = load_csv_and_fix_id(path).drop(columns=["id"])
        df = pd.concat([df, df_temp], axis=1, join='outer')
        print(s, df.shape)

df.to_csv(out_path)


short_dict2 = short_dict.copy()
short_dict2.pop('total_rate')
# Make tables of the FFS experiments.
exps = ['a', 'b', 'c']
t_ms = pd.DataFrame(index=[x for x in short_dict2.values() if x.isnumeric()], 
                    columns=pd.MultiIndex.from_product([['Oob error', 'Validation error'], 
                                                        ["a","b","c"]]))
                    
t_cr = pd.DataFrame(index=[x for x in short_dict2.values() if not x.isnumeric()], 
                    columns=pd.MultiIndex.from_product([['Oob error', 'Validation error'], 
                                                        ["a","b","c"]]))
rv = 4
rv_cr = 2
for s in short_dict2.values():
    if s.isnumeric():  

        for i in exps:
            path = os.path.join(pp['data_out_a_exp'], 'metrics_'+i+'_ps_'+s+'.csv')
            df_temp = pd.read_csv(path)
            t_ms.loc[s, ('Oob error', i)] = df_temp.loc[df_temp['metric']=='mse',
                                                      'oob_mean'].values[0].round(rv)
            t_ms.loc[s, ('Validation error',i)] = df_temp.loc[df_temp['metric']=='mse',
                                                            'val_mean'].values[0].round(rv)
            path = os.path.join(pp['data_out_a_exp'], 'var_imp_'+i+'_ps_'+s+'.csv')
            df_temp = pd.read_csv(path).sort_values('avrg_var_imp', 
                                                     ascending=False)
            
            
    else: 

        for i in exps:
            path = os.path.join(pp['data_out_a_exp'], 'mse_'+i+'_me_'+s+'.csv')
            df_temp = pd.read_csv(path)
            t_cr.loc[s, ('Oob error', i)] = df_temp.loc[
                    df_temp['Unnamed: 0']=='avrg. oob MSE', 'x'].values[0].round(rv_cr)
            t_cr.loc[s,('Validation error',i)] = df_temp.loc[
                    df_temp['Unnamed: 0']=='avrg. Validation MSE', 'x'].values[0].round(rv_cr)
            
            path = os.path.join(pp['data_out_a_exp'], 'var_imp_'+i+'_me_'+s+'.csv')
            df_temp = pd.read_csv(path).sort_values('avrg_var_imp', 
                                                     ascending=False)
with pd.option_context("max_colwidth", 1000):
    a = t_ms.to_latex()
    b = t_cr.to_latex()

a = a.replace("{} &    ", " Setup: & ").replace('lllllll','lrrrrrr').replace('{l}','{c}')                                               
b = b.replace("{} &    ", " Setup: & ").replace('lllllll','lrrrrrr').replace('{l}','{c}')                                               
                    
b = b.replace("FR", "Fatal rate").replace("SR"," Severe rate").replace("LR"," Light rate")                                                      
  
a = re.sub(r'(?<=\.)(\d+)(?=\s)',lambda m:m.group(1)+'0'*(rv-len(m.group(1))),a)
b = re.sub(r'(?<=\.)(\d+)(?=\s)',lambda m:m.group(1)+'0'*(rv_cr-len(m.group(1))),b)
                    
                                    
with open(os.path.join(pp['out_tables'], 'prop_score_experiment.tex'), 'w') as tf:
    tf.write(a)

with open(os.path.join(pp['out_tables'], 'main_effect_experiment.tex'), 'w') as tf:
    tf.write(b)
    

