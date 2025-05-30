"""

runtime: 
"""
import pandas as pd
import os 
import re
import json

from project_paths import get_project_paths
import src.data_analysis.functions.visualize_functions as visual
import src.data_analysis.functions.tables_functions as tables
from src.process_numdata.functions.concat_functions import load_csv_and_fix_id

pp = get_project_paths(laptop=False)


p_in = pp['data_out_a_grf']
p_out = pp['fin_data']

s = "" # "", "r_", "r_sl_", "r_ps_"
exps = ["a", "b", "ba"]

in_p_me = os.path.join(pp['fin_data'], 'df_max_ext_imputed.csv')
df_me = load_csv_and_fix_id(in_p_me)

in_p = os.path.join(pp['fin_data'], 'analysis_df.csv')
df = load_csv_and_fix_id(in_p)
crash_counts = ['total', 'fatal', 'severely_injured', 'lightly_injured']
df = df.join(df_me[crash_counts])

pp_frame = os.path.join(pp['data_out_km'], 'frame.shp')
cates_d = json.load(open(os.path.join(pp['lib'], 'cates_cutoffs.json')))

limits = ["maxspeed_100",  "maxspeed_120", "maxspeed_130"]
outcomes = ['total_rate', 'fatal_rate', 'severe_rate', 'light_rate']

# Make a dictionary of rates and of cates and derive CI bounds and semi 
# elasticities. 
# Also add info to the df.
cd, df = visual.make_data_dicts(exps, limits, outcomes, df, "cates", 
                                  cates_d, p_in, s)
ra, df = visual.make_data_dicts([x for x in exps if x in ['a','b', 'ba']], limits, outcomes, df, "rate_q", 
                                  cates_d, p_in, s)
#ra_cr, df = visual.make_data_dicts([x for x in exps if x in ['a','b', 'ba']], limits, outcomes, df, "rate_q_cr", 
#                                  cates_d, p_in, s)

# Get dictionary with saved crashes. 
sc = tables.make_save_crash_dict_setup_b(limits, outcomes, pp, df_me)

      
##  Load data from the amtliche Statistik and derive the estimated reduction of
#   crashes based on the estimated ATEs. 
saved_crashes_old, n_people = tables.get_official_stat_values(pp, exps, limits, 
                                                          outcomes,df_me,s,cd)
# Get mean and std by speed limit table
ft = tables.get_stats_by_lim(df, df_me, pp, limits)
tables.ft_to_table(ft, os.path.join(pp['out_tables'], 'stats_by_lim.tex'))
# Make plots to compare principal components of the unrestricted analysis 
# sample and the unrestricted excluded sample.
tables.plot_var_by_sl_oos(df, df_me, pp)


####################################
# Create tables. ## 
outcomes_s = [x.replace("_rate","") for x in outcomes]
outcomes_sc = []
round_val = 3 # Round value
for x in [i for i in outcomes_s if i!="fatal"]:
    outcomes_sc.extend([x+" high",x+" low"])


############ Transofrm SE for crash rates to SE for crash numbers. 

count_eff_d = tables.get_count_effect_dict(s, exps, limits, outcomes, p_in, df,
                                           cates_d, cd, crash_counts)



# Make tables of ATEs
r = 3 # round value
ate_col_names = [x+"100" for x in outcomes_s]+[x+"120" for x in outcomes_s] + \
                [x+"130" for x in outcomes_s]

for exp in exps:
    if (exp =="b" and s==""):   
        ate_df = pd.DataFrame(index=["limit", "cr", "ATE", r"\textit{se}", 
                                     r"$\bar{SE}$", r"$\bar{Y}$",
                                     "RATE",r"\textit{se} ", 
                                     "ATE$_{c}$", r"\textit{se$_{c}$}",
                                     r"$\bar{SE}_{c}$"," "],
                              columns=ate_col_names)
    elif (exp in ["a", 'ba'] or (exp=="b" and s=="r_")):
        ate_df = pd.DataFrame(index=["limit", "cr", "ATE", r"\textit{se}", 
                                     r"$\bar{SE}$", r"$\bar{Y}$", "RATE", 
                                     r"\textit{se} ", "ATE$_{c}$", 
                                     r"\textit{se$_{c}$}", r"$\bar{SE}_{c}$"," "],
                              columns=ate_col_names)   
    else:
        ate_df = pd.DataFrame(index=["limit", "cr", "ATE", r"\textit{se}", 
                                     r"$\bar{SE}$", r"$\bar{Y}$","ATE$_{c}$", 
                                     r"\textit{se$_{c}$}", r"$\bar{SE}_{c}$", " "],
                              columns=ate_col_names)           
    n_cols = len(ate_df) - 1
    # Add basic information on limit and outcome.
    ate_df.loc['cr',:] = outcomes_s*3 
    ate_df.loc['limit',:] = ["100"]*4 + ["120"]*4 +["130"]*4 

    
    for limit in limits:
        for outcome in outcomes:
            e, ey, lim, out = visual.get_short_names(exp, limit, outcome)
            N =  df['th_'+exp+'_'+lim+'_'+out].notnull().sum()
            oc = outcome.replace("_rate","")
            j = oc+lim
            
            v = cd[exp][limit][outcome]["ate"]
            ate_df.loc["ATE", j] = str(round(v.estimate, r))+visual.stars(v.estimate,v["std.err"],N)
            ate_df.loc[r"\textit{se}", j] = tables.get_se(v["std.err"],r)
            # ate_df.loc["t", oc+lim] = round(v['estimate']/v["std.err"], r)
            ate_df.loc[r"$\bar{SE}$", j] = r"\textbf{"+str(round(v["elasticities"]))+"$\%$}"
            ate_df.loc[r"$\bar{Y}$", j] = round(v["Y_mean"],r)
            if exp in ['a', 'b', 'ba']:
                rate = ra[exp][limit][outcome]
                # rate_cr = ra_cr[exp][limit][outcome]
                
                if rate.loc["std.err"].values[0] == 0:
                    ate_df.loc["RATE", j] = "-"
                    ate_df.loc[r"\textit{se} ", j]  = "-"
                else:
                    ate_df.loc["RATE", j] = str(round(float(rate.loc["estimate"].values[0]),r))+\
                        visual.stars(float(rate.loc["estimate"].values[0]),float(rate.loc["std.err"].values[0]),N)
                    ate_df.loc[r"\textit{se} ", j] = tables.get_se(rate.loc["std.err"].values[0],r)
                    # ate_df.loc["t ", oc+lim] = round(blp.Estimate[1]/blp["Std. Error"][1],r)
                
                # ate_df.loc["RATE$_{cr}$", oc+lim] = str(round(float(rate_cr.loc["estimate"].values[0]),r))+\
                #     visual.stars(float(rate_cr.loc["estimate"].values[0]),float(rate_cr.loc["std.err"].values[0]),N)
                # ate_df.loc[r"\textit{se}  ", oc+lim] = round(float(rate_cr.loc["std.err"].values[0]),r)
                
            # Add info for crash counts. 
            out_count = [x for x in crash_counts if x.startswith(outcome[:5])][0]
            cv = count_eff_d[exp][limit][out_count]["ate"]
            ate_df.loc["ATE$_{c}$", j] = str(round(cv.estimate, r))+visual.stars(cv.estimate,cv["std.err"],N)
            ate_df.loc[r"\textit{se$_{c}$}", j] = tables.get_se(cv["std.err"],r)
            ate_df.loc[ r"$\bar{SE}_{c}$", j] = r"\textbf{"+str(round(cv["elasticities"]))+"$\%$}"
                
            if outcome == outcomes[-1]:
                 ate_df.loc[" ", oc+lim] = "\multicolumn{"+str(n_cols)+"}{l}{N="+str(N)+r" \quad $\bar{W}$="+str(round(v["W_share"],r))+"}"
                 
    
    f = ate_df.T
    # Set double index. 
    f.set_index(['limit', 'cr'],inplace=True,drop=False)
    
    # Put the multicolumn argument where it belongs and clean up. 
    f.loc[(f[' '][('100','light')]+"{removeme}", ' '),:] = ['101.123456'] + [' ']*n_cols # weird numbers to make them replacable
    f.loc[(f[' '][('120','light')]+"{removeme}", ' '),:] = ['121.123456'] + [' ']*n_cols
    f.loc[(f[' '][('130','light')]+"{removeme}", ' '),:] = ['131.123456'] + [' ']*n_cols
    f.rename(columns={'limit':'limit_col'},inplace=True)
    f.sort_values('limit_col',inplace=True)
    f.drop(columns=[' ','limit_col', 'cr'], inplace=True)
    
    
    with pd.option_context("max_colwidth", 1000):
        a = f.to_latex()
    a = a.replace("101.123456", " ").replace("121.123456", " ").replace("131.123456", " ")
        
    a = a.replace(r"\textbackslash ", "\\").replace(r"\{","{").replace(r"\}", 
          "}").replace(r"\\%", " \%").replace("\\$", 
          "$").replace(r"\textasciicircum ", "^").replace('l'*int(n_cols), {"b":"llrrrrrrrrr", "a":"llrrrrrrr", "ab":"llrrrrrr", "ba":"llrrrrrrr"}[exp]
          ).replace("\\_", "_")
    a = a.replace("\n100", 
                  "\n \multirow{4}{*}{\\rotatebox[origin=c]{90}{100}} ").replace("\n120", 
                  "\n \multirow{4}{*}{\\rotatebox[origin=c]{90}{120}} ").replace("\n130", 
                  "\n \multirow{4}{*}{\\rotatebox[origin=c]{90}{130}} ")
    # Insert decimal commas.
    a = re.sub("[0-9](?=(?:[0-9]{3})+(?![0-9]))", "\g<0>,", a)
    
                
    # Insert a same number of zeros for all decimal numbers.     
    a = re.sub(r'(?<=\.)(\d+)(?=[\s|$])',lambda m:m.group(1)+'0'*(r-len(m.group(1))),a)
    a = a.replace("( ", "(").replace(" )", ")")
        
    a = re.sub("limit & cr[^>]+midrule","midrule", a)
    a = a.replace("midrule", "\\midrule")

    a = re.sub("{removeme}[^>]+?\\\\\\n","\\\\\\\\\n", a)
    a = a.replace("\\multirow", "\\hline \\multirow")
    
    with open(os.path.join(pp['out_tables'], s+'grf_ates_'+exp+'.tex'), 'w') as tf:
        tf.write(a)
 

cates_p = ['AADT', 'HT_share','ramps']
# Make latex tables of CATES. 
for exp in [x for x in exps if x in ['a','b']]:
    for limit in limits:
        lim = limit[-3:]
        out_p=os.path.join(pp['out_tables'],s+'grf_cates_'+exp+'_'+lim+'.tex')
        tables.make_cate_table(cates_p, exp, limit, outcomes, df, cd, 
                               round_val, out_path = out_p)
        
    
# Save out analysis dataframe with all individual estimators etc. 
df.to_csv(os.path.join(pp['fin_data'], s+"all_results_df.csv"))


# Next, print out how large the share of segments is that has an exact counting
# station value. 
tables.print_share_of_sec_with_counting_station(pp)









