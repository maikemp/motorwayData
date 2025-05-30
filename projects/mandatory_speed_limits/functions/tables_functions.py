# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 16:42:42 2022

@author: metzpmrh
"""

import os
import pandas as pd
import numpy as np
import json
# import arcpy
import re
from src.process_numdata.functions.concat_functions import load_csv_and_fix_id
import src.data_analysis.functions.visualize_functions as visual
import matplotlib.pyplot as plt

def plot_var_by_sl_oos(df, df_me, pp):
    colors = ['black', 'lightgray']
    for var in [x for x in df.columns if 'p_comp' in x]:
        x = [df[df['maxspeed_none']==1][var], df_me[var]]
    
        plt.hist(x, 15, density=True, histtype='bar', color=colors,
                 label=["Unrestricted analysis sample", "Unrestricted excluded sample"])
        plt.xlabel(var, fontsize=14)
        plt.ylabel("Density", fontsize=14)
        plt.legend(prop={'size': 12})
        plt.title('Density of {} by speed limit.'.format(var.replace("pc_", "p_comp_")), fontsize=16)
        plt.savefig(os.path.join(pp['out_graphics'], var + "_by_sl_oos"), dpi=400)
        plt.show()    
    return()

def print_share_of_sec_with_counting_station(pp):
    in_p = os.path.join(pp['fin_data'], 'analysis_df.csv')
    df = load_csv_and_fix_id(in_p)
    df_1 = load_csv_and_fix_id(os.path.join(pp['data_out_df'],
                                            'counts_interpolated2019.csv'), 
                               encoding='latin1')
    df = df.join(df_1['Sector'])[['Sector']]
    
    
    for year in ["2017", "2018", "2019"]:
    	# Load data and merge information of sector to it.
    	df_zs = load_csv_and_fix_id(os.path.join(pp['data_out_df'], 
                                              'zahlst_dir_'+year+'.csv'), 
                                    encoding='latin1')
    	df_1 = load_csv_and_fix_id(os.path.join(pp['data_out_df'], 
                                             'counts_interpolated'+year+'.csv'), 
                                   encoding='latin1')
    	df_zs = df_zs.join(df_1['Sector'])
    	df_zs = df_zs[['DZNr', 'AbschnittA', 'Hi', 'TVKfzms', 'Sector']]
    
    	# Only keep first row per sector and make sure it contains non-nan 
        # TVKfzms if available.
    	df_s = df_zs.sort_values(["Sector","TVKfzms"])
    	df_s = df_s.sort_values(by=['Sector', 'TVKfzms'], na_position='first')
    
    	df_f = df_s.groupby('Sector').first()
    	n_group = df_zs.groupby('Sector').size().reset_index(name='n_group')
    	df_f = df_f.merge(n_group, on='Sector', how='left')
    
    	df_f['has_val'] = df_f['TVKfzms'].notnull()*1
    
    	print(year, ":", 
    		  round(sum(df_f['has_val']*df_f['n_group']) / sum(df_f['n_group']) * 100,1),
    		  "% have an exact value of the counting station in full sample.")
    
    	df_2 = pd.merge(df, df_f[['Sector', 'has_val']], on='Sector', how='left')
    	n_group = df_2.groupby('Sector').size().reset_index(name='n_group')
    	df_2 = df_2.merge(n_group, on='Sector', how='left')
    	df_2 = df_2.sort_values(by=['Sector'], na_position='first')
    	df_2 = df_2.groupby('Sector').first()
    	print(year, ":", 
    	  round((df_2['has_val']*df_2['n_group']).sum() / df_2['n_group'].sum() * 100,1),
    	  "% have an exact value of the counting station in analysis sample.")
    	print("------------------------")
    return()

def get_se(se, r):
    
    s = str(round(float(se),r))
    nc = s.split(".")[1]
    s = s + (3-len(nc))*'0'
    
    return("("+s+")")

def get_count_effect_dict(s, exps, limits, outcomes, p_in, df, cates_d, cd, 
                          crash_counts):
    dd = dict(zip(exps,[dict(zip(limits,
                                 [dict(zip(crash_counts,[0]*len(crash_counts)))
                                  ]*len(limits)))]*len(exps)))
    for exp in exps:
        dd[exp] = {}
        for limit in limits:
            dd[exp][limit] = {}
            for outcome in outcomes:
        
                # print("-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- ")
                out_count = [x for x in crash_counts if x.startswith(outcome[:5])][0]
                
                e, ey, lim, out = visual.get_short_names(exp, limit, outcome)
                
                results_count = get_count_effects(p_in, df, s, exp, limit, 
                                                  outcome, out_count, cd, 
                                                  cates_d)
                dd[exp][limit][out_count] = results_count.copy()
    return(dd)

def get_count_effects(p_in, df, s, exp, limit, outcome, out_count, cd, cates_d):
    # Dataframe with individual results and all relevant information. 
    i_res, files_d = get_i_res(p_in, s, exp, limit, outcome, df,out_count)
    
    
    cates = cd[exp][limit][outcome]
    
    # Get AIPW
    aipw_cr = aipw(y=i_res[outcome], m=i_res["Y_hat"], w=i_res[limit], 
                   e=i_res["W_hat"], ite=i_res['tau_hat_oob'])
    #print(round(aipw_cr,6))
    
    i_res["DRs"] = get_DR_scores(y=i_res[outcome], m=i_res["Y_hat"], 
                                 w=i_res[limit], e=i_res["W_hat"], 
                                 ite=i_res['tau_hat_oob'])     
    std_err = get_std_err(i_res, aipw_cr, exp)
    #print(round(std_err,6))
    
    
    if round(aipw_cr,6) != round(cates.loc["estimate","ate"],6):
        raise ValueError("Something is wrong with manual AIPW calculation!")
    if round(std_err,6) != round(cates.loc["std.err","ate"],6):
        raise ValueError("Something is wrong with manual Std.Err calculation!")  
        
    # Transform tau_i to crash counts.
    i_res['test_tau'] = i_res['mu_hat_1']- i_res['mu_hat_0']
    if sum(abs(i_res['test_tau'] - i_res['tau_hat_oob']) < 0.0001) != len(i_res):
        raise ValueError("Something wrong with tau_hat derivation?!")
    
    i_res, files_d = get_i_res(p_in, s, exp, limit, outcome, df, out_count)
    # Transform m_0 and m_1 individually and from this get ite. 
    m_0 = (i_res["mu_hat_0"] * i_res["AADT"] * 3)/1000000
    m_1 = (i_res["mu_hat_1"] * i_res["AADT"] * 3)/1000000
    ite = m_1 - m_0
    
    # Transform m to crash counts
    m = (i_res["Y_hat"] * i_res["AADT"] * 3)/1000000
    
    i_res["DRs"] = get_DR_scores(y=i_res[out_count], m=m, w=i_res[limit], 
                                 e=i_res["W_hat"], ite=ite)     

    
    aipw_c = aipw(y=i_res[out_count], m=m, w=i_res[limit], e=i_res["W_hat"], 
                  ite=ite)
    se_c = get_std_err(i_res, aipw_c, exp)
    
    semi_elast_c = aipw_c/m_0.mean() * 100
    results_count = pd.DataFrame({'ate': [aipw_c, se_c, semi_elast_c, aipw_c/se_c]}, 
                                  index=['estimate', 'std.err', 'elasticities', "t"])
    
    # print(limit, "and outcome", outcome, ":", round(semi_elast_c,2),"% crash counts")
    # print("t = ", round(aipw_c/se_c, 2))
    
    for k,v in cates_d.items():
        if k in ["AADT", "HT_share", "ramps"]:
            i_res, files_d = get_i_res(p_in,s,exp,limit,outcome,df,out_count)
            i_res["DRs"] = get_DR_scores(
                y=i_res[out_count], m=m, w=i_res[limit], e=i_res["W_hat"], 
                ite=ite) 

            cates = get_subgroup_effects(i_res, var=k,cutoff=v,outcome=outcome, 
                                         out_count=out_count, limit=limit, 
                                         exp=exp)
            results_count = results_count.join(cates)
    
    return(results_count)


def get_subgroup_effects(i_res, var, cutoff, outcome, out_count, limit, exp):
        

    i_res_a_h = i_res[i_res[var]>cutoff].copy()
    m = (i_res_a_h["Y_hat"] * i_res_a_h["AADT"] * 3)/1000000
    m_0 = (i_res_a_h["mu_hat_0"] * i_res_a_h["AADT"] * 3)/1000000
    m_1 = (i_res_a_h["mu_hat_1"] * i_res_a_h["AADT"] * 3)/1000000
    ite = m_1 - m_0
    aipw_h = aipw(y=i_res_a_h[out_count], m=m, w=i_res_a_h[limit], 
                  e=i_res_a_h["W_hat"], ite=ite)
    se_h = get_std_err(i_res_a_h, aipw_h, exp)
    semi_elast_h = aipw_h/m_0.mean() * 100
    # print("--- HIGH",var,"---")
    # print(limit, "and outcome", outcome, ":", round(semi_elast_h,2),"% crash counts")
    # print("t = ", round(aipw_h/se_h, 2))
    
    i_res_a_l = i_res[i_res[var]<=cutoff].copy()
    m = (i_res_a_l["Y_hat"] * i_res_a_l["AADT"] * 3)/1000000
    m_0 = (i_res_a_l["mu_hat_0"] * i_res_a_l["AADT"] * 3)/1000000
    m_1 = (i_res_a_l["mu_hat_1"] * i_res_a_l["AADT"] * 3)/1000000
    ite = m_1 - m_0
    aipw_l = aipw(y=i_res_a_l[out_count], m=m, w=i_res_a_l[limit], 
                  e=i_res_a_l["W_hat"], ite=ite)
    se_l = get_std_err(i_res_a_l, aipw_l, exp)
    semi_elast_l = aipw_l/m_0.mean() * 100
    # print("--- LOW",var,"---")
    # print(limit, "and outcome", outcome, ":", round(semi_elast_l,2),"% crash counts")
    # print("t = ", round(aipw_l/se_l, 2))
    
    return(
        pd.DataFrame({'cate_l_'+var: [aipw_l, se_l, semi_elast_l, aipw_l/se_l], 
                      'cate_h_'+var: [aipw_h, se_h, semi_elast_h, aipw_h/se_h]}, 
                     index=['estimate', 'std.err', 'elasticities', 't'])
        )

def get_std_err(i_res, aipw_cr, exp):
    '''Get standard errors for AIPW of ATEs and GATEs. 
    
    Parameters
    ----------
    i_res : Data Frame of individual results
    aipw_cr:  cluster robust aipw
    exp : Setup, defining if clustered or not. 
    '''

    # De-mean DR scores.
    i_res["de_m"] = i_res["DRs"]-aipw_cr
    
    cl_df = pd.DataFrame({"tau_j": i_res.groupby("location").mean()["DRs"],
                          "n_j": i_res['location'].value_counts(),
                          "de_m": i_res.groupby("location").mean()["de_m"]})
    
    # aipw_cr also = sum(cl_df['tau_j']*cl_df['n_j'])/sum(cl_df['n_j'])
    J = len(cl_df)
    n = len(i_res)
    scale_J = J/(J-1)
    scale_n = n/(n-1)
    
    # Cluster robust standard errors.
    if exp in ["b", "ab"]:
        var = sum((cl_df["de_m"]*cl_df["n_j"])**2)/n**2 * scale_J
    # Non robust standard errors.
    elif exp in ["a", "ba"]:
        var = sum((i_res["de_m"])**2)/n**2 * scale_n
    else:
        raise ValueError("exp=", exp, "is not defined")
    
    std_err = np.sqrt(var)


    return(std_err)

def aipw(y, m, w, e, ite):

    DR_scores = get_DR_scores(y, m, w, e, ite)
    
    aipw = np.mean(DR_scores)
    
    return(aipw)

def get_DR_scores(y, m, w, e, ite):
    
    w_e = w - e 
    y_m =  y - m
    _1_e = 1 - e
    e1_ee = e * _1_e
    
    DR_scores = ite + ((w_e)/e1_ee)*(y_m - w_e * ite)
    
    return(DR_scores)
    

def get_i_res(p_in, s, exp, limit, outcome, df, out_count):
    
    e, ey, lim, out = visual.get_short_names(exp, limit, outcome)

    files_d = visual.get_files_dict(p_in, s, exp+'_', lim, out)
    # Merge all required files.
    i_res = load_csv_and_fix_id(files_d["i_res"])
    i_res = i_res.join(
        load_csv_and_fix_id(files_d["oob_main_eff"])["Y_hat_"+e+out]
        ).rename(columns={"Y_hat_"+e+out:"Y_hat"})
    i_res = i_res.join(
        load_csv_and_fix_id(files_d["oob_prop_score"])["W_hat_"+e+lim]
        ).rename(columns={"W_hat_"+e+lim:"W_hat"})
    i_res = i_res.join(df[[limit, outcome, out_count, "AADT", "location",
                           "HT_share", "ramps"]])
    
    return(i_res, files_d)

def ft_to_table(ft, out_path):
    ft = ft.round(3)
    # For large values, remove part after decimal point by converting to int and 
    # then to str. 
    ft.loc[ft.mean(axis=1)>=100] = ft.loc[ft.mean(axis=1)>=100].astype(int).astype(str)
    a = ft.to_latex()
    
    # Insert decimal commas.
    a = re.sub("[0-9](?=(?:[0-9]{3})+(?![0-9]))", "\g<0>,", a)
    a = re.sub(r'(?<=\.)(\d+)(?=[\s|$])',lambda m:m.group(1)+'0'*(3-len(m.group(1))),a)
    a = a.replace("lllllllllll", "lrrrrrrrrrr")
    with open(out_path, 'w') as tf:
        tf.write(a)
    return()


def get_stats_by_lim(df, df_me, pp, limits):
    with open(os.path.join(pp['lib'], 'variable_descriptions.json'), 'r') as f:
      var_descr = json.load(f)
    c =[x.replace("$^{*}$","").replace("$^{+}$","").replace("(",
                                                            "").replace(")","") 
        for x in var_descr.keys()]
    # cols = [x.replace(" ","") for x in c[:c.index('lightly_injured')+1] + \
    #         c[c.index('tunnel'):c.index('n_lanes')+1] + \
    #         c[c.index('AADT'): c.index('asphalt')+1] +['sub_max', 'perf_max']+\
    #         c[c.index('down_change'): c.index('max_slope')] ['pop_dens']]
    cols = [x for x in c if x in df.columns] + [x for x in df_me.columns if 'comp' in x]
    
    df = df.rename(columns={'HT_share_MSV50': 'HT_share_MSV'})
    df_me = df_me.rename(columns={'HT_share_MSV50': 'HT_share_MSV'})

    ft = df[df['maxspeed_none']==1][cols].mean()
    ft = pd.DataFrame(ft, columns=[('Maxspeed None', 'Mean')]).join(
        pd.DataFrame(df[df['maxspeed_none']==1][cols].std(), 
                     columns=[('Maxspeed None', 'Std.')]))
    for lim in limits: 
        l = lim.capitalize().replace('_',' ')
        ft = ft.join(pd.DataFrame(df[df[lim]==1][cols].mean(), columns=[(l,'Mean')]))
        ft = ft.join(pd.DataFrame(df[df[lim]==1][cols].std(), columns=[(l,'Std.')]))
    
    df_me_excl = df_me.join(df['holiday'], rsuffix='_')
    df_me_excl = df_me_excl[(df_me_excl['holiday_'].isna())]
    lim = 'maxspeed_none'
    me_cols = [x for x in cols if x not in ['total', 'fatal', 'severely_injured', 'lightly_injured']]
    ft = ft.join(pd.DataFrame(df_me_excl[df_me_excl[lim]==1][me_cols].mean(), columns=[("Maxspeed None excl.",'Mean')]))
    ft = ft.join(pd.DataFrame(df_me_excl[df_me_excl[lim]==1][me_cols].std(), columns=[("Maxspeed None excl.",'Std.')]))

    ft.columns = pd.MultiIndex.from_tuples(ft.columns, names=['Limit',' '])
    
    # Sort
    ft = ft[ft.columns[2:8].append(ft.columns[:2]).append(ft.columns[8:])]
    
    return(ft)


def get_official_stat_values(pp, exps, limits, outcomes, df_me, s, cd):
    
    pp_stat = os.path.join(pp['lit_cs'], 'Unfallstatistik.xlsx')
    cs = pd.read_excel(pp_stat)
    cs = cs[(cs['Unnamed: 1']=="Auf anderen Strecken")&
            (cs['Unnamed: 2']=="Ohne Geschwindigkeitsbegrenzung am Unfallort")]
    cs = cs[[x for x in cs.columns if not x.startswith("Unnamed")]]           
    
    # Add predictions to max_ext data. 
    pp_in_oob_me = os.path.join(pp["data_out_a_oob"], "max_ext")
    df_me = add_predictions_to_me(exps, limits, outcomes, df_me, pp, s)
    
    # Load ate for test data.
    path_ate_trans = os.path.join(pp['data_out_a_grf'], s+"ate_transport")
    
    saved_crashes = get_nr_annual_saved_crashes(['b'], limits, outcomes, 
                                                       df_me, path_ate_trans, s, 
                                                       pp, cs, cd, pp_in_oob_me)
    
    # To also derive number of injured/killed/lightly injured. 
    n_people = get_n_people(cs)   
    
    return(saved_crashes, n_people)

def make_save_crash_dict_setup_b(limits, outcomes, pp, df_me):
    
    sc_dict = {}
    
    exp='b'
    for limit in limits:
        sc_dict[limit] = {}
        for outcome in outcomes:
    #            if outcome in ['total_rate', 'light_rate']:
    #                rv = round_val - 1
            e, ey, lim, out = visual.get_short_names(exp, limit, outcome)
            
            a = df_me[['maxspeed_none','AADT']]
            a = a[a['maxspeed_none']==1]
            
            df_ir = load_csv_and_fix_id(os.path.join(pp['data_out_a_grf'], 
                                                     'ind_res_'+exp+'_'+lim+'_'+out+'.csv'))
            df_ir_me = load_csv_and_fix_id(os.path.join(pp['data_out_a_grf'], 
                                                        'ind_res_max_ext_'+exp+'_'+lim+'_'+out+'.csv'))
            a = a.join(df_ir['tau_hat_oob'])
            a = a.join(df_ir_me['tau_hat_i'])
            
            # Take tau_hat_oob if available, else take tau_hat_i. 
            a['tau'] = a['tau_hat_oob'].fillna(0) + 1*a['tau_hat_oob'].isna()*a['tau_hat_i']
            
            # Transform to counts
            a['tau_count'] = a['tau']*a['AADT']*3/1000000
            # Transform to annual values
            sc_dict[limit][outcome] = a['tau_count'].sum()/3
            
            print('--------------------------------------')
            print(limit)
            print(outcome)
            print(a['tau_count'].sum()/3) # Devide by 3 again to get annual values. 
            
    return(sc_dict)

  

def get_n_people(cs):
    
    dd = {}
    
    dd['fatal'] = cs['Getötete'].sum()/cs["mit Getöteten"].sum()
    dd["severe"] = cs['Schwerverletzte'].sum()/cs['mit Schwerverletzten'].sum()
    dd["light"] = cs['Leichtverletzte'].sum()/cs['mit Leichtverletzten'].sum()
   
    return(dd)
    

def add_predictions_to_me(exps, limits, outcomes, df_me, pp, s):
    for exp in exps:
        for limit in limits:
            for outcome in outcomes:
                lim, out, out_cr_nr, cs_col = get_short_names(limit, outcome)
                app = '_'+exp+'_'+lim+"_"+out
                
                df_temp = load_csv_and_fix_id(
                              os.path.join(pp["data_out_a_grf"],
                              s+"ind_res_max_ext"+app+".csv"))
                df_me = df_me.join(df_temp[['tau_hat_i', 'simga_hat_i']])
                # Rename according to variable combination.
                df_me = df_me.rename(columns={'tau_hat_i': 'tau_hat'+app, 
                                              'simga_hat_i': 'sigma_hat'+app})
        
                df_me['th_num'+app] = (df_me['tau_hat'+app]*df_me['AADT'])/10000
    return(df_me)




def get_short_names(limit, outcome):
    
    
    lim = limit.split("_")[1]

    if outcome == "total_rate":
        out_short = "t"
    else:
        out_short = "".join([x[:1] for x in outcome.split("_")])

    if outcome in ["total_rate", "fatal_rate"]:
        out_cr_nr = outcome.split("_")[0]
    else:
        out_cr_nr = outcome.split("_")[0]+'ly_injured'
     
    # Also define column of the official statistic. 
    if outcome=="total_rate":
        os_col = 'Unfälle mit Personenschaden insgesamt'
    elif outcome=="light_rate":
        os_col = 'mit Leichtverletzten'
    elif outcome=="severe_rate":
        os_col = 'mit Schwerverletzten'
    elif outcome=='fatal_rate':
        os_col = 'mit Getöteten'
    
    
    return(lim, out_short, out_cr_nr, os_col)



def get_nr_annual_saved_crashes(exps, limits, outcomes, df_me, path_ate_trans,
                                s, pp, cs, cd, pp_in_oob_me): 
    
    df_me_nsl = df_me[df_me['maxspeed_none']==1]

    # Make a dictionary for saved crahses. 
    sc = dict(zip(exps,[dict(zip(limits,
                                 [dict(zip(outcomes,[0]*len(outcomes)))
                                  ]*len(limits)))]*len(exps)))
    
    for exp in exps:
        sc[exp] = {}
        for limit in limits:
            sc[exp][limit] = {}
            for outcome in outcomes:   
                lim, out, out_cr_nr, cs_col = get_short_names(limit, outcome)
              
                ate_trans = pd.read_csv(
                                os.path.join(path_ate_trans,  
                                             "ate_t_"+exp+"_"+lim+"_"+out+".csv"),
                                index_col=0)
                
                i_res = load_csv_and_fix_id(
                            os.path.join(path_ate_trans,  
                                         "ate_t_i_res_"+exp+"_"+lim+"_"+out+".csv"))
                
                # Subset to only get obs that were not in analysis data.
                a_df = load_csv_and_fix_id(
                           os.path.join(pp['data_out_a_grf'], 
                                        s+"ind_res_"+exp+'_'+lim+'_'+out+'.csv')) 
                transp_ids = [i for i in df_me_nsl.index if i not in a_df.index]
                ateut_ids = [i for i in df_me_nsl.index if i in a_df.index]
                del(a_df)
                
                # collect data required to calculate elasticity for test data.
                i_res = i_res.loc[transp_ids]
                i_res_w_hat = load_csv_and_fix_id(
                                  os.path.join(pp_in_oob_me, 
                                               "max_ext_predictions_"+exp+'_'+lim+".csv")
                              ).rename(columns={"predictions": "W_hat_"+lim})
                i_res_y_hat = load_csv_and_fix_id(
                                  os.path.join(pp_in_oob_me, 
                                               "max_ext_predictions_"+exp+'_'+out+".csv")
                              ).rename(columns={"predictions": "Y_hat_"+out})
                i_res = i_res.join(i_res_w_hat).join(i_res_y_hat)
                
                i_res['mu_hat_0'] = i_res["Y_hat_"+out] - (i_res["W_hat_"+lim] * i_res["tau_hat_all"])
                elast = ate_trans["x"].loc["ATE_test"]/np.mean(i_res['mu_hat_0'])
                print("Elasticity of ATE_transport for", exp, limit, outcome,":", round(elast*100))
                print("t-test for", exp, limit, outcome,":", round(ate_trans.loc["ATE_test"]/ate_trans.loc["std_err"],3)[0])
                
                # First, get elasticity of ATEUT and multiply with observed crashes. 
                nc_train = df_me.loc[ateut_ids][out_cr_nr].sum()
                e_ateut = cd[exp][limit][outcome]["ateut"].loc["elasticities"]/100
                # Derive number of saved crashes. 
                saved_cr = e_ateut * nc_train + elast * (cs[cs_col].sum()-nc_train)
                
                sc[exp][limit][outcome] = saved_cr/3
                print("Number of saved crashes:", round(saved_cr/3))
                print("------------------------------------------------------------")
                
    return(sc)

def make_cate_table(cates_p, exp, limit, outcomes, df, cd, round_val, out_path):        
    cates_df = pd.DataFrame(index=pd.MultiIndex.from_tuples(
                   [(y,x) for y in cates_p for x in 
                       ["CATE",  r"\textit{se}", r"$\bar{SE}$", 
                        r"$\bar{Y}$", r"$\bar{W}$", "Diff."]], 
                names=[" ",""]),
                columns=pd.MultiIndex.from_product([
                        ['total', 'severe', 'light'], ['Low','High']]))
        
    for outcome in outcomes:    
        e, ey, lim, out = visual.get_short_names(exp, limit, outcome)

        if outcome != 'fatal_rate':
            for c in cates_p:
                # Collect values to test significance of difference.
                ests = []
                N = df['th_'+exp+'_'+lim+'_'+out].notnull().sum()
                for lev in ['Low','High']:
                    values = cd[exp][limit][outcome
                                ]["cate_"+lev.lower()[:1]+'_'+c]
                    oc = outcome.replace("_rate","")
                    rv = round_val
                    cates_df.loc[(c,"CATE"), 
                                 (oc, lev)] = str(
                            round(values.estimate, rv))+visual.stars(
                                values.estimate,values["std.err"],N)
                    cates_df.loc[(c, r"\textit{se}"), 
                                 (oc, lev)] = round(values["std.err"],rv)
                    cates_df.loc[(c,r"$\bar{SE}$"), 
                                 (oc, lev)] =  r"\textbf{"+str(
                                     round(values["elasticities"]))+"$\%$}"
                    cates_df.loc[(c,r"$\bar{Y}$"), 
                                 (oc, lev)] = round(values["Y_mean"],rv)
                    cates_df.loc[(c,r"$\bar{W}$"), 
                                 (oc, lev)] = round(values["W_share"],rv)
                    # Collect for difference test.
                    ests.append((values.estimate, values["std.err"]))
                
                diff_se = np.sqrt(ests[0][1]**2+ests[1][1]**2)
                diff = str(round(ests[0][0]-ests[1][0],rv)) + visual.stars(
                    ests[0][0]-ests[1][0], diff_se, N) +\
                    " (" + str(round(diff_se,rv))+")"
                cates_df.loc[(c,"Diff."),
                             (oc,lev)] = "\multicolumn{2}{c}{"+diff+"}"
        
        print("N =", N)
    
    with pd.option_context("max_colwidth", 1000):
        a = cates_df.to_latex()
    a = a.replace(r"\textbackslash ", "\\").replace("NaN &", 
                 "" ).replace(r"\{","{").replace(r"\}","}").replace(r"\\%",
                 " \%").replace("\\$","$").replace(r"\textasciicircum ", 
                 "^").replace("\\bottomrule",
                 "\\bottomrule \n \\multicolumn{8}{r}{*: p=0.1, **: p=0.05, ***: p=0.01}"
                 ).replace("{2}{l}", "{2}{c}")
    a = re.sub("[0-9](?=(?:[0-9]{3})+(?![0-9]))", "\g<0>,", a)
    for x in cates_p[1:]:
        xi = x.replace('_', '\\_')
        a=a.replace(xi,r"\hline "+xi)
        
    # Insert a same number of zeros for all decimal numbers.     
    a = re.sub(r'(?<=\.)(\d+)(?=[\s|$])',lambda m:m.group(1)+'0'*(rv-len(m.group(1))),a)
    # Insert a tiny space before a number starting with 0., if it is positive to align with negative numbers.
    # a = re.sub(r'\s(?=0\.)',lambda m: " \hspace{.01em} "+m.group(0), a)
    # # Insert little space if a semi elasticity has only one digit.
    # a = re.sub(r'(?<=-|{)(\d{1})(?=\$ )',lambda m: "\;"+m.group(0), a)      
    
    # Add more descriptions in the table. 
    a = re.sub('(?<=\\n)(\s+&\s+&\s+)(?=Low)', " & Subgroup & ", a)
    a = re.sub('(?<=\\n)(\s+&\s+&\s+)', " & Crash severity & ", a)
    a = a.replace("& {}",  " Variable & {}")
    
    with open(out_path, 'w') as tf:
        tf.write(a)

    return()


def make_spatial_df_all_res(df, pp_spatial_out, pp_frame, limits):
    frame_df = pd.DataFrame.spatial.from_featureclass(pp_frame)
    # Join maximum extent data to the  frame. 
    frame_df.set_index('ID', drop=False, inplace=True)
    
    # Add lines with 0 and 1 for W_hat to easier visualize the prop scores. 
    row1 = df.mean()
    row2 = df.mean()
    row1.loc
    row1[[x for x in df.columns if x.startswith("W_hat")]] = 1
    row2[[x for x in df.columns if x.startswith("W_hat")]] = 0
    dr = []
    for i in range(191160060,191160060+76):
        r1 = row1.copy()
        r2 = row2.copy()
        r1.name = '0' + str(int(i))
        r2.name = '0' + str(int(i+76))
        dr.extend([r1.name, r2.name])
        df = df.append(r1).append(r2)
         
    frame_df = frame_df.join(df, how="right")
    
    sl_cols = limits + ['maxspeed_none']
    # Create one variable for maxspeed for visualization.
    frame_df['maxspeed'] = pd.DataFrame(frame_df[sl_cols].idxmax(axis=1)).apply(
                                            lambda x: '_'.join(
                                                    str(x[0]).split('_')[1:]),
                                            axis=1)
    # Sort frame right    
    frame_df = frame_df[[x for x in frame_df.columns if x != 'SHAPE']+['SHAPE']]
    ms_cols = [x for x in frame_df.columns if 'maxspeed' in x]
    frame_df.rename(columns=dict(zip(ms_cols, 
                                     [x.replace('maxspeed','ms') for x in ms_cols])),
                    inplace=True)
    # 
    for x in [x for x in frame_df.columns if x!='SHAPE']:
        frame_df.fillna({x:-3.33333},inplace=True)
    
    frame_df.spatial.to_featureclass(r"in_memory\data")
    arcpy.management.CopyFeatures(r"in_memory\data", pp_spatial_out)
    df = df.drop(dr)
    return(df)


