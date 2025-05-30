# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 11:08:12 2022

@author: metzpmrh
"""


import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from adjustText import adjust_text

from scipy.stats import t
from src.process_numdata.functions.concat_functions import load_csv_and_fix_id

# import arcpy
# arcpy.env.overwriteOutput = True
# from arcgis import GIS


def plot_ateut_and_ateut_tranport(cd, s, exp, outcomes, limits, pp, color_list):

    fig, ax = plt.subplots(figsize=(10,4))
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.3)
    texts = []
    
    cats = ["ATEUT", "ate_transp"]

    for c in cats:
        va = c.lower()
        name = "the untreated in and out-of analysis sample."
        
                      
        if c in ['ate_transp']: # These are the darker colors
            #color_list = ['grey', 'darkgoldenrod', 'deepskyblue']  
            #['grey','indianred', 'goldenrod', 'skyblue'] 
            marker = "s"
            outcome_lab=[' ']*len(outcomes)
            v = 1.5
            lw = 1
        else: # Note: dardkgrey is lighter than grey!! 
            # color_list = ['darkgrey', 'burlywood', 'skyblue']
            #['darkgrey', 'pink', 'wheat', 'powderblue'] 
            outcome_lab=outcomes
            v = 0
            marker = 'd'
            lw = 0.5
    
        for i in range(len(outcomes)):
            plot_col = color_list[i]
            # d = df[df['outcome']==outcomes[i]]
            b = [x + i*4 + v for x in [0,20,40]] 
            
            elast = [cd[exp][l][outcomes[i]][va].elasticities for l in limits]
            me_upper = [cd[exp][l][outcomes[i]][va].ME_upper for l in limits]
            me_lower = [cd[exp][l][outcomes[i]][va].ME_lower for l in limits]
           
            ax.plot(b, elast, marker=marker, linestyle='', 
                    color=plot_col, label=outcome_lab[i])
            ax.plot(b, me_upper,marker='_',linestyle='',color=plot_col)
            ax.plot(b, me_lower,marker='_',linestyle='',color=plot_col)
            for j in b:
                idx = list(b).index(j)
                texts.append(plt.text(j, elast[idx]+2, str(round(elast[idx]))+"%", fontsize=12))
                # ax.annotate(str(round(elast[idx]))+"%", 
                #             (j, elast[idx]+2), fontsize=12)
                ax.plot([j,j],[me_lower[idx], me_upper[idx]],
                         color=plot_col, linewidth=lw)
    
    adjust_text(texts, only_move={'points':'y', 'text':'xy', 'objects':'xy'}, 
                arrowprops=dict(arrowstyle="->", color='silver', lw=0.5))
    
    ax.legend(loc='upper right', bbox_to_anchor=(1.1, 1.1),
              ncol=2, fancybox=True, shadow=True, 
              title=cats[0]+" "*12+"ATEUT$_{transp.}$")
    if s == 'r_':
        temp_str = ', restricted sample'
    else:
        temp_str = ''
    if exp=="b":
        ts = "Semi-elasticities and CIs for {}{} \n".format(name,
                                                            temp_str)
    else:
        ts = "Semi-elasticities and CIs for setup {}, {}{} \n".format(
             exp, name, temp_str)
    ax.set_title(ts, fontsize=14.5)
    # Make axis description without tickmarks. 
    ax.tick_params(axis='x', bottom=False, left=False) 
    plt.xticks([7.5, 27.5, 47.5], [100, 120, 130], size='large')
    ax.set_xlim(-2, max(b)+2+v)
    ax.set_ylim(-95, 50)
    #ax.get_xaxis().tick_bottom()
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    ax.set_ylabel("Semi-elasticity", fontsize=13)
    ax.set_xlabel("Speed limit", fontsize=13)

    fig.savefig(os.path.join(pp['out_graphics'], #'new',
                             s+'grf_CATES_ATEUT_'+exp+'.png'), 
                dpi=400)
    return()


def add_ate_transp_to_dict(cd, pp, limits, outcomes, s, df, df_me, exp):
    for l in limits:
        for o in outcomes:
            e, ey, lim, out = get_short_names(exp, l, o)
            cd[exp][l][o]['ate_transp'] = np.nan
            at = pd.read_csv(os.path.join(pp['data_out_a_grf'], s+'ate_transport', 
                                          'ate_t_'+exp+'_'+lim+'_'+out+'.csv'),
                             index_col=0)
            cd[exp][l][o]['ate_transp'].loc['estimate'] = at.loc['ATE_test'
                                                                 ].values[0]
            cd[exp][l][o]['ate_transp'].loc['std.err'] = at.loc['std_err'
                                                                ].values[0]
            cd[exp][l][o]['ate_transp'].loc['W_share'] = 0
            
            i_res = load_csv_and_fix_id(
                os.path.join(pp['data_out_a_oob'], 'max_ext', 
                             "max_ext_predictions_"+exp+'_'+out+".csv"))
            # Only get columns of i_res, that are not in df and that are unrestricted.
            i_res = i_res.join(df_me['maxspeed_none']).join(df['holiday'])
            i_res = i_res[(i_res['maxspeed_none']==1)&(i_res['holiday'].isna())]
            # Add individual tau_hat and W_hat
            i_tau = load_csv_and_fix_id(os.path.join(pp['data_out_a_grf'], 
                                        s+"ind_res_max_ext_"+exp+"_"+lim+"_"+out+".csv"))        
            i_e = load_csv_and_fix_id(os.path.join(pp['data_out_a_oob'], 'max_ext', 
                                                    "max_ext_predictions_"+exp+'_'+lim+".csv"))
            i_res = i_res.join(i_tau["tau_hat_i"]).join(
                            i_e.rename(columns={'predictions':'W_hat'}))
            del(i_tau,i_e)
            # mu_hat_0 <- Y.hat - W.hat * tau_hat_oob
            i_res['Y_hat_0'] = i_res['predictions'] - i_res['W_hat']*i_res['tau_hat_i']
            
            # Get elasticities.
            Y0_mean = i_res['Y_hat_0'].mean()
            elast = cd[exp][l][o]['ate_transp'].loc['estimate']/Y0_mean*100
            cd[exp][l][o]['ate_transp'].loc['elasticities'] = elast
            
            
            # Get CI bounds.
            half_CIwidth =  1.96*at.loc['std_err'].values[0]
            
            cd[exp][l][o]['ate_transp'
                          ].loc['CI_upper'
                                ] = at.loc['ATE_test'].values[0] + half_CIwidth
            cd[exp][l][o]['ate_transp'
                          ].loc['CI_lower'
                                ] = at.loc['ATE_test'].values[0] - half_CIwidth
            # Transform CI bounds to elasticity. 
            t_df = cd[exp][l][o]['ate_transp']
            cd[exp][l][o]['ate_transp'
                          ].loc['ME_upper'] = t_df.loc['CI_upper']/Y0_mean*100
            cd[exp][l][o]['ate_transp'
                          ].loc['ME_lower'] = t_df.loc['CI_lower']/Y0_mean*100
            cd[exp][l][o]['ate_transp'].loc['N'] = len(i_res['Y_hat_0'])
    
    return(cd)

def plot_ates(ate_list, cd, s, exps, outcomes, limits, pp, color_list):
    
    for ate in ate_list:
        for exp in exps:
            fig, ax = plt.subplots(figsize=(9,4.8))
            ax.axhline(y=0, color='black', linestyle='-', linewidth=0.3)
            texts = []
            for i in range(len(outcomes)):
                plot_col = color_list[i]
                
                b = [x + i*1.5 for x in [0,8,16]]
                print(b)
                #
                elast = [cd[exp][l][outcomes[i]][ate].elasticities for l in limits]
                me_upper = [cd[exp][l][outcomes[i]][ate].ME_upper for l in limits]
                me_lower = [cd[exp][l][outcomes[i]][ate].ME_lower for l in limits]
                ax.plot(b, elast, marker='o', linestyle='', color=plot_col, 
                        label=outcomes[i])
                ax.plot(b, me_upper, marker='_', linestyle='', color=plot_col)
                ax.plot(b, me_lower, marker='_', linestyle='',color=plot_col)  
                
                for j in b:
                    idx = list(b).index(j)
                    texts.append(plt.text(j, elast[idx]+2, str(round(elast[idx]))+"%", fontsize=12, ha='right'))
                    # ax.annotate(str(round(elast[idx]))+"%", 
                    #             (j, elast[idx]+2), fontsize=12)
                    ax.plot([j,j],[me_lower[idx], me_upper[idx]],
                             color=plot_col, linewidth=0.6)
            adjust_text(texts, only_move={'points':'y', 'text':'xy', 'objects':'xy'}, 
                        arrowprops=dict(arrowstyle="->", color='silver', lw=0.5))
            ax.legend(loc='upper left', bbox_to_anchor=(-0.1, 1.1),
                          ncol=2, fancybox=True, shadow=True)
            ax.set_xlim(-2, max(b)+2)
            ax.set_ylim(-70, 40)
            b = ""
            if ate == "ate_ov":
                b="for overlap \n"
            
            if exp=="b":
                if s == 'r_':
                    title_str = "ATE, restricted sample \n" + b
                else:
                    title_str = "ATE \n" + b
            else:
                if s == 'r_':
                    title_str = "ATE for setup "+exp+", restricted sample \n " + b
                else:
                    title_str = "ATE for setup "+exp+" \n" + b
                    
            ax.set_title(title_str, fontsize=15)
            # Make axis description without tickmarks. 
            ax.tick_params(axis='x', bottom=False, left=False) 
            plt.xticks([2.3, 10.4, 18.5], [100, 120, 130], size='large')
            # ax.get_xaxis().tick_bottom()
            ax.yaxis.set_major_formatter(mtick.PercentFormatter())
            ax.set_ylabel("Semi-elasticity", fontsize=14)
            ax.set_xlabel("Speed limit", fontsize=14)
            fig.savefig(os.path.join(pp['out_graphics'], #'new',
                                     s+'grf_'+ate+'_'+exp+'.png'), dpi=400)
    return()


def plot_cates(cates_d, cd, s, exps, outcomes, limits, pp, color_list):

    for exp in exps:
        print("---------------------------",exp,"---------------------------")
        for var in list(cates_d.keys())+["atet"]:
            print("-------------------------",var,"-------------------------")
            fig, ax = plt.subplots(figsize=(9,3))
            ax.axhline(y=0, color='black', linestyle='-', linewidth=0.3)
            texts = []
            if var == 'atet':
                cats = ["ATET", "ATEUT"]
                w = ""
                d = 31
            elif var == 'ramps':
                cats = ["Without","With"]
                w = var
                d = 18
            else:
                cats = ["Low","High"]
                w = var
                d = 24
            for c in cats:
                if var.startswith("ate"):
                    va = c.lower()
                    name = "untreated vs. treated"
                
                elif var == "ramps":
                    va = "cate_" + ["l" if c=="Without" else "h"][0] + "_" + var
                    name = "variable " + var
    
                else:
                    va = "cate_" + c[0].lower() + "_" + var
                    name = "variable " + var
                              
                if c in ['High', 'ATEUT', 'With']: # These are the darker colors
                    #color_list = ['grey', 'darkgoldenrod', 'deepskyblue']  
                    #['grey','indianred', 'goldenrod', 'skyblue'] 
                    marker = "s"
                    outcome_lab=[' ']*len(outcomes)
                    v = 1
                    lw = 1
                else: # Note: dardkgrey is lighter than grey!! 
                    # color_list = ['darkgrey', 'burlywood', 'skyblue']
                    #['darkgrey', 'pink', 'wheat', 'powderblue'] 
                    outcome_lab=outcomes
                    v = 0
                    marker = 'd'
                    lw = 0.5
            
                for i in range(len(outcomes)):
                    plot_col = color_list[i]
                    # d = df[df['outcome']==outcomes[i]]
                    b = [x + i*2 + v for x in [0,8,16]] 
                    
                    elast = [cd[exp][l][outcomes[i]][va].elasticities for l in limits]
                    me_upper = [cd[exp][l][outcomes[i]][va].ME_upper for l in limits]
                    me_lower = [cd[exp][l][outcomes[i]][va].ME_lower for l in limits]
                   
                    ax.plot(b, elast, marker=marker, linestyle='', 
                            color=plot_col, label=outcome_lab[i])
                    ax.plot(b, me_upper,marker='_',linestyle='',color=plot_col)
                    ax.plot(b, me_lower,marker='_',linestyle='',color=plot_col)
                    for j in b:
                        idx = list(b).index(j)
                        texts.append(plt.text(j, elast[idx]+2, str(round(elast[idx]))+"%", fontsize=13.5))
                        # ax.annotate(str(round(elast[idx]))+"%", 
                        #             (j, elast[idx]+2), fontsize=12)
                        ax.plot([j,j],[me_lower[idx], me_upper[idx]],
                                 color=plot_col, linewidth=lw)
            
            adjust_text(texts, only_move={'points':'y', 'text':'xy', 'objects':'xy'}, 
                        arrowprops=dict(arrowstyle="->", color='silver', lw=0.5))
            
            ax.legend(loc='upper right', bbox_to_anchor=(1.05, 1.25),
                      ncol=2, fancybox=True, shadow=True, 
                      title=cats[0]+" "+w+" "*(d-2*len(var))+cats[1]+" "+w)
            if s == 'r_':
                temp_str = ', restricted sample'
                font_size_title = 14
            else:
                temp_str = ''
                font_size_title = 14
            if exp=="b":
                ts = "Semi-elasticities and CIs for {}{} \n".format(name,
                                                                    temp_str)
            else:
                ts = "Semi-elasticities and CIs for setup {}, {}{} \n".format(
                     exp, name, temp_str)
            ax.set_title(ts, fontsize=font_size_title, x=0.35, y=0.95)
            # Make axis description without tickmarks. 
            ax.tick_params(axis='x', bottom=False, left=False) 
            plt.xticks([2.3, 10.4, 18.5], [100, 120, 130], size='large')
            ax.set_xlim(-1*v, max(b)+v)
            ax.set_ylim(-70, 30)
            #ax.get_xaxis().tick_bottom()
            ax.yaxis.set_major_formatter(mtick.PercentFormatter())
            ax.set_ylabel("Semi-elasticity", fontsize=13)
            ax.set_xlabel("Speed limit", fontsize=13)
            #fig.subplots_adjust(top=3)
            fig.savefig(os.path.join(pp['out_graphics'], #'new',
                                     s+'grf_CATES_'+var+'_'+exp+'.png'), 
                        bbox_inches="tight", dpi=400)
    return()


def make_geodf(pp, df, limits):

    pp_frame = os.path.join(pp['data_out_km'], 'frame.shp')
    pp_spatial_out = os.path.join(pp['fin_data'], 'analysis_data.shp')
    frame_df = pd.DataFrame.spatial.from_featureclass(pp_frame)
    # Join maximum extent data to the  frame. 
    frame_df.set_index('ID', drop=False, inplace=True)
    frame_df = frame_df.join(df)
    
    sl_cols = limits + ['maxspeed_none']
    # Create one variable for maxspeed for visualization.
    frame_df['maxspeed'] = pd.DataFrame(frame_df[sl_cols].idxmax(axis=1)).apply(
                                            lambda x: '_'.join(
                                                    str(x[0]).split('_')[1:]),
                                            axis=1)
    
    
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
    
    return()


def make_data_dicts(exps, limits, outcomes, df, name, cates_d, p_in, s):
    dd = dict(zip(exps,[dict(zip(limits,
                                 [dict(zip(outcomes,[0]*len(outcomes)))
                                  ]*len(limits)))]*len(exps)))
    for exp in exps:
        dd[exp] = {}
        for limit in limits:
            dd[exp][limit] = {}
            for outcome in outcomes:
                e, ey, lim, out = get_short_names(exp, limit, outcome)
                files_d = get_files_dict(p_in, s, exp+'_', lim, out)
                data = pd.read_csv(files_d[name], index_col=0)
                
                # For cates some additional data is needed.
                if name == "cates":
                    i_res = load_csv_and_fix_id(files_d['i_res']).join(df[
                                    ["Y_hat_"+ey+out, "W_hat_"+e+lim, limit, outcome]+\
                                    list(cates_d.keys())
                                ])
                    df = df.join(i_res[['tau_hat_oob', 'sigma_hat']].rename(columns={
                                                'tau_hat_oob':  'th_'+exp+'_'+lim+'_'+out,
                                                'sigma_hat':'sh_'+exp+'_'+lim+'_'+out}))
                    data = add_CI_bounds_and_elast(i_res, data, cates_d, e, 
                                                           ey, lim, out)
                    
                dd[exp][limit][outcome] = data.copy()

    return(dd, df)

def stars(estimate, std_err, N):
    
    if std_err != 0:
        t_stat = abs(estimate/std_err)
        if t_stat >= t.ppf(0.995, N):
            s = "$^{***}$"
        elif t_stat >= t.ppf(0.975, N):
            s = "$^{**}$"
        elif t_stat >= t.ppf(0.95, N):
            s = "$^{*}$"
        else:
            s = ""
    else:
        s=""
    
    return(s)

def get_cates_cutpoints_dict(df, cate_vars):
    
    cates_d = {}
    for var in cate_vars:
        
        if var == 'AADT':
            b = pd.DataFrame(df[var].sort_values().cumsum())
        else:
            b= pd.DataFrame(df[var].sort_values().cumsum()).join(df['AADT'])
            b['AADT'] = b['AADT'].cumsum()
        a = round(df['AADT'].sum()/2)
        
        lower_id = b[b['AADT']<a][var].idxmax()
        upper_id = b[b['AADT']>a][var].idxmin()
        
        cut_var = df[var].loc[[lower_id, upper_id]].mean()
        if cut_var == 0:
            prox_min = 0.01
            id_min = b[b[var]>=prox_min][var].idxmin()
            
            print(var,' -> cut_var=0, instead take', prox_min, 
                  'covers {}% of traffic'.format(
                          round((1-b.loc[id_min]["AADT"]/df['AADT'].sum())*100)))
            cates_d[var] = prox_min
        else:
            print(var, round(cut_var,5))
            cates_d[var] = cut_var
        
    return(cates_d)


def get_files_dict(p_in, s, exp_, lim, out):
    
    oob_files = os.path.join(p_in, "..", "oob_predictions")
    oob_files_me = os.path.join(p_in, "..", "oob_predictions", "max_ext")

    files_d = {
            "cates":   os.path.join(p_in, s+"cates_"+exp_+lim+'_'+out+".csv"),
            "omn_t":   os.path.join(p_in, s+"omn_test_"+exp_+lim+'_'+out+".csv"),
            "i_res_me":os.path.join(p_in, 
                                    s+"ind_res_max_ext_"+exp_+lim+'_'+out+".csv"),
            "i_res":   os.path.join(p_in, s+"ind_res_"+exp_+lim+'_'+out+".csv"),
            "var_imp": os.path.join(p_in, s+"var_imp_"+exp_+lim+'_'+out+".csv"),
            "oob_main_eff": os.path.join(oob_files, "oob_predictions_"+out+".csv"),
            "oob_prop_score":os.path.join(oob_files, "oob_predictions_"+lim+".csv"),
            "main_eff_me":os.path.join(oob_files_me, "max_ext_predictions_"+exp_+out+".csv"),
            "prop_score_me":os.path.join(oob_files_me, "max_ext_predictions_"+exp_+lim+".csv")
        }  
    if exp_ not in ['ab_']:
        files_d.update({
                "rate_q":    os.path.join(p_in, '..', 'rate', 
                                          s+'rate_qini_'+exp_+lim+'_'+out+".csv"),
                "rate_a":    os.path.join(p_in, '..', 'rate', 
                                          s+'rate_autoc_'+exp_+lim+'_'+out+".csv"),
                "rate_q_cr":    os.path.join(p_in, '..', 'rate', 
                                          s+'rate_qini_'+exp_+"cr_"+lim+'_'+out+".csv"),
                "rate_a_cr":    os.path.join(p_in, '..', 'rate', 
                                          s+'rate_autoc_'+exp_+"cr_"+lim+'_'+out+".csv")
            })
    
    return(files_d)  

def add_CI_bounds_and_elast(i_res, cates, cates_d, e, ey, lim, out):

    i_res["y_hat_0"] = i_res["Y_hat_"+ey+out]-i_res["W_hat_"+e+lim]*i_res["tau_hat_oob"]
    
    cates.loc['CI_upper'] = cates.loc["estimate"] + 1.96* cates.loc["std.err"]
    cates.loc['CI_lower'] = cates.loc["estimate"] - 1.96* cates.loc["std.err"]
    
    y_0 = i_res["y_hat_0"]
    ate_cols = ["ate", "atet", "ateut", "ate_ov"]
    
    elasticities = [elast(cates,x,"estimate",y_0,i_res,lim) for x in ate_cols]
    me_upper = [elast(cates,x,"CI_upper",y_0,i_res,lim) for x in ate_cols]
    me_lower = [elast(cates,x,"CI_lower",y_0,i_res,lim) for x in ate_cols]
    
    limit = 'maxspeed_'+lim
    n = [len(y_0), len(i_res[i_res[limit]==1]["y_hat_0"]),  
         len(i_res[i_res[limit]!=1]["y_hat_0"]), len(y_0)]
    
    for col in list(cates.columns)[4:]:
        k = "_".join(col.split("_")[2:])
        pos = col.split("_")[1]
        
        if pos == "h":
            i_res_temp = i_res[i_res[k]>cates_d[k]]
        elif pos == "l":
            i_res_temp = i_res[i_res[k]<=cates_d[k]]
            
        else:
            raise(ValueError("Something is wrong with the column names!"))
        y_0 =  i_res_temp["y_hat_0"]
        
        me_upper.append(elast(cates,col, "CI_upper", y_0))
        me_lower.append(elast(cates,col, "CI_lower", y_0))    
        elasticities.append(elast(cates,col, "estimate", y_0))
        n.append(len(y_0))
    cates.loc["elasticities"] = elasticities
    cates.loc["ME_upper"] = me_upper
    cates.loc["ME_lower"] = me_lower
    cates.loc["N"] = n
    
    return(cates)
    

def elast(cates, col, row, y_hat_0, i_res=None, lim=None):
    
    if col == "atet":
        limit = 'maxspeed_'+lim
        y_hat_0 = i_res[i_res[limit]==1]["y_hat_0"]
        
    if col == "ateut":
        limit = 'maxspeed_'+lim
        y_hat_0 = i_res[i_res[limit]!=1]["y_hat_0"]

    elast = cates[col][row]/np.mean(y_hat_0)*100
    
    
    return(elast)
    
    
def doubly_robust_ate(df, tau_col, y_col, y_hat_col, w_col, w_hat_col):
    
    prop_sc_weight = (df[w_col]-df[w_hat_col])/(df[w_hat_col]*(1-df[w_hat_col]))
    
    mu_xw = df[y_hat_col] + (df[w_col]-df[w_hat_col])*df[tau_col]
    
    correction = prop_sc_weight*(df[y_col] - mu_xw)
    
    
    dr_ate = np.mean(df[tau_col] + correction)
    
    return(prop_sc_weight, mu_xw, correction, dr_ate)


def get_short_names(exp, limit, outcome):
    
    if exp in ["b",'ab1', 'ab2', "ba"]:
        e = ""
    elif exp == "ab":
        e="a_"
    else:
        e=exp+"_"
    
    lim = limit.split("_")[1]
    
    if outcome == "total_rate":
        out_short = "t"
    else:
        out_short = "".join([x[:1] for x in outcome.split("_")])
    
    
    if exp == 'c':
        ey = ""
    else:
        ey = e
    
    return(e,ey,lim,out_short)






