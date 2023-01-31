"""
"""

import arcpy
import pandas as pd





def merge_krs_to_bkg_and_clean(df_krs, bkg_path, sd_lib):
    
    # Take bkg_shps[1] and dissolve by Kreis ID.
    arcpy.management.Dissolve(bkg_path, r"in_memory\kreise", "AGS")
    sdf_krs = pd.DataFrame.spatial.from_featureclass(r"in_memory\kreise")
    sdf_krs["AGS"] = sdf_krs["AGS"].apply(int)
    
    df = sdf_krs.merge(df_krs, left_on=["AGS"], right_on=['Kennziffer'])
    df.rename(columns=sd_lib["rename_dict"], inplace=True)
            
    # Drop due to missing entries?/variance. 
    for var in df.columns:
        try:
            a = round(df[var].var(),2)
            if a==0:
                df.drop(columns=[var], inplace=True)
                print("Dropped Variable: {} due to missing variance.".format(
                        var))
        except:
            pass
    
    return(df)


def get_krs_df(aggregates, inkar_files):
    '''
    '''    
        
    df_dict = dict(zip(aggregates, [pd.DataFrame()]*3))
    for f in inkar_files: 
        df = pd.read_csv(f, sep=";", thousands=".", decimal=",")
        df = df.iloc[1:]
        for aggr in aggregates:
            if sum(df["Aggregat"]==aggr)>0:
                if len(df_dict[aggr])==0:
                    df_dict[aggr] = df[df["Aggregat"]==aggr]
                else:
                    df_dict[aggr] = df_dict[aggr].merge(
                            df[df["Aggregat"]==aggr], 
                            on=['Kennziffer', 'Raumeinheit', 'Aggregat'])
    
    df_krs = df_dict['krsfr. Stadt'].append(df_dict['Landkreis'])
    df_krs['Kennziffer'] = df_krs['Kennziffer'].apply(int)
    
    return(df_krs)

