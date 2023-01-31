# -*- coding: utf-8 -*-
"""
Add regional sociodemographic variables.

Merge to zeb data frame since they contain sk_kennung on frame ID.

Runtime: approx. 20 minutes. 

"""

import arcpy
import os
import pandas as pd
import numpy as np
import json
import random

from project_paths import get_project_paths, get_params
import src.process_geodata.functions.f14_sociodem_functions as sociodem

from arcgis import GIS

random.seed(42)

arcpy.env.overwriteOutput = True


pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))
pp_zeb = os.path.join(pp["data_out_df"], "zeb.csv")

# Load mapping from original variables names to new names. 
sd_lib = json.load(open(os.path.join(pp['lib'], 'sociodemography_lib.json'), 
                        encoding='UTF-8'))

sr = arcpy.SpatialReference(param['spacial_reference'])

pp_inkar_out = os.path.join(pp["data_out_inkar"],"kreise_inkar.shp")
inkar_files = [os.path.join(pp['data_in_inkar'],x) for x in 
                     os.listdir(pp['data_in_inkar']) if x.endswith('.csv')]
bkg_shps = [os.path.join(pp['data_in_BL'],x) for x in 
                     os.listdir(pp['data_in_BL']) if x.endswith('.shp')]
aggregates = ['Landkreis', 'Oberzentrum', 'krsfr. Stadt']


# Construct the "Kreise" dataframe from all inkar files.
df_krs = sociodem.get_krs_df(aggregates, inkar_files)
# Generate a random variable on the unit interval. 
df_krs['random'] = np.random.uniform(0, 1, len(df_krs))
sd_lib['rename_dict'].update({'random':'random'})
# Merge to shape file and adjust spatial reference system.
df = sociodem.merge_krs_to_bkg_and_clean(df_krs, bkg_shps[1], sd_lib)
df.spatial.to_featureclass(r"in_memory\kreise")
arcpy.management.Project(r"in_memory\kreise", pp_inkar_out, sr)
arcpy.management.Delete([r"in_memory\kreise"])

# Load in BISSTra frame.
pp_in_bisstra = os.path.join(pp["data_temp"], "BISSTra_lage.shp")
arcpy.management.Dissolve(pp_in_bisstra, r"in_memory\bisstra", "sk_kennung")

# Buffer BISStra
arcpy.analysis.Buffer(r"in_memory\bisstra", r"in_memory\bisstra_Buffer", 
                      "25 Kilometers", "FULL", "ROUND", "NONE", None, "PLANAR")

# Clip Buffer with union of kreise. 
arcpy.management.Dissolve(pp_inkar_out, r"in_memory\ger")
arcpy.analysis.Clip(r"in_memory\bisstra_Buffer", r"in_memory\ger", 
                    r"in_memory\bisstra_Buffer_Clip", None)

# Take average of var in buffer (weighted by area)
arcpy.management.CopyFeatures(pp_inkar_out, r"in_memory\inkar")
arcpy.analysis.Intersect(r"in_memory\bisstra_Buffer_Clip #;in_memory\inkar #", 
                         r"in_memory\intersect", "ALL", None, "INPUT")
arcpy.management.Dissolve(r"in_memory\intersect", r"in_memory\intersect_diss", 
                          "sk_kennung;FID_inkar", 
                          "emp_quo FIRST;pop_18_25 FIRST;pop_ol_65 FIRST;fem_share FIRST;hh_inc FIRST;pop_dens FIRST;pcar_dens FIRST;bip_p_cap FIRST;rurality FIRST;random FIRST", 
                          "MULTI_PART", "DISSOLVE_LINES")
# Add shape area in km^2. 
arcpy.management.AddGeometryAttributes(r"in_memory\intersect_diss", "AREA",
                                       Area_Unit="SQUARE_KILOMETERS")
# Make numerical data out of it 
df = pd.DataFrame.spatial.from_featureclass(r"in_memory\intersect_diss")
# Cleanup
arcpy.management.Delete([r"in_memory\inkar", r"in_memory\bisstra_Buffer_Clip",
                         r"in_memory\bisstra_Buffer",
                         r"in_memory\bisstra", r"in_memory\intersect",
                         r"in_memory\intersect_diss"])

# Sum Area by sk_kennung
df.drop(columns="SHAPE", inplace=True)
df.rename(columns=dict(zip([x for x in df.columns], 
                           [x.replace("FIRST_","") for x in df.columns])),
          inplace=True)

df['sk_area'] = df['POLY_AREA'].groupby(df['sk_kennung']).transform('sum')
df['krs_share'] = df['POLY_AREA']/df['sk_area']
# Derive shape of 
cols = ["emp_quo", "pop_18_25", "pop_ol_65", "fem_share", "hh_inc", 
        "pop_dens", "pcar_dens", "bip_p_cap", "rurality", "random"]

for col in cols:
    # Make auxiliary weighted cols
    df[col+'_w'] = df[col] * df['krs_share']
    # Sum it up
    df[col+'_out'] = df[col+'_w'].groupby(df['sk_kennung']).transform('sum')

df = df[['sk_kennung']+[x for x in df.columns if x.endswith('_out')]]
df.drop_duplicates(inplace=True)

df.rename(columns=dict(zip([x for x in df.columns if x.endswith('_out')], 
                            [x[:-4] for x in df.columns if x.endswith('_out')])),
    inplace=True)
# Merge it to numerical frame via the frame with BISStra Sektion numbers. 
df_zeb = pd.read_csv(pp_zeb)
df = df_zeb.merge(df, on="sk_kennung", how="left")

df.set_index('id', inplace=True)
df.index.name = None
df.to_csv(os.path.join(pp['data_out_df'],"zeb_inkar.csv"))

