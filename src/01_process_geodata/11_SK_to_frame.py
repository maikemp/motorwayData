# -*- coding: utf-8 -*-
"""
Assign a unique value for each BASt Sector of the motorway and derive node 
area as well as a location id. 

Runtime: Approx. 20 Minutes.

"""

import arcpy
import os
import pandas as pd

from project_paths import get_project_paths, get_params
from arcgis import GIS
import src.process_geodata.functions.f11_sk_functions as sk

arcpy.env.overwriteOutput = True


pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))

sr = arcpy.SpatialReference(param['spacial_reference'])


pp_frame = os.path.join(pp['data_out_km'], 'frame.shp')
pp_sk = os.path.join(pp['data_out_bast'], 'BFStr_Netz_SKPolyline.shp')
pp_nk = os.path.join(pp['data_out_bast'], 'BFStr_Netz_NKPoint.shp')
pp_links = os.path.join(pp['data_out_links'], 'main_links.shp')

pp_out_sk = os.path.join(pp['data_out_df'], 'sk.csv')
pp_out = os.path.join(pp['data_out_df'], 'sk_node_loc.csv')

arcpy.management.CopyFeatures(pp_sk, r"in_memory\bast_sk")
arcpy.management.CopyFeatures(pp_frame, r"in_memory\frame")


arcpy.analysis.Select(r"in_memory\bast_sk", r"in_memory\bast_sk_a", 
                      "Str_Klasse = 'A' And Sk_Subtyp_ = 'Abschnitt'")
arcpy.analysis.Buffer(r"in_memory\bast_sk_a", r"in_memory\sk_buffer", 
                      "50 Meters", "FULL", "FLAT", "LIST", "Sk_Kennung", 
                      "PLANAR")
arcpy.analysis.SpatialJoin(r"in_memory\frame", r"in_memory\sk_buffer", 
                           r"in_memory\frame_sk", "JOIN_ONE_TO_ONE", 
                           "KEEP_ALL", match_option="CLOSEST")

cursor = arcpy.da.SearchCursor(r"in_memory\frame_sk", ['ID', 'Sk_Kennung'])
df_dict = {}
SK_id_dict = {}
i=0

for r in cursor:
    df_dict.update({r[0]: r[1]})
    if r[1] not in list(SK_id_dict.keys()):
        SK_id_dict.update({r[1]:i})
        i += 1
 
# Replace the complicated value by a new id.       
for i in df_dict.keys():
    df_dict.update({i: [df_dict[i], SK_id_dict[df_dict[i]]]})

df = pd.DataFrame.from_dict(df_dict, orient='index', columns=['sk_kennung', 
                                                              'Sector'])
df.sort_index(inplace=True)
df.to_csv(pp_out_sk)


# Continue with "Einflussbereiche" arcpy.management.CopyFeatures(pp_frame, r"in_memory\frame")
arcpy.management.CopyFeatures(pp_links, r"in_memory\links")
df['node_area']= 0

for a in list(set(r[0] for r in arcpy.da.SearchCursor(r"in_memory\frame", "ref"))):
    for d in [0,1]: 
        # Select links and frame parts. 
        arcpy.analysis.Select(r"in_memory\links", r"in_memory\links_select", 
                              "one_direct = '{}' And ref = '{}'".format(d, a))
        arcpy.analysis.Select(r"in_memory\frame", r"in_memory\frame_select", 
                              "one_direct = '{}' And ref = '{}'".format(d, a))     
        # Devide by entry/exit. 
        arcpy.analysis.Select(r"in_memory\links_select", r"in_memory\ls_Entry", 
                              """MERGE_SRC = 'in_memory\\auffahrt'""")        
        arcpy.analysis.Select(r"in_memory\links_select", r"in_memory\ls_Exit", 
                              """MERGE_SRC = 'in_memory\\abfahrt'""")   
        # Einflussbereich: 500 m vor Knotenpunkt.
        arcpy.analysis.Buffer(r"in_memory\ls_Exit", r"in_memory\Exit_buff",
                              "750 Meters")
        # Einflussbereich: 300 m nach Knotenpunkt.
        arcpy.analysis.Buffer(r"in_memory\ls_Entry", r"in_memory\Entry_buff",
                              "550 Meters")
        # Merge back together and dissolve.
        arcpy.management.Merge([r"in_memory\Entry_buff",r"in_memory\Exit_buff"],
                               r"in_memory\nodes_merged")
        arcpy.management.Dissolve(r"in_memory\nodes_merged",r"in_memory\nodes",
                                  None, None, "SINGLE_PART", "DISSOLVE_LINES")
        
        arcpy.analysis.Clip(r"in_memory\frame_select", r"in_memory\nodes", 
                            r"in_memory\frame_Clip", None)
        arcpy.management.AddGeometryAttributes(r"in_memory\frame_Clip", 
                                               "LENGTH", "METERS")
        
        df_temp = pd.DataFrame.spatial.from_featureclass(r"in_memory\frame_Clip")
        df_temp = df_temp.set_index('ID')[['LENGTH']]
        df_temp.index.name = None
        df_temp['node_area'] = df_temp['LENGTH']/500
        
        df.loc[df_temp.index, 'node_area'] = df_temp['node_area']

# Add location_ids.Define locations by borders while dissolving city states and
# keeping whole motorways of less than 30km.
df_bi = sk.get_locations(
        pp_sk, pp_nk, pp_frame, 
        pp_temp=os.path.join(pp['data_temp'], 'location.shp'),
        pp_temp_NK=os.path.join(pp['data_temp'], 'NK_.shp'), how='border', 
        dissolve_city_states=True, dissolve_short_A=30000
    )

# Join to frame df via SK_Kennung. 
df = df.reset_index().merge(df_bi, left_on='sk_kennung', right_on='Sk_Kennung', 
        how='left').set_index('index')
df.drop(columns='Sk_Kennung', inplace=True)
df.to_csv(pp_out)

