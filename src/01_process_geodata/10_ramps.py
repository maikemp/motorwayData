# -*- coding: utf-8 -*-
"""
Runtime: approx: 2 minutes. 
"""

import arcpy
import os
import pandas as pd

import src.process_geodata.functions.f10_ramps_functions as ramps
from project_paths import get_project_paths, get_params

arcpy.env.overwriteOutput = True


pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))

sr = arcpy.SpatialReference(param['spacial_reference'])

# 1. Select and meger all junctions
pp_temp = os.path.join(pp['data_temp'], 'links_temp.shp')
pp_full_temp = os.path.join(pp['data_temp'], 'links.shp')
#BL_path_list = [x[0] for x in os.walk(pp['data_in_gf'])][1:]
pp_frame = os.path.join(pp['data_out_km'], 'frame.shp')
pp_control = os.path.join(pp['data_out_bast'], 'BFStr_Netz_SKPolyline.shp')

pp_out = os.path.join(pp['data_out_links'], 'main_links.shp')
pp_out_sec = os.path.join(pp['data_out_links'], 'secondary_links.shp')
pp_out_df_main = os.path.join(pp['data_out_df'], 'main_links.csv')
pp_out_df_sec = os.path.join(pp['data_out_df'], 'sec_links.csv')

pp_geojson = os.path.join(pp['data_in_ot'], 
                          'motorway_links_'+param['reference_date']+'.geojson')
pp_geojson_out = os.path.join(pp['data_temp'], 
                          'motorway_links_'+param['reference_date']+'.geojson')
pp_ot_out = os.path.join(pp['data_out_otshp'], 'links'+param['reference_date'])

in_fields = ['lanes']

ramps.get_links_from_ot_data(pp_in=pp_geojson,
                             pp_temp=pp_geojson_out, 
                             pp_out=pp_ot_out,
                             in_fields=in_fields, 
                             sr=sr)

# Replaced by ot data: see above.
#merge_road_types(BL_path_list, sr, 'gis_osm_roads_free_1.shp',
#                 "\"fclass\"='motorway_link'", pp_temp, pp_full_temp)


# Find primary link points.
ramps.find_linkpoints_start_and_end(
                       pp_in=pp_ot_out+'.shp', pp_control=pp_control,
                       control_query="Str_Klasse = 'A' And Sk_Subtyp_ = 'Ast'",
                       invert='NOT_INVERT', pp_out=pp_out, 
                       pp_frame=pp_frame,
                       intersect_tolerance="0.4 Meters")
# Find secondary link points. (Not in Bast Netz, thus only for Raststätten etc)
ramps.find_linkpoints_start_and_end(
                       pp_in=pp_ot_out+'.shp', pp_control=pp_control,
                       control_query="Str_Klasse = 'A' And Sk_Subtyp_ = 'Ast'",
                       invert='INVERT', pp_out=pp_out_sec, 
                       pp_frame=pp_frame,
                       intersect_tolerance="0.4 Meters")

# Next, get position of link point along the frame segment. 
df_main = ramps.locate_points_and_to_df(pp_frame, tolerance="0,4 Meters", 
                                        pp_points=pp_out, 
                                        in_fields=in_fields+['MERGE_SRC'])
df_sec = ramps.locate_points_and_to_df(pp_frame, tolerance="0,4 Meters", 
                                       pp_points=pp_out_sec, 
                                       in_fields=in_fields+['MERGE_SRC'])
### had to do a manual work-around with arcgis pro due to an apparent bug 
### (or mistake?) in the above funtion. thus here: a manual import of the
### intermediary results. 
#df_main = pd.read_csv('...\\main_link_temp_df.csv')
#df_sec = pd.read_csv('...\\sec_link_temp_df.csv')


# Delete links at the beginning of the whole road segment
#df_main = df_main.drop(df_main[(df_main['ID'].str.endswith('0000')) 

#                             & (round(df_main['dist'])==0)].index)
# Aggregate to IDs: link1 (type, level, lanes, dist), link2 (1,2,3,4) ...
df_main = ramps.reshape_to_wide(df_main)
df_sec = ramps.reshape_to_wide(df_sec)

df_main.to_csv(pp_out_df_main)
df_sec.to_csv(pp_out_df_sec)








# Write them to a dataframe with ids from the frame. 
#data = {}
#columns = ['main_auffahrt', 'main_abfahrt', 'sec_auffahrt', 'sec_abfahrt']
#with arcpy.da.SearchCursor(pp_frame, ['ID']) as c:
#    for r in c:
#        data.update({r[0]:[0, 0, 0, 0]})
#
#with arcpy.da.SearchCursor(pp_out, ['ID', 'MERGE_SRC']) as c:
#    for r in c:
#        if r[1] == r'in_memory\auffahrt':
#            data[r[0]] = list(np.array(data[r[0]]) + np.array([1,0,0,0]))
#        if r[1] == r'in_memory\abfahrt':
#            data[r[0]] = list(np.array(data[r[0]]) + np.array([0,1,0,0]))
#            
#with arcpy.da.SearchCursor(pp_out_sec, ['ID', 'MERGE_SRC']) as c:
#    for r in c:
#        if r[1] == r'in_memory\auffahrt':
#            data[r[0]] = list(np.array(data[r[0]]) + np.array([0,0,1,0]))
#        if r[1] == r'in_memory\abfahrt':
#            data[r[0]] = list(np.array(data[r[0]]) + np.array([0,0,0,1]))
#
## Set all values >1 to 1.
#for key in data.keys():
#    data[key] = [x if x<=1 else 1 for x in data[key]]
## To df and save out.
#df = pd.DataFrame.from_dict(data, orient='index', columns=columns)
#df.sort_index(inplace=True)
#df.to_csv(pp_out_df)





