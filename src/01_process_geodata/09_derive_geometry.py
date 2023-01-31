# -*- coding: utf-8 -*-
"""
Derive addiotional geometrical and locational properties of the 
road. 

Curvature, Number  of curves, length of straight section, 
length of entry/exit lane, connector/junction area, 
'Einflussbereich', number of digits of motorway.

Runtime: approx. 45 Minutes
"""

import arcpy
import os
from arcgis import GIS

import src.process_geodata.functions.f09_geometry_functions as geometry

from project_paths import get_project_paths, get_params
from src.process_geodata.functions.f10_ramps_functions import locate_points_and_to_df

arcpy.env.overwriteOutput = True


pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))

# Input, temporary and output paths. 
pp_frame = os.path.join(pp['data_out_km'], 'frame.shp')
pp_frame_split = r"in_memory\frame_split"
pp_join_pts = r"in_memory\join_pts"
pp_out_csv = os.path.join(pp['data_out_df'], 'geometry.csv')


# Dissolve frame, cut into 100m pieces and match linear directional mean to it. 
geometry.segment_frame_and_get_ldm(pp_frame, pp_frame_split, 
                                   split_length="100 Meters")
# Get vertices, join them and derived difference between them.
geometry.get_and_join_vertices(pp_frame_split, pp_join_pts)

# Derive where the point lies along the routes (by ref and one_direct)
# Locate zeb points along the segments. 
df = locate_points_and_to_df(pp_frame, tolerance="0,1 Meters", 
                             pp_points=pp_join_pts, 
                             in_fields=['ref','one_direct','angle_diff'], 
                             clean_auf_ab=False, route_locations="ALL")

# Sort values, select columns and reshape. 
df_wide, cols, sharp_turn = geometry.sort_select_and_reshape(df)

# Derive relevant variables. 
df_wide['right_turn'] = df_wide[cols].apply(lambda x: sum([abs(i) for i in list(x.values) if i >0]), axis=1)
df_wide['left_turn'] = df_wide[cols].apply(lambda x: sum([abs(i) for i in list(x.values) if i <0]), axis=1)
df_wide['total_turn'] = df_wide['right_turn'] + df_wide['left_turn']
df_wide['sharp_turns'] = df_wide[cols].apply(lambda x: len([i for i in list(x.values) if abs(i) > sharp_turn]), axis=1)

df_wide[['right_turn', 'left_turn', 'total_turn', 'sharp_turns']].to_csv(pp_out_csv)

