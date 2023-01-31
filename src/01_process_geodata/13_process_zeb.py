# -*- coding: utf-8 -*-
"""
Process zeb data.
Runs roughly 1h.
"""

import arcpy
import os
import pandas as pd

from project_paths import get_project_paths, get_params
from src.process_geodata.functions.f10_ramps_functions import locate_points_and_to_df
from src.process_numdata.functions.f02_collapse_functions import long_to_wide

import src.process_geodata.functions.f13_zeb_functions as zeb

from arcgis import GIS

arcpy.env.overwriteOutput = True


## Enter and remove credentials. 
## Without this, it is not possible to use spatially enabled dataframes.
#g2 = GIS("https://www.arcgis.com", "username", "password")

pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))
sr = arcpy.SpatialReference(param['spacial_reference'])


zeb_files, pp_frame, pp_bisstra, pp_nk, pp_temp, pp_out_pts, pp_out_csv = (
                                                             zeb.get_paths(pp))





# Subset bisstra datasets as needed. 
arcpy.management.MakeFeatureLayer(pp_bisstra, r"in_memory\A_line", 
                              "Str_Klasse = 'A' And Sk_Subtyp_ = 'Abschnitt'")
# Load and subset network nodes. 
arcpy.management.MakeFeatureLayer(pp_nk, r"in_memory\NK_Points", 
                                  "NK_GeoStrK LIKE '%A%'")
# Somehow A_line gets lost through function take_lines_start_end_join_pts... 
# Thus copy fist. 
arcpy.management.CopyFeatures(r"in_memory\A_line", r"in_memory\A_line_")
# Get start and end vertice of bisstra network.
sdf_startpts, sdf_endpts = zeb.take_line_start_end_and_join_pts(    
    in_line=r"in_memory\A_line", join_pts= r"in_memory\NK_Points",
    search_radius="71 Meters")
# Identify if segment lies left or right, looking into stationary direction.
sdf_bisstra = zeb.find_side_and_join_to_bisstra(sdf_startpts=sdf_startpts, 
                                                sdf_endpts=sdf_endpts,
                                                bisstra_in=r"in_memory\A_line_")
# Save out to pp_temp.
sdf_bisstra.spatial.to_featureclass(pp_temp)





# Select bisstra separated by Left and Right side (in stationary direction).
arcpy.management.MakeFeatureLayer(pp_temp, r"in_memory\A_Right", "l_r = 'R'")
arcpy.management.MakeFeatureLayer(pp_temp, r"in_memory\A_Left", "l_r = 'L'")

# Required because of a weird bug, that didn't exist before and hopefully goes 
# away in the future:
arcpy.management.CopyFeatures(r"in_memory\A_Right", r"in_memory\A_R")
arcpy.management.CopyFeatures(r"in_memory\A_Left", r"in_memory\A_L")


# For 'Left' sections, reverse direction. Transforms input object.
arcpy.edit.FlipLine(r"in_memory\A_L")

# Generate Points along lines and add distance along segment to it. 
zeb.generate_and_locate_points(r"in_memory\A_R", r"in_memory\right_pts")
zeb.generate_and_locate_points(r"in_memory\A_L", r"in_memory\left_pts")

# Flip A_left again, to get directional mean right. 
arcpy.edit.FlipLine(r"in_memory\A_L")

# Split lines, get directional mean of each cut segement and merge to segment. 
zeb.split_line_at_points_and_add_dir_mean(r"in_memory\A_R", 
                                          r"in_memory\right_pts", 
                                          r"in_memory\right_split")
zeb.split_line_at_points_and_add_dir_mean(r"in_memory\A_L", 
                                          r"in_memory\left_pts", 
                                          r"in_memory\left_split")
# Join all info to startpoints of split lines. 
zeb.join_to_startpts(r"in_memory\right_split", r"in_memory\right_pts",
                     r"in_memory\right_pts_dir")
zeb.join_to_startpts(r"in_memory\left_split", r"in_memory\left_pts",
                     r"in_memory\left_pts_dir")

# Next, snap points to frame.  somehow... 

sdf_frame = zeb.load_frame_add_dir_by_direction(
        pp_frame, os.path.join(pp['data_out_df'], 'zahlst_dir_2018.csv'))

sdf_r = zeb.find_nearpoints_to_sdf_by_A(r"in_memory\right_pts_dir", 
                                       {"0": r"in_memory\frame_0", 
                                        "1": r"in_memory\frame_1"},
                                        "99 Meters", sdf_frame)
sdf_l = zeb.find_nearpoints_to_sdf_by_A(r"in_memory\left_pts_dir", 
                                        {"0": r"in_memory\frame_0", 
                                         "1": r"in_memory\frame_1"},
                                         "99 Meters", sdf_frame)
# Select the coordinates on the frame, that the points should be projected to.
# and subset columns to remove old coordinates and info not needed anymore. 
df_r = zeb.choose_coordinates_and_subset_df(sdf_r)
df_l = zeb.choose_coordinates_and_subset_df(sdf_l)
# Save out to temp. 
df_r.to_csv(os.path.join(pp['data_temp'], "zeb_frame_r.csv"))
df_l.to_csv(os.path.join(pp['data_temp'], "zeb_frame_l.csv"))



# For each row in sdf_r and sdf_l, choose point 0 or 1 on frame
# First step is to take closest point (NEAR_DIST_0 or NEAR_DIST_1)
# Check whether direction fits okay 
# Else check whether other direction fits well. 
# If not, assign no near point. 
df_r = pd.read_csv(os.path.join(pp['data_temp'], "zeb_frame_r.csv"))
df_l = pd.read_csv(os.path.join(pp['data_temp'], "zeb_frame_l.csv"))
# Load zeb data as df.
zeb_df = pd.concat([pd.read_excel(zeb_files[0]), pd.read_excel(zeb_files[1])],
                    ignore_index=True)



df_r['fmeas'] = df_r['fmeas'].round(0).apply(int)
df_l['fmeas'] = df_l['fmeas'].round(0).apply(int)
# Remove 0 points since they are artefact from overlapping BISSTRa shapes. 
df_l = df_l[df_l['fmeas']!=0]
# Same for endpoints of df_r that are not 0 or devisible by 100.
df_r = df_r[(df_r['fmeas']==0)|(df_r['fmeas'].apply(str).apply(lambda x: x[-2:])=='00')]


# Sort data from zeb_df too df_r and df_l, now. 
ident_cols = ['vnk', 'nnk', 'vst', 'bst', 'lage', 'fs']
value_cols = ['bauw_3', 'geb_15', 'sub_15', 'gw_15']
zeb_df = zeb_df[ident_cols + value_cols]

col_map_r = {'vnk':'vnk', 'nnk':'nnk', 'l_r':'lage', 'fmeas':'vst'}
col_map_l = {'vnk':'vnk', 'nnk':'nnk', 'l_r':'lage', 'fmeas':'bst'}
# Nochmal checken, ob die dtypes übereinstimmen oder irgendwas noch in int/str
# umgewandelt werden muss!

zeb_df_r, zeb_df_l, df_r, df_l = zeb.get_and_clean_dfs(zeb_df, df_r, df_l, 
                                              col_map_r, col_map_l, ident_cols)
# Merge value_cols of zeb to point data using the columns defined in col_map.
df_r = zeb.merge_zeb_to_df(zeb_df_r, df_r, col_map_r, value_cols)
df_l = zeb.merge_zeb_to_df(zeb_df_l, df_l, col_map_l, value_cols)

df_r = df_r[[x for x in df_r.columns if not x.endswith(('_x', '_y', ": 0"))]]
df_l = df_l[[x for x in df_l.columns if not x.endswith(('_x', '_y', ": 0"))]]



pp_temp = os.path.join(pp['data_temp'], "zeb_pts.csv")


df = df_r.append(df_l)
df.to_csv(pp_temp)


arcpy.management.XYTableToPoint(pp_temp, r"in_memory\proj_pts", 
                                'Near_X', 'Near_Y', None, sr)
arcpy.management.CopyFeatures(r"in_memory\proj_pts", pp_out_pts)





# Transform to numerical data. 
arcpy.management.CopyFeatures(pp_out_pts, r"in_memory\zeb_pts")
arcpy.management.AlterField(r"in_memory\zeb_pts", 'fmeas', 'sk_dist')

in_fields = [f.name for f in arcpy.ListFields(r"in_memory\zeb_pts") if 
             f.name.startswith(('geb', 'sub','bauw','gw', 'sk_kennung'))]
id_fields = ['vnk','nnk','l_r', 'sk_dist']
# Locate zeb points along the segments. 
df = locate_points_and_to_df(pp_frame, tolerance="0,1 Meters", 
                             pp_points=r"in_memory\zeb_pts", 
                             in_fields=in_fields+id_fields, 
                             clean_auf_ab=False, route_locations="ALL")
# Add startpoints (FMEAS=0) to frame segments.
df = zeb.add_startpoints(df, id_fields, set_to_start=10)

df_save = df.copy()

# Next, reshape to wide df.
# Set to unique index.   
df = df.reset_index().drop(columns='index')
df['id'] = df['ID']
df['num'] = df.sort_values(['FMEAS']).groupby(['id']).cumcount()
df.set_index(['num','id'], inplace=True)

df_wide = long_to_wide(df, stubnames=[x for x in df.columns if x not in 
                                      ['sk_dist', 'idx']])

unique_cols = ['vnk', 'nnk', 'sk_kennung', 'l_r']
df_wide.rename(columns={x: x[:-2] for x in [y+'_0' for y in unique_cols]}, 
                        inplace=True)
df_wide.drop(columns=[x for x in df_wide.columns if '_' in x and 
                      x.startswith(tuple(unique_cols)) and x[-1].isnumeric()],
             inplace=True)



dc = [x for x in df_wide.columns if x.startswith('FMEAS')]
 
# At the moment, maximum length is set to 120. Maybe change this, and consider 
# as best imputed value.  
df_wide[['length_'+str(i) for i in range(len(dc))]] = df_wide[dc].apply(
        zeb.get_length, args=[len(dc)], axis=1)

# Asphalt umwandeln in Asphalt share an gesamter bekannter Bauweise.
df_wide = df_wide
d = {}
for i in range(1,5):
    d.update({i:{'A': (df_wide[[x for x in df_wide.columns if 
                                      x.startswith('bauw_3_'+str(i))]]=='A')*1,
                       'B': (df_wide[[x for x in df_wide.columns if 
                                      x.startswith('bauw_3_'+str(i))]]=='B')*1}})

l = df_wide[[x for x in df_wide.columns if x.startswith('length')]].values
df_wide['A_l'] =  sum([d[i]['A'].mul(l).sum(axis=1) for i in range(1,5)])
df_wide['B_l'] =  sum([d[i]['B'].mul(l).sum(axis=1) for i in range(1,5)])
df_wide['asphalt_share'] = df_wide['A_l'] / (df_wide['A_l']+df_wide['B_l'])


# Aggregate values.
df, newvars = zeb.aggregate(df_wide, lanes=4,
                            newvar_dict={'bauw': ['mode', 'het'], 
                                         'sub': ['mean', 'max', 'var', 'num'], 
                                         'geb': ['mean', 'max', 'var', 'num']})


df = df[['sk_kennung']+newvars]
df = df.join(df_wide['asphalt_share'])

df.to_csv(pp_out_csv)









