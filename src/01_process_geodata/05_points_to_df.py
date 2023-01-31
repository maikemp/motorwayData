# -*- coding: utf-8 -*-
"""
Runtime: without making sideroads approx. 80 minutes
With making side roads: + ~ 50 minutes. 
"""
import ast
import os 
import arcpy 
import re
import pandas as pd

import src.process_geodata.functions.f05_pts_to_df_functions as pts_to_df
from project_paths import get_project_paths, get_params

arcpy.env.overwriteOutput = True


pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))
sr = arcpy.SpatialReference(param['spacial_reference'])


pp_frame = os.path.join(pp['data_out_km'], 'frame.shp')
# First, create a complete network of streets, exclude only motorways from it,
# and clip it with a generous buffer from the frame. 
buff_width_wide = "40 Meters"
buff_width_narr = "25 Meters"
pp_buffer_wide = r"in_memory\buffer_wide"
pp_buffer_narr = r"in_memory\buffer_narr"

arcpy.analysis.Buffer(pp_frame, pp_buffer_wide, buff_width_wide, 
                      "FULL", "ROUND", "NONE", None, "PLANAR")

arcpy.analysis.Buffer(pp_frame, pp_buffer_narr, buff_width_narr, 
                      "FULL", "FLAT", "NONE", None, "PLANAR")

dates = param['dates']

pp_side_roads = os.path.join(pp['data_temp'], 'side_roads.shp')
# Only make side roads again if something has changed in the OSM files.  
#pp_temp_BL = os.path.join(pp['data_temp'], 'BL_temp.shp')
#pts_to_df.make_side_roads(pp_buff=r"in_memory\buffer_wide", 
#                          pp_out=pp_side_roads, 
#                          pp_temp_BL=pp_temp_BL, 
#                          BL_path_list=[x[0] for x in 
#                                         os.walk(pp['data_in_gf'])][1:],
#                          filesname='gis_osm_roads_free_1.shp',
#                          sr=sr)

# Process accidents and put labels on them.

with arcpy.da.SearchCursor(pp_frame, ['ID']) as cursor:
        id_list = sorted({row[0] for row in cursor})

for date in dates: 
    print(date)

    # Get names of the paths etc. 
    yearmatch = re.match('\d{4}', date)
    year = yearmatch[0]
    pp_unf_in = os.path.join(pp['data_in_unfall'], 
                             'Unfallorte_'+year+'_Shapefile', 
                             'Unfallorte'+year+'_LinRef.shp')
    pp_unf_temp = os.path.join(pp['data_temp'], 'proje_unf'+year+'.shp')
    pp_unf_temp_gdb = os.path.join(pp['data_temp_gdb'], 'proje_unf'+year)
    pp_unf_out = os.path.join(pp['data_out_unf'], 'A_unfaelle_'+year+'.shp')
    pp_unf_df_out = os.path.join(pp['data_out_df'],'unf_df_'+date+'.csv')
    
    # Change the projection of the accident data
    arcpy.management.Project(pp_unf_in, pp_unf_temp, sr)
    
    print('Start with accidents in', year)
    pts_to_df.Select_and_label_Points(
            path_in_pts = pp_unf_temp,
            path_temp_pts = pp_unf_temp_gdb,
            path_out_pts = pp_unf_out,
            path_frame = pp_frame,
            path_buffer_wide = pp_buffer_wide,
            path_buffer_narr = pp_buffer_narr,
            path_side_roads = pp_side_roads,
            buff_width_narr = buff_width_narr,
            sr = sr,
            info_freq = 1000)

    # Next, make dictionary out of it.    
    order_var = 'from_start'
    excl_fields_list = ['FID','Join_Count','TARGET_FID', 'JOIN_FID', 'LINREFX', 
                        'LINREFY', 'NEAR_FID', 'NEAR_DIST', 'NEAR_X', 'NEAR_Y', 
                        'NEAR_FC', 'in_seg_ID', 'ref', 'seg_ID', 'par_seg_ID', 
                        'one_direct', 'Shape', 'UJAHR']
    
    # Make the wanted dictionary. 
    unf_dict_sort, varlist = pts_to_df.Points_to_frame_ID_dict(pp_unf_out, 
                                                excl_fields_list, order_var)
    
    unf_df = pts_to_df.point_dict_to_df(unf_dict_sort, varlist[1:], 
                                        'from_start',id_list)
    
    unf_df.to_csv(pp_unf_df_out)


    


# Make BL and RBZ data frame
# May have to delete outputs in pp_BL_out before, due to bug in arcpy (?).
dfs = {}
for level in ['LAN', 'KRS']:
    
    pp_BL_in = os.path.join(pp['data_in_BL'], 'VG250_'+level+'.shp')
    pp_BL_out = os.path.join(pp['data_out_BL'], level+'_polygon')
    
    arcpy.management.Project(pp_BL_in, pp_BL_out, sr)
    
    BL_dict = pts_to_df.Polygon_to_frame_ID_dict(pp_BL_out+'.shp', 
                                                 pp_frame, ['GEN'])
    df = pd.DataFrame.from_dict(BL_dict, orient='index',columns=['GEN_'+level])
    dfs.update({level:df})

full_data = dfs['LAN'].join(dfs['KRS'], how='outer')

full_data.to_csv(os.path.join(pp['data_out_df'],'BL.csv'))






    


