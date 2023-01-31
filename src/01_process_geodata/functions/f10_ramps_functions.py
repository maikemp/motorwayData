# -*- coding: utf-8 -*-
"""

"""


import os
import arcpy
from geojson import dump
import json
import pandas as pd
import numpy as np 

from src.process_geodata.functions.f01_process_ot_functions import (
                                                            cr_uniform_geojson,
                                                            cr_prop_list)

arcpy.env.overwriteOutput = True


def reshape_to_wide(df):
    
    n_rep_cols = max(df['ID'].value_counts())
    cols = [x+'_'+str(i) for i in range(n_rep_cols) for x in df.columns if x !='ID']
    
    data = {}
    
    for i in range(len(df)):
        row = df.iloc[i]
        if row.ID in list(data.keys()):
            data[row.ID].update({
                row.dist: list(row[[x for x in df.columns if x != 'ID']])})
        else:
            data.update({row.ID: {row.dist: list(row[[
                                       x for x in df.columns if x != 'ID']])}})
    
    for key in data.keys():
        a = dict(sorted(data[key].items()))
        b = []
        for k in a.keys():
            b.extend(a[k])
        data.update({key: b + list(np.repeat(np.nan, len(cols)-len(b)))})
    df = pd.DataFrame.from_dict(data, orient='index', columns=cols)
    
    return(df)
    

def locate_points_and_to_df(pp_lines, tolerance, pp_points, in_fields, 
                            id_field="ID", clean_auf_ab=True, 
                            route_locations="FIRST"):
    
    arcpy.CopyFeatures_management(pp_lines, r"in_memory\frame")
    
    arcpy.AddField_management(r"in_memory\frame", "F", "Double")
    arcpy.AddField_management(r"in_memory\frame", "T", "Double")
    
    c = arcpy.da.UpdateCursor(r"in_memory\frame", ['SHAPE@',"F", "T"])
    for r in c:
        r[1] = 0
        r[2] = r[0].length
        c.updateRow(r)
    
    arcpy.lr.CreateRoutes(r"in_memory\frame", id_field, r"in_memory\routes", 
                          "TWO_FIELDS", "F", "T", "UPPER_LEFT", 1, 0, "IGNORE",
                          "INDEX")
    
    arcpy.LocateFeaturesAlongRoutes_lr(pp_points, 
                                       r"in_memory\routes", 
                                       id_field, 
                                       tolerance, 
                                       r"in_memory\location_table", 
                                       id_field+" POINT FMEAS",
                                       route_locations=route_locations)

    # Select columns.
    columns = [x for x in [id_field]+[id_field+'2']+['FMEAS']+in_fields
               if x in [f.name for f in arcpy.ListFields(r"in_memory\location_table")]]
    # Extract data.
    data = [r for r in arcpy.da.SearchCursor(r"in_memory\location_table", 
                                             columns)]
    df = pd.DataFrame(data, columns=columns)

    if clean_auf_ab==True:
        # Drop duplicates and falsely matched links.
        df.drop_duplicates(inplace=True)
        df = df[df['ID']==df['ID2']]
        
        # Clean 
        df['Type'] = 'Entry'
        df.loc[df['MERGE_SRC']==r"in_memory\abfahrt",'Type']='Exit'
        
        df.drop(['ID2', 'MERGE_SRC'], axis=1, inplace=True)
        df.rename({'FMEAS':'dist'}, axis=1, inplace=True)
    #    
    #    df.set_index('ID', inplace=True)
    #    df.sort_index(inplace=True)
    
    return(df)


def get_links_from_ot_data(pp_in, pp_temp, pp_out, in_fields, sr):
    
    # Process data from overpass query.
    geojson = json.load(open(pp_in, encoding="utf8"))
    in_list, out_list, prop_count = cr_prop_list(geojson,  50000, in_fields)
    
    geojson = cr_uniform_geojson(geojson, in_list, out_list, miss_val='missing')
    
    with open(pp_temp, 'w') as f:
        dump(geojson, f)
        
    arcpy.JSONToFeatures_conversion(pp_temp, r"in_memory\links", 'Polyline')
    
    if arcpy.Exists(pp_out+'.shp'):
        files = [x for x in os.walk(os.path.join(pp_out,'..'))][0][2]
        for file in files:
            if file.startswith(pp_out.split('\\')[-1]):
                os.remove(os.path.join(pp_out,'..', file))
        
    arcpy.management.Project(r"in_memory\links", pp_out, sr)


def find_linkpoints_start_and_end(pp_in, pp_control, control_query, invert, 
                                  pp_out, pp_frame, intersect_tolerance):
    
    # Dissolve them to aviod confounding start and endpoints
    arcpy.management.Dissolve(pp_in, r"in_memory\links", None, None, 
                              "MULTI_PART", "UNSPLIT_LINES")

    # Load the relevant elements from the BASt network to select relevant links
    #    from the OSM links. 
    arcpy.Select_analysis(pp_control, r"in_memory\bast_links", control_query)
    arcpy.Buffer_analysis(r"in_memory\bast_links", r"in_memory\bast_buffer", 
                          "20 Meters", "FULL", "ROUND", "NONE", None, "PLANAR")
    # Buffer the control layer and select either only the links inside of it 
    # or their inverse. 
    links = arcpy.management.SelectLayerByLocation(r"in_memory\links", 
                "INTERSECT", r"in_memory\bast_buffer", "1 Meters", 
                "NEW_SELECTION", invert)
    arcpy.CopyFeatures_management(links, r"in_memory\sel_links")

    # Select startpoints of link segments and choose those that lie within
    # intersect_tolerance of the frame to select abfahrten.
    arcpy.management.FeatureVerticesToPoints(r"in_memory\sel_links", 
                                             r"in_memory\startpts", "START")
    arcpy.Intersect_analysis([r"in_memory\startpts", pp_frame], 
                             r"in_memory\abfahrt", "ALL", intersect_tolerance, 
                             "POINT")
    # Do the same for endpoints and make them the auffahrten.
    arcpy.management.FeatureVerticesToPoints(r"in_memory\sel_links", 
                                             r"in_memory\endpts", "END")
    arcpy.Intersect_analysis([r"in_memory\endpts", pp_frame], 
                             r"in_memory\auffahrt", "ALL", intersect_tolerance, 
                             "POINT")
    # Merge auf- und abfahrten and save out to pp_out
    arcpy.Merge_management([r"in_memory\auffahrt", r"in_memory\abfahrt"], 
                           r"in_memory\links", add_source='ADD_SOURCE_INFO')
    
    if arcpy.Exists(pp_out):
        files = [x for x in os.walk(os.path.join(pp_out,'..'))][0][2]
        for file in files:
            if file.startswith(pp_out[:-4].split('\\')[-1]):
                os.remove(os.path.join(pp_out,'..', file))
    # First argument in spatial join can't be in in_memory. 
    arcpy.CopyFeatures_management(r"in_memory\links", pp_out)
    # Spatial join with input to get attributes of the link.
    arcpy.SpatialJoin_analysis(pp_out, pp_in, 
                               r"in_memory\links_info", "JOIN_ONE_TO_ONE", 
                               "KEEP_ALL", match_option="CLOSEST")
    
    arcpy.CopyFeatures_management(r"in_memory\links_info", pp_out)


def merge_road_types(path_list, sr, files_name, type_query, pp_temp, pp_out):
    i=0
    for BL in path_list:
        path = os.path.join(BL, files_name)
        arcpy.Select_analysis(path, r"in_memory\BL", type_query)
        arcpy.management.Project(r"in_memory\BL", pp_temp, sr)
        
        if i == 0:
            arcpy.CopyFeatures_management(pp_temp, pp_out)
        else:
            arcpy.Merge_management([pp_out, pp_temp], 
                                   r"in_memory\links")
            arcpy.CopyFeatures_management(r"in_memory\links", pp_out)
        i+=1
        print(i)

