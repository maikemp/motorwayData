# -*- coding: utf-8 -*-
"""
Warning: due to arcpy bugs, only works after a fresh kernel restart and when
 none of output files exists before.  


Runtime: approx. 30 minutes. 
"""

import arcpy

from geojson import dump
import json
import os

import src.process_geodata.functions.f01_process_ot_functions as process_ot
from project_paths import get_project_paths, get_params

arcpy.env.overwriteOutput = True


pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))
sr = arcpy.SpatialReference(param['spacial_reference'])
dates = param['dates']


proc_ot_param = json.load(open(os.path.join(pp['lib'], 'proc_ot_lib.json')))
add_list = proc_ot_param["add_list"]
add_sw = proc_ot_param["add_sw"]
excl_list = proc_ot_param["excl_list"]
excl_sw = proc_ot_param["excl_sw"]
replace_dict = proc_ot_param['repl_if_miss']
min_num = proc_ot_param['min_keep_feature_num']


# Load reference date geojson as feature name reference.
geojson = json.load(open(os.path.join(pp['data_in_ot'],
    'motorways_DE_'+param['reference_date']+'.geojson'), encoding="utf8"))

# Create lists of all wanted and unwanted properties.
in_list, out_list, prop_count = process_ot.cr_prop_list(
        geojson, min_num, add_list, add_sw, excl_list, excl_sw
    )

##############################################################################
# For Testing:
in_list_line = ['@id', 'ref', 'official_ref', 'maxspeed', 'bridge', 'tunnel', 
                'lit',  'under_construction', 'width', 
                'shoulder', 'shoulder:right', 'shoulder:left', 
                'hgv', 'maxspeed:hgv', 'overtaking:hgv', 'overtaking:trailer',
                'maxspeed:variable', 'maxspeed:wet', 'maxspeed:variable:default',
                'maxspeed:conditional', 'maxspeed:reason', 
                'maxspeed:hgv:conditional', 'overtaking:conditional', 
                'overtaking:hgv:conditional', 'overtaking:trailer:conditional', 
                'temporary:maxspeed', 
                'fixme', 'lanes', 'construction', 'note:maxspeed', 
                'maxspeed:note', 'note', 'note:de']
in_list_point = ['highway', 'ref']

in_list_new = list(set(in_list_line).union(set(in_list_point)))

out_list_new = set(in_list).union(set(out_list)) - set(in_list_new)
in_list = in_list_new
out_list = out_list_new
# These are the only fields that will be created for Points, and only for Points. 
# Execute geojson manipulation.
# May also write a function for this.
for date in dates:
    
    print('Start harmonizing geojson for date:',date)
    # Set paths and load geojson data from overpass turbo. 
    pp_in = os.path.join(pp['data_in_ot'], 'motorways_DE_'+date+'.geojson')
    pp_out = os.path.join(pp['data_temp'],'motorways_DE_'+date+'.geojson' )
    geojson = json.load(open(pp_in,encoding="utf8"))
    
    # Make the geojson uniform by setting miss_val for all in_list and removing
    # all out_list.
    geojson = process_ot.cr_uniform_geojson(geojson, in_list, 
                                            out_list, miss_val='missing')
    
    # Defines all fields that have to correspond to the respective expression.
    del_dict={'ref':r'A{1}\s?\d{1,3}'}   
    # Keep only featuers for which regular expression in respective field is 
    # matched, or for which the value is missing.
    geojson = process_ot.del_features(geojson, del_dict, miss_val='missing')
    
    # Rename some fields so they are still recognizable after turning them into
    # a shapefile, where field names have only up to 10 characters. 
    geojson = process_ot.rename_fields(geojson, proc_ot_param['rename_dict'])
    
    # Define fields as keys which will be changed if they meet the second list
    # entry instead of the first, according to the third list entry expression
    # First entry: enter space if missing, second entry: delete everything after
    # a semicolon so that A 1; A 61 is now counted for A 1.
    change_dict = {
            'ref': [[r'A{1}\s{1}\d{1,3}', r'A{1}\d{1,3}',
                    [r'(A{1})(\d{1,3})',r'\1 \2']],
                   [r'A{1}\s?\d{1,3}$', r'A{1}\s?\d{1,3};',
                    [r';.*$', ""]]]}   
    geojson = process_ot.change_features(geojson, change_dict)

    geojson = process_ot.replace_if_missing(geojson, replace_dict, 
                                            miss_val="missing")
    
    with open(pp_out, 'w') as f:
        dump(geojson, f)


# Covnert json to shapefile for each geometry type. 
for date in dates:
    for geomt in ['Point', 'Polyline']:
        print('Start with:', geomt, 'in', date)
        # Define required paths.
        pp_in = os.path.join(pp['data_temp'], 'motorways_DE_'+date+'.geojson')
        pp_out = os.path.join(pp['data_out_otshp'], 'motorways_DE_'+date+geomt)
        
        # Produce shapefiles and store them in memory.
        arcpy.JSONToFeatures_conversion(pp_in, r"in_memory\features", geomt)

        # Assign projection and save out to pp_out.       
        arcpy.management.Project(r"in_memory\features", pp_out, sr)
        print(geomt, date)
        if geomt == 'Polyline':
            remove_fields = [x[:10].replace(':','_') for x in in_list if x not in in_list_line]
            arcpy.DeleteField_management(pp_out+'.shp', remove_fields)
        else:
            remove_fields = [x[:10].replace(':','_') for x in in_list if x not in in_list_point]
            arcpy.DeleteField_management(pp_out+'.shp', remove_fields)
        arcpy.Delete_management(r"in_memory\features")

# Only for the reference file: impute missings from field "ref" to get uniform 
# container file and fully named momtorway network.
# Inplace operation on the file in path directly. Delete all the  
date = param['reference_date']
geomt = 'Polyline'
path = os.path.join(pp['data_out_otshp'], 'motorways_DE_'+date+geomt+'.shp')

process_ot.impute_miss_from_neighbor_line(
        path, 'ref', 'missing', delete_miss=True, delete_isol=True
    )
    





# Get namelist as they are in the shapefile due to field name restrictions 
# (max of 10 characters, "-" and ":" replaced by "_" and "@" by "F_")
#name_list = []
#for item in  in_list:
#    name_list.append(item[0:10].replace(":","_").replace("-","_").replace("@","F_"))
#print("number of duplicates: ", len(name_list)-len(list(set(name_list))))

#cursor = arcpy.da.SearchCursor(pp_in,
#    ['FID', 'Shape', 'OBJECTID', 'F_id', 'bdouble', 'bridge', 'highway','ref'],
#    "\"ref\"='missing'")
#print(cursor.fields)
#
#for row in cursor:
#    print(row)










