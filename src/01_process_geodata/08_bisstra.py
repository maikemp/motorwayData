# -*- coding: utf-8 -*-
"""
Runtime: approx. 2 Minutes
"""

import arcpy
import os
from geojson import dump
import json 

from src.process_geodata.functions.f01_process_ot_functions import cr_uniform_geojson
from project_paths import get_project_paths, get_params

arcpy.env.overwriteOutput = True


pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))

files = [x for x in os.listdir(pp['data_in_bast']) if x.endswith('.geojson')] 

arcpy.env.overwriteOutput = True
sr = arcpy.SpatialReference(param['spacial_reference'])

for file in files:
    
    path = os.path.join(pp['data_in_bast'], file)
    temp_path = os.path.join(pp['data_temp'], file)
    geojson = json.load(open(path, encoding="utf8"))
    
    properties = []
    geometries = []
    for feature in geojson['features']:
        for key in feature['properties']:
            if key not in properties:
                properties.append(key)
        a = feature['geometry']['type']
        if a == 'LineString':
            a = 'Polyline'
        if a not in geometries:
            geometries.append(a)
        
#    geojson = cr_uniform_geojson(geojson, properties, [])
#    with open(temp_path, 'w') as f:
#        dump(geojson, f)
    
    for geom in geometries:
        
        out_path = os.path.join(pp['data_out_bast'], file[:-8]+geom)
        arcpy.JSONToFeatures_conversion(path, r"in_memory\features", geom)
        # if outpath exists, delete file . 
        
        arcpy.management.Project(r"in_memory\features", out_path, sr)
    
    
    
    

