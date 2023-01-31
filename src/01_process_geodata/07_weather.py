# -*- coding: utf-8 -*-
"""
Runtime: approx. 5 minutes.
"""

import arcpy
arcpy.CheckOutExtension('Spatial')
import arcpy.sa
import os
from arcgis import GIS


import src.process_geodata.functions.f07_weather_functions as weather
from project_paths import get_project_paths

arcpy.env.overwriteOutput = True


pp = get_project_paths()

# Define paths. 
path_frame = os.path.join(pp['data_out_km'], 'frame.shp')
dict_path =  os.path.join(pp['lib'], 'weather_files.json')
file_dict, year_dict = weather.load_file_dict_and_make_year_dict(dict_path)
out_coor_sys, in_coor_sys = weather.get_in_and_out_coor_systems_dwd()

variables = [x for x in os.listdir(pp['data_in_dwd']) 
             if '.' not in x and x!='Wind']

# Define paths for the wind part.
wind_file = 'wind_wdat_geo_10m_BRD_200m.asc' 
p_in_w = os.path.join(pp['data_in_dwd'], 'Wind', wind_file,wind_file,wind_file)
p_out_w = os.path.join(pp['data_out_weather'], 'wind_proj')
p_df_out_w = os.path.join(pp['data_out_df'], 'wind.csv')

# Make frame midpoints as target for all joins. 
arcpy.management.FeatureVerticesToPoints(path_frame,r"in_memory\midpts","MID")


for year in year_dict.keys():
    
    print(year)
    
    pp_df_out = os.path.join(pp['data_out_df'], 'weather_'+year+'.csv')
    i=0
    
    for var in variables:
        
        print(var, year)
        
        path_in = os.path.join(pp['data_in_dwd'], var, year_dict[year][var])
        path_out = os.path.join(pp['data_out_weather'], var+'_'+year+'_proj')
        # res_path = os.path.join(pp['data_out_weather'], var+year+'frame')
        
        # Project raster, save to path_out, find closest raster value to frame
        # midpoints and return a data frame with the variable.
        df_v = weather.project_and_join_nearest_raster_pt(path_in, path_out, 
                           r"in_memory\midpts", in_coor_sys, out_coor_sys, var)
        
        if i == 0: 
            df = df_v[[var]]
        else:
            df = df.join(df_v)
        i+=1
    
    # Save out years specific data frame.
    df.to_csv(pp_df_out)


# Do the same for the wind grid.
df_wind = weather.project_and_join_nearest_raster_pt(
        p_in_w, p_out_w, r"in_memory\midpts", in_coor_sys, out_coor_sys,'wind')
df_wind.to_csv(p_df_out_w)
