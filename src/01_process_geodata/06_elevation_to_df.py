# -*- coding: utf-8 -*-
"""
Runtime: approx 5 hours. 
"""

import arcpy
arcpy.CheckOutExtension('Spatial')
arcpy.CheckOutExtension('3D')
import arcpy.sa
import os
import pandas as pd
import numpy as np
import glob
from arcgis import GIS

import src.process_geodata.functions.f06_elevation_functions as elevation
from project_paths import get_project_paths, get_params
from src.process_numdata.functions.f02_collapse_functions import long_to_wide

arcpy.env.overwriteOutput = True


pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))

sr = arcpy.SpatialReference(param['spacial_reference'])

pp_tempraster_elev = os.path.join(pp['data_temp'], 'elev.tif')
pp_tempraster_idw = os.path.join(pp['data_temp'], 'idw.tif')
pp_tempraster_slope = os.path.join(pp['data_temp'], 'slope.tif')

pp_ot = os.path.join(pp['data_out_otshp'], 'motorways_DE_2019-07-01Polyline.shp')
pp_inframe = os.path.join(pp['data_out_km'], 'frame.shp')
pp_outdf = os.path.join(pp['data_out_df'], 'elevation.csv')

tiles = glob.glob(os.path.join(pp['data_in_elev'], '*.hgt'))
arcpy.MosaicToNewRaster_management(tiles, pp['data_temp'], 'elev.tif', sr, 
                                   '16_BIT_SIGNED', number_of_bands=1)

# Make the points along the frame line.
arcpy.management.GeneratePointsAlongLines(pp_inframe, r"in_memory\frame_pts",
                                          "PERCENTAGE", None, 5, "END_POINTS")

arcpy.analysis.Select(pp_ot, r"in_memory\bridges_tunnel",
                      "bridge = 'yes' Or bridge = 'viaduct' Or tunnel = 'yes' Or tunnel = 'covered' Or tunnel = 'avalanche_protector'")
#arcpy.analysis.buffer(r"in_memory\bridges_tunnel", r"in_memory\small_cutout", 
#                      "36 meters", "full", "round", "none", none, "planar")
arcpy.analysis.Buffer(r"in_memory\bridges_tunnel", r"in_memory\large_cutout",
                      "100 Meters", "FULL", "FLAT", "NONE", None, "PLANAR")
arcpy.analysis.Buffer(pp_ot, r"in_memory\motorways", "90 Meters", "FULL", 
                      "ROUND", "ALL", None, "PLANAR")

arcpy.analysis.Erase(r"in_memory\motorways", r"in_memory\large_cutout", 
                     r"in_memory\motorways_Erase", None)
#arcpy.analysis.Erase(r"in_memory\motorways_Erase", r"in_memory\small_cutout", 
#                     r"in_memory\motorways_Erase1", None)
#arcpy.management.CopyFeatures(r"in_memory\motorways_Erase1", r"in_memory\motorways_Erase")
# Clip raster. 
arcpy.management.Clip(pp_tempraster_elev, 
                      "32290119.1163 5265946.0749 32920173.149 6073534.7723", 
                      r"in_memory\elev_Clip", r"in_memory\motorways_Erase",
                      "32767", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
arcpy.conversion.RasterToPoint(r"in_memory\elev_Clip", r"in_memory\elev_points",
                               "Value")

outIDW = arcpy.sa.Idw(r"in_memory\elev_points", "grid_code", 25, 2, 
                      "FIXED 8000", None)
outIDW.save(pp_tempraster_idw)


arcpy.sa.ExtractValuesToPoints(r"in_memory\frame_pts", pp_tempraster_idw, 
                               r"in_memory\elev_points", "NONE", "VALUE_ONLY")
arcpy.AlterField_management(r"in_memory\elev_points", "RASTERVALU", "elevation")
arcpy.ddd.Slope(pp_tempraster_idw, pp_tempraster_slope, "DEGREE", 1, "PLANAR", 
                "METER")
# Add the slope information to the points with the elevation data. 
arcpy.sa.ExtractValuesToPoints(r"in_memory\elev_points", pp_tempraster_slope, 
                               r"in_memory\inf_points", "NONE", "VALUE_ONLY")
arcpy.AlterField_management(r"in_memory\inf_points", "RASTERVALU", "slope")

# Transform to dataframe.
df = pd.DataFrame.spatial.from_featureclass(r"in_memory\inf_points")
df = df[['ID', 'OID', 'elevation', 'slope']]

# Order and create num to prepare reshape. 
df.sort_values(['ID', 'OID'], inplace=True)
df['num'] = np.arange(len(df.index))
df['temp'] = df.groupby('ID')['num'].transform('min')
df['num'] = df['num']-df['temp']
df['id'] = df['ID']
# Reshape
df.set_index(['id', 'num'], inplace=True)
df = long_to_wide(df, ['elevation', 'slope'])

for i in range(int(len(df.columns)/2)-1):
    df['change_'+str(i)] = df['elevation_'+str(i+1)] - df['elevation_'+str(i)]
change_cols = [x for x in df.columns if x.startswith('change')]
df['down_change'] = df[change_cols].where(df[change_cols]<0).sum(axis=1)
df['up_change'] = df[change_cols].where(df[change_cols]>0).sum(axis=1)
df['net_change'] = df['up_change'] + df['down_change']
## Here starts the feature construction 
df['max_slope'] = df[[x for x in df.columns if x.startswith('slope')]].max(axis=1)
df['mean_slope'] = df[[x for x in df.columns if x.startswith('slope')]].mean(axis=1)
df['elevation'] = df[[x for x in df.columns if x.startswith('elev')]].mean(axis=1)

df = df[['down_change', 'up_change', 'net_change', 'max_slope', 'mean_slope', 
         'elevation']]

df.to_csv(pp_outdf)
