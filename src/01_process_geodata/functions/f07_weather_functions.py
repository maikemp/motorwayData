"""
"""
import json
import re
import os
import arcpy
import shutil
import pandas as pd
arcpy.CheckOutExtension('Spatial')
from arcpy.sa import *


def project_and_join_nearest_raster_pt(path_in, path_out, frame_pts, 
                                       in_coor_sys, out_coor_sys, var): 
    
    # Assign the right projection to dwd files.
    arcpy.management.DefineProjection(path_in, in_coor_sys)
    # Project the raster into the used projection.        
    arcpy.management.ProjectRaster(path_in, path_out, out_coor_sys, 
            "NEAREST", "1000 1000", "DHDN_To_ETRS_1989_8_NTv2", None, 
            in_coor_sys, "NO_VERTICAL")
    
    # Transform grid to points. 
    arcpy.conversion.RasterToPoint(path_out, r"in_memory\raster_pts")
    arcpy.AlterField_management(r"in_memory\raster_pts", "grid_code", var)
    # Spatial join 
    arcpy.analysis.SpatialJoin(frame_pts, r"in_memory\raster_pts", 
                               r"in_memory\joinpts", "JOIN_ONE_TO_ONE", 
                               "KEEP_ALL", match_option="CLOSEST", 
                               search_radius="1500 Meters")
    
    df = pd.DataFrame.spatial.from_featureclass(r"in_memory\joinpts")
    df = df[['ID', var]]
    df.set_index('ID', inplace=True)
    
    df.index.name = None
    
    return(df)


def get_in_and_out_coor_systems_dwd():
    '''in_coor_sys info comes from DWD website: 
    https://opendata.dwd.de/climate_environment/CDC/help/gk3.prj
    '''
    
    out_coor_sys = '''PROJCS['ETRS_1989_UTM_Zone_32N_8stellen',
                      GEOGCS['GCS_ETRS_1989',
                          DATUM['D_ETRS_1989',
                          SPHEROID['GRS_1980',6378137.0,298.257222101]],
                      PRIMEM['Greenwich',0.0],
                      UNIT['Degree',0.0174532925199433]],
                      PROJECTION['Transverse_Mercator'],
                      PARAMETER['False_Easting',32500000.0],
                      PARAMETER['False_Northing',0.0],
                      PARAMETER['Central_Meridian',9.0],
                      PARAMETER['Scale_Factor',0.9996],
                      PARAMETER['Latitude_Of_Origin',0.0],
                      UNIT['Meter',1.0]]'''
    in_coor_sys = '''PROJCS["DHDN_3_Degree_Gauss_Zone_3",
                     GEOGCS["GCS_Deutsches_Hauptdreiecksnetz",
                         DATUM["D_Deutsches_Hauptdreiecksnetz",
                             SPHEROID["Bessel_1841",6377397.155,299.1528128]],
                         PRIMEM["Greenwich",0.0],
                         UNIT["Degree",0.0174532925199433]],
                     PROJECTION["Gauss_Kruger"],
                     PARAMETER["False_Easting",3500000.0],
                     PARAMETER["False_Northing",0.0],
                     PARAMETER["Central_Meridian",9.0],
                     PARAMETER["Scale_Factor",1.0],
                     PARAMETER["Latitude_Of_Origin",0.0],
                     UNIT["Meter",1.0]]'''
    return(out_coor_sys, in_coor_sys)


def load_file_dict_and_make_year_dict(file_dict_path):
    
    file_dict = json.load(open(file_dict_path))
    pattern = re.compile(r'20\d{2}')
    year_dict = {}
    for key in file_dict:
        for file in file_dict[key]:
            year  = pattern.search(file)[0]
            if year in year_dict.keys():
                year_dict[year].update({key: file})
            else:
                year_dict.update({year: {key: file}})
            year_dict.update()
    return(file_dict, year_dict)


def copy_files_and_make_dict(variables, dict_path, pp):

    file_names = {}
    
    for folder in variables:
        years = [x for x in os.listdir(os.path.join(pp['data_in_dwd'], \
                                folder)) if '.pdf' not in x and '.gz' not in x]
        a = []
        for year in years:
            path = os.path.join(pp['data_in_dwd'], folder, year)
            subfiles = os.listdir(path)
            for file in subfiles:
                shutil.copy(os.path.join(path, file), 
                            os.path.join(pp['data_in_dwd'], folder))
                a.append(file)
        file_names.update({folder:a})
        
        json.dump(file_names, open(dict_path, 'w')) 



