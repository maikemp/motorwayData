# -*- coding: utf-8 -*-
"""
Runtime: approx. 5 Minutes
"""
import arcpy
import os 
import pandas as pd

import src.process_geodata.functions.f02_cr_km_functions as cr_km
from project_paths import get_project_paths, get_params

arcpy.env.overwriteOutput = True


pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))

date = param['reference_date']
pdist = param['section_length_m']

min_length = 1.6
buffer_width =   '35 Meters'

pp_in = os.path.join(pp['data_out_otshp'], 'motorways_DE_'+date+'Polyline.shp')
pp_temp =  os.path.join(pp['data_temp'], 'motorways_DE'+date+'dissolved.shp')
pp_out = os.path.join(pp['data_out_km'], 'frame.shp')
pp_out_points = os.path.join(pp['data_out_km'], 'cutpoints.shp')

# Create some statistics on motorways:
len_by_A = cr_km.create_length_by_ref_dict(date, "ref", pp_in)
len_by_maxspeed = cr_km.create_length_by_ref_dict(date, "maxspeed", pp_in)

# Dissolve lines from pp_in by 'ref' and remove any cutbacks with an angle <45.
# and save out to pp_temp.
cr_km.dissolve_lines_no_cutbacks(pp_in, pp_temp, dissolve_field='ref', 
                                 angle=45)


# Delete all segments that are shorten than min_length, assign a permanent ID,
# make a field with the ID of parallel road, and keep a dictionary thereof as 
# id_dict, and make a field that contains a 1 to one road and a 0 to all 
# parallel running roads.
id_dict, length_dict = cr_km.del_short_and_assign_id_of_opposing_road(
        min_length, pp_temp, buffer_width
    )

# Cut the polylines into pdist long pieces, and save them to  pp_out, save
# the respective cutpoints to pp_out_points. Change number for other length.
cr_km.cut_sections(pp_temp, pp_out, pp_out_points, 
                   param['spacial_reference'], pdist)

# Create a unique ID for each 1 km Segment. 
# Delete objects that are not equal 1 km. Count how many are lost.
cr_km.make_id_and_remove_short_sec(pp_out, pdist)

# Create numeric dataframe with all and only the attributes (columns) of frame.
arcpy.conversion.TableToTable(pp_out, pp['data_out_df'], "frame.csv")

df_f = pd.read_csv(os.path.join(pp['data_out_df'], "frame.csv"), index_col=6, sep=";")
df_f.to_csv(os.path.join(pp['data_out_df'], "frame.csv"))
