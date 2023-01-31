# -*- coding: utf-8 -*-
"""

"""


import os
import json
import sys

#import arcpy
from datetime import datetime
print(datetime.now().strftime("%H:%M:%S"))



def get_project_paths():

    pp = {}
    
    # DEFINE PATHS HERE.
    # Root of the directory the Repo is downloaded to. 
    pp['root'] = r"..."
    # Directory, the data folders (data_in, data_out and data_temp) are in.
    pp['data'] = r"..."
    
    # Structure data storage.
    pp['data_in'] = os.path.join(pp['data'], 'data_in')

    pp['data_in_gf'] = os.path.join(pp['data_in'], 'OSM_shape_BL')
    pp['data_in_unfall'] = os.path.join(pp['data_in'], 'unfallatlas')
    pp['data_in_ot'] = os.path.join(pp['data_in'], 'overpass_turbo')
    pp['data_in_BL'] = os.path.join(pp['data_in'],'BKG','vg250-ew_ebenen_1231')
    pp['data_in_elev'] = os.path.join(pp['data_in'], 'Nasa_srtm30_tiles')
    pp['data_in_count'] = os.path.join(pp['data_in'], 'Zaehlstellen')
    pp['data_in_bast'] = os.path.join(pp['data_in'], 'BISStra')
    pp['data_in_dwd'] = os.path.join(pp['data_in'], 'DWD')
    pp['data_in_zeb'] = os.path.join(pp['data_in'], 'ZEB')
    pp['data_in_inkar'] = os.path.join(pp['data_in'], 'INKAR')

    pp['data_temp'] = os.path.join(pp['data'],"data_temp")
    pp['data_temp_gdb'] = os.path.join(pp['data_temp'], "temp.gdb")
    
    pp['data_out'] = os.path.join(pp['data'], 'data_out')
    
    pp['data_out_km'] = os.path.join(pp['data_out'], 'km_sections')
    pp['data_out_otshp'] = os.path.join(pp['data_out'], 'overpass_shapefile')
    pp['data_out_otsecpt'] = os.path.join(pp['data_out'], 'ot_on_sec_points')
    pp['data_out_df'] = os.path.join(pp['data_out'], 'data_frames')
    pp['data_out_BL'] = os.path.join(pp['data_out'], 'BL')
    pp['data_out_unf'] = os.path.join(pp['data_out'], 'unfaelle')
    pp['data_out_eleva'] = os.path.join(pp['data_out'], 'elevation')
    pp['data_out_bast'] = os.path.join(pp['data_out'], 'bast_netz')
    pp['data_out_links'] = os.path.join(pp['data_out'], 'links')
    pp['data_out_weather'] = os.path.join(pp['data_out'], 'weather.gdb')
    pp['data_out_zeb'] = os.path.join(pp['data_out'], 'ZEB')
    pp['data_out_inkar'] = os.path.join(pp['data_out'], 'INKAR')
    
    # Outputs for analysis data
    pp['data_out_a'] = os.path.join(pp['data_out'], 'analysis')
    pp['data_out_a_exp'] = os.path.join(pp['data_out_a'], 'exp_results')
    pp['data_out_a_oob'] = os.path.join(pp['data_out_a'], 'oob_predictions')
    pp['data_out_a_grf'] = os.path.join(pp['data_out_a'], 'grf_results')
    
    # Structure workspace.
    pp['process_geo'] = os.path.join(pp['root'], 'src', '01_process_geodata')
    pp['geo_functions'] = os.path.join(pp['process_geo'], 'functions')
    
    pp['process_num'] = os.path.join(pp['root'], 'src', '02_process_numdata')
    pp['num_functions'] = os.path.join(pp['process_num'], 'functions')
    
    pp['fin_data'] = os.path.join(pp['root'], 'out','data')
    pp['analysis_data'] = os.path.join(pp['root'], 'analysis_data')
    
    # Outputs
    pp['out'] = os.path.join(pp['root'], 'out')
    pp['out_tables'] = os.path.join(pp['out'], 'tables')
    pp['out_graphics'] = os.path.join(pp['out'], 'graphics')
    pp['out_regressions'] = os.path.join(pp['out'], 'regression_results')
    
    pp['lib'] = os.path.join(pp['root'], 'Libraries')
    pp['lit_cs'] = os.path.join(pp['root'], 'Literature', 'data', 'Unfallstatistik')    
    
    return(pp)

pp = get_project_paths()

def get_params(path):
    
    param = json.load(open(path))
    
    return(param)

# Set workspace and load required parameters.
sys.path.insert(1, pp['root'])


