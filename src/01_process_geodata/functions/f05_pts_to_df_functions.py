"""
"""


import arcpy
import os
import numpy as np
import pandas as pd

import src.process_geodata.functions.f00_general_functions as general



def point_dict_to_df(pt_dict, varlist, order_var, id_list):
    '''Take a nested point dictionary with first layer of keys being the 
    IDs, second layer being the values of the (numerical) variable, that
    the elements should be ordered by, and varlist being the variables
    that the nested dictionaries values refer to. Note: varlist must be 
    ordered such that the order_var is the first variable and the other 
    variables are such as the ordering in the lists in the nested dicts.
    
    '''
    
    if order_var != varlist[0]:
        raise ValueError(
                'False ordering of varlist: {} must be the first element.\
                 Check other values, too!'.format(order_var))
        
    pt_dict_long = {}
    for key in pt_dict.keys():
        longlist = []
        for dist in pt_dict[key].keys():
            longlist.extend([dist] + pt_dict[key][dist])
        pt_dict_long.update({key: longlist})
    
    unf_num_dict = dict(zip(id_list, [0]*len(id_list)))    
    b = 0
    for key in pt_dict.keys():
        a = 0
        for dist in pt_dict[key].keys():
            a += 1
        if a > b:
            b = a
        unf_num_dict.update({key: a})

    # Fill up dictionary with nans for missing values.
    for key in list(pt_dict_long.keys()):
         pt_dict_long[key] = pt_dict_long[key] + list(
                np.repeat(np.nan, len(varlist)*b - len(pt_dict_long[key])))
    # Create column names for changepoints...
    unf_cols = []
    for i in range(b):
        for var in varlist:
            unf_cols.append( var+'_'+'0'*(len(str(b-1))-len(str(i))) + str(i))
            
    # Append the elements from pt_dict_long to the unf_num_dict.
    for key in unf_num_dict:
        if key in pt_dict_long.keys():
            pt_dict_long.update({key: [unf_num_dict[key]] + pt_dict_long[key]})
        else:
            pt_dict_long.update({key: [unf_num_dict[key]] + [np.nan]*b*len(varlist)})
        
    # ordered changept dict to dataframe.
    pt_df = pd.DataFrame.from_dict(pt_dict_long, 
                            orient='index', columns=['unf_num']+unf_cols)
    return(pt_df)



def Polygon_to_frame_ID_dict(path_Polygon, path_frame, varlist):
    '''Take a Polygon from path_Polygon, match its attributets to the elements 
    in the frame in path_frame, and create and return a dictionary with of the 
    form {ID: var_list}.
    
    '''
    
    BL_dict = {}

    # BL Polygone verarbeiten... gar nicht nötig
    
    arcpy.SpatialJoin_analysis(path_frame, path_Polygon, r"in_memory\frame_BL")
    
    cursor = arcpy.da.SearchCursor(r"in_memory\frame_BL", ['ID'] + varlist)
    if len(varlist) == 1:
        for r in cursor:
            BL_dict.update({r[0]:r[1]})
    else:
         for r in cursor:
             BL_dict.update({r[0]: list(r[1:])})
        
    return(BL_dict)
    

def Points_to_frame_ID_dict(path_pts, excl_fields_list, order_var):
    '''Take the points in path_pts that have already been preprocessed to 
    contain the IDs and the distance to startpoint. Create a dictionary with
    the desired fields (all except excl_fields_list) of the form:
    {ID: {order_var: varlist}, {order_var: varlist}}
    # Note difference to Polygon_to_frame_ID_dict: more preprocessing required
    with Select_and_label_Points!
    
    '''
    
    varlist = general.make_fieldlist(path_pts, excl_fields_list, 
                                     make_type_dict=False, 
                                     first_fields=['ID'] + [order_var])
    
    point_dict = {}
    
    cursor = arcpy.da.SearchCursor(path_pts, varlist)
    for r in cursor:
        if r[0] in point_dict.keys():
            point_dict[r[0]].update({r[1]: list(r[2:])})
        else:
            point_dict.update({r[0]:{r[1]: list(r[2:])}})
            
    for key in point_dict.keys():
        point_dict[key] = dict(sorted(point_dict[key].items()))
    
    return(point_dict, varlist)
    


def Select_and_label_Points(
    path_in_pts, path_temp_pts, path_out_pts, path_frame, path_buffer_wide, 
    path_buffer_narr, path_side_roads, buff_width_narr, sr, info_freq
):
    '''Select the points that are considered to be on motorways since they are 
    in close proximity to them and closer to them than to any side roads.
    
    path_temp_pts must be in a geodatabase!
    '''
    
    # Select only the point in the narrow buffer.
    arcpy.Clip_analysis(path_in_pts, path_buffer_narr, r"in_memory\clip_pts")
    
    # Check if the points are closer to the frame or to the other side roads.
    arcpy.analysis.Near(r"in_memory\clip_pts", [path_frame, path_side_roads], 
                        buff_width_narr, "LOCATION", "NO_ANGLE", "PLANAR")
    # Select only the points closer to the frame and save to .
    arcpy.Select_analysis(r"in_memory\clip_pts", r"in_memory\temp_pts", 
                          "\"NEAR_FC\"='{}'".format(path_frame))
    
    #arcpy.CopyFeatures_management(path_in_pts, r"in_memory\temp_pts")
    # suck information from frame onto the points and store in path_temp_points.
    arcpy.analysis.SpatialJoin(
             r"in_memory\temp_pts", path_frame, path_in_pts, 
             "JOIN_ONE_TO_MANY", "KEEP_ALL", match_option = "CLOSEST") 
    
    # Make a feature layer with 'NEAR_X' and 'NEAR_Y' as coordinates for cuts
    # and save out to r"in_memory\temp_pts"
    arcpy.management.XYTableToPoint(
            path_in_pts, path_temp_pts, 'NEAR_X', 'NEAR_Y', None, sr)    
    
    # Maybe easier to just copy this into memory space....
    arcpy.CopyFeatures_management(path_temp_pts, r"in_memory\temp_pts")
    # Add field for distance of start to  accident along frame line to accidents.
    arcpy.AddField_management(r"in_memory\temp_pts", 'from_start', 'FLOAT')
    
    acc_c = arcpy.da.UpdateCursor(r"in_memory\temp_pts", 
                                  ['SHAPE@', 'ID', 'SHAPE@XY', 'from_start'])
    arcpy.CopyFeatures_management(path_frame, r"in_memory\frame")
    j, start = general.track_time_start()
    
    for acc in acc_c:
        j = general.track_time(j, info_freq, start)
        # Select frame segment with the respective ID.
        arcpy.MakeFeatureLayer_management(r"in_memory\frame", 
                r"in_memory\frame_id", "\"ID\" = '{}'".format(acc[1]))
        # Split frame segement at accident point.
        arcpy.SplitLineAtPoint_management(r"in_memory\frame_id", 
                                          acc[0], r"in_memory\split_frame")
        mini_dict = {}
        with arcpy.da.SearchCursor(r"in_memory\split_frame", ['SHAPE@']) as split_c:
            for split in split_c:
                startpt = split[0].firstPoint
                startpoint = (startpt.X, startpt.Y)
                length = split[0].length
                mini_dict.update({startpoint: length})
        
        # Select the key that is further away to the accident point.
        keys = list(mini_dict.keys())
        # Distance from the accident point
        dist_list = [general.euclidean(acc[2], key) for key in keys]
        key = keys[dist_list.index(max(dist_list))]
        
        acc[3] = mini_dict[key]
        acc_c.updateRow(acc)
    
    arcpy.CopyFeatures_management(r"in_memory\temp_pts", path_out_pts)
    # Cleanup
    for memo in [r"in_memory\clip_pts", r"in_memory\temp_pts",
                 r"in_memory\frame", r"in_memory\frame_id",
                 r"in_memory\split_frame"]:
        arcpy.Delete_management(memo)

        
        

def make_side_roads(pp_buff, pp_out, pp_temp_BL, BL_path_list, filesname, sr): 
    '''Take the paths in BL_path_list, and thereof the files with name
    filesname, project them into the spatial reference defined in sr and 
    save out intermediary file to pp_temp_BL. Then take the file and cut out
    only part inside the buffer in pp_buff and concatenate for all BLs, and
    save out to pp_temp.
    
    '''
    
    i=0 
    for BL in BL_path_list:
        print(i)
        path = os.path.join(BL, filesname)
        arcpy.Select_analysis(path, r"in_memory\BL", "fclass <> 'motorway'")
        arcpy.management.Project(r"in_memory\BL", pp_temp_BL, sr)
        
        arcpy.Clip_analysis(pp_temp_BL, pp_buff, r"in_memory\clip")
        
        if i == 0:
            arcpy.CopyFeatures_management(r"in_memory\clip", pp_out)
        else:
            arcpy.Merge_management([pp_out, r"in_memory\clip"], r"in_memory\temp")
            arcpy.CopyFeatures_management(r"in_memory\temp", pp_out)
        i+=1
        




