# -*- coding: utf-8 -*-
"""

"""

import arcpy
import glob
import os
import time
import numpy as np

from src.process_geodata.functions.f00_general_functions import *

    
  
    
def add_segm_length_to_points(
    pp_buffer, pp_ot, pp_changepts, pp_cutpts,
    pdist, length_div_tolerance, info_freq,
    
):
    """Take a polygon layer buffer_paths, and delete all elements, for which 
    the length of the polyline features in pp_ot inside of the polygon
    element diverts from pdist by more than length_div_tolerance.
    
    """
            
    startpt_dict, endpt_dict, length_dict, tot_length_dict, id_list, \
        unclear_list = _make_dictionaries_and_id_list(
                pp_buffer, pp_ot, pdist, length_div_tolerance, info_freq
            )

    print('Done making dictionaries and lists.')
    
    # Delete all change and cutpoints elements that are not in id_list 
    with arcpy.da.UpdateCursor(pp_changepts, 'ID') as uc:
        for r in uc:
            if r[0] not in id_list:
                uc.deleteRow()
                
    with arcpy.da.UpdateCursor(pp_cutpts, 'ID') as uc:
        for r in uc:
            if r[0] not in id_list:
                uc.deleteRow()
                   
    # Next, add field 'length' to pp_cutpts, and populate via dictionary. 
    for field in ['length', 'starttoend', 'tot_length']:
        arcpy.AddField_management(pp_cutpts, field, 'FLOAT')
    
    j, start = track_time_start()
    # Add Euclidean distance between start and endpoint. 
    with arcpy.da.UpdateCursor(pp_cutpts, ['ID', 'length','starttoend', 'tot_length']) as uc:
        for r in uc:
            # Add length of the first segment f
            if r[0] in list(startpt_dict.keys()):
                j = track_time(j, info_freq, start)
                if startpt_dict[r[0]] in list(length_dict[r[0]].keys()):
                    # segment length.
                    r[1] = length_dict[r[0]][startpt_dict[r[0]]]
                    # Distance between start and endpoint. 
                    r[2] = euclidean(startpt_dict[r[0]], endpt_dict[r[0]])
                    r[3] = tot_length_dict[r[0]]
                    
                else:
                    r[1] = -9999
                uc.updateRow(r)
                    #length_dict[r[0]].pop(startpt_dict[r[0]])
    print('Done adding length to cutpoints..')

    # Next, add field for ot-segment length and  to changepoints.
    arcpy.AddField_management(pp_changepts, 'length', 'FLOAT')
    arcpy.AddField_management(pp_changepts, 'dist', 'FLOAT')
    
    j, start = track_time_start()
    with arcpy.da.UpdateCursor(pp_changepts, ['ID','length','SHAPE@XY', 'dist']) as uc:
        for r in uc:
            if r[0] in list(length_dict.keys()):
                
                # Calculate dist as distance from changept to segment start.
                r[3] = euclidean(r[2], startpt_dict[r[0]]) 
                usedict = length_dict[r[0]]
                j = track_time(j, info_freq, start)
                if usedict != {}:
                    key_list = list(usedict.keys())
                    match_key = key_list[min(
                        range(len(key_list)), 
                        key = lambda i: euclidean(r[2], key_list[i])
                    )]
                    r[1] = usedict[match_key]
                else:
                    r[1] = -9999
                uc.updateRow(r)
    print('Done adding length to changepoints...')    
    return(id_list)    



def _make_dictionaries_and_id_list(
    pp_buffer, pp_ot, pdist, length_div_tolerance, info_freq
):
    '''Create the defined dictionaries and lists...
    
    '''
    start_num = int(arcpy.GetCount_management(pp_buffer)[0])    

    j, start = track_time_start()

    # initiate lists and counters. 
    delete_count = 0
    
    id_list = []
    # id: coordinates
    startpt_dict = {}
    endpt_dict = {}
    # Length_dict for all lengthes of all segments in one buffer. 
    length_dict = {}
    tot_length_dict = {}
    # Keep track of buffers with two full segments in the defined boundaries.
    unclear_list = []

    # Initiate cursor and start loop
    buff_update_c = arcpy.da.UpdateCursor(pp_buffer, ['SHAPE@', 'ID'])

    for buff_row in buff_update_c:
        # keeping track of time
        j = track_time(j, info_freq, start)

        # Select all lines in the polygon
        arcpy.Clip_analysis(
            pp_ot, buff_row[0], r"in_memory\ot_clip"
        )
        # Dissolve lines that share start/end points.
        arcpy.Dissolve_management(
            r"in_memory\ot_clip", r"in_memory\ot_clip_diss",
            unsplit_lines = 'UNSPLIT_LINES'
        )
        # Make a cursor with the lines in the buffer. 
        diss_search_c = arcpy.da.SearchCursor(
            r"in_memory\ot_clip_diss", ['SHAPE@LENGTH', 'SHAPE@']
        )
        # Initiate tracker for if a buffer is kept
        c = 0
        for diss_row in diss_search_c:
            # Check if length is within the defined boundaries. 
            if ((diss_row[0] > (pdist-length_div_tolerance) and \
                 diss_row[0] < (pdist+length_div_tolerance))
            ):
                # if so, add id to id_list and update the startpt_dict end 
                # endpt_dict with the ID and the respective coordinates.
                c += 1
                startpt = diss_row[1].firstPoint
                startpt_dict.update({buff_row[1]: (startpt.X, startpt.Y)})
                endpt = diss_row[1].lastPoint
                endpt_dict.update({buff_row[1]: (endpt.X, endpt.Y)})
                tot_length_dict.update({buff_row[1]: diss_row[0]})
        # If the point was added to the dictionaries, proceed updating the 
        # changept_dict and save startpoints and length to the changedict. 
        if c == 1:
            id_list.append(buff_row[1])
            # Update the length dictionary with initiated sub-dict. 
            length_dict.update({buff_row[1]: {}})
            # Fill empty sub-dict with startpt coordinates and object length.
            with arcpy.da.SearchCursor(
                r"in_memory\ot_clip", ['SHAPE@LENGTH', 'SHAPE@']
            ) as sc:
                for r in sc:
                    startpt = r[1].firstPoint
                    length_dict[buff_row[1]].update({
                        (startpt.X, startpt.Y): r[0]
                    })
        elif c > 1:
            unclear_list.append(buff_row[1])
            buff_update_c.deleteRow()
        else:
            buff_update_c.deleteRow()
            delete_count +=1  
  
    # Cleanup
    arcpy.Delete_management(r"in_memory\ot_clip")
    arcpy.Delete_management(r"in_memory\ot_clip_diss")
    
    # End report
    end_num = int(arcpy.GetCount_management(pp_buffer)[0])
    print('Deleted objects due to false length in buffer:',(start_num-end_num))
    print('should be:', delete_count)
    print('Remaining number of objects:', end_num)
    print('Rows deleted because of multiple full segments:', len(unclear_list))
    
    return(startpt_dict, endpt_dict, length_dict, tot_length_dict, id_list, unclear_list)




def string_to_numeric_attr(featclass, orig_name, new_name):
    """Take a featureclass featclass, for which orig_name is a numeric string.
    Create a new attribute new_name and make it the numeric value of the 
    string.
    
    """
    
    a = []
    fs = arcpy.ListFields(featclass)
    for f in fs:
        a.append(f.name)
    if new_name not in a:
        arcpy.AddField_management(featclass, new_name, "Double")
        arcpy.management.CalculateField(
            featclass, new_name, "int(!{}!)".format(orig_name), "PYTHON3", '', 
            "TEXT"
        )
    else:
        print(new_name, "allready in", featclass)


def make_featureclass_with_fields(
    workspace, name, f_type, sr, field_names, field_type_dict
):
    """Create a new empty feature class with name in workspace, of type f_type,
    with spatial reference sr, that has field_names, the types of which are
    defined in field_type_dict.
    
    """
    
    arcpy.CreateFeatureclass_management(workspace, name, f_type)
    arcpy.DefineProjection_management(workspace + '\\' + name, sr)
    
    for field in field_names:
        arcpy.AddField_management(
            workspace + '\\' + name, field, field_type_dict[field]
        )



def line_attr_to_close_points(lines_path, points_path, out_path, buffer_width,
        keep_buffer = False, buffer_path = r"in_memory\buffer"
    ):
    """Take a feature class of polylines, make a buffer of buff_width widt 
    around them, and then join all attributes from the line to all points that 
    are in the buffer. Save the results to ouput_path and keeps the buffer in 
    r"in_memory\buffer" if keep_buffer = True.
    
    """

    arcpy.CopyFeatures_management(lines_path, r"in_memory\frame")
    arcpy.analysis.Buffer(
        r"in_memory\frame", r"in_memory\buffer", buffer_width, "FULL", 
        "FLAT", "NONE", None, "PLANAR"
    )
    
    # Populate the change points with the IDs from the buffer.
    arcpy.SpatialJoin_analysis(
            points_path, r"in_memory\buffer", out_path, 
            "JOIN_ONE_TO_ONE", match_option="WITHIN" 
        )
    arcpy.Delete_management(r"in_memory\frame")
    
    if (keep_buffer == True and buffer_path != r"in_memory\buffer"):
        arcpy.CopyFeatures_management(r"in_memory\buffer", buffer_path)
        
    if buffer_path != r"in_memory\buffer":
        arcpy.Delete_management(r"in_memory\buffer")




def ot_polylines_to_change_points(
    pp_ot, pp_change_points, keep_fields, type_dict, sr
):
    """Take a polyline layer from pp_ot, take only the start points of them 
    and save a point layer of those start points to pp_change_points. They 
    contain all attributes from the polyline layer, except for 'FID', 'Shape', 
    'OBJECTID', and the once defined in lose_fields. It gets the spacial 
    reference assigned in sr.
    
    """
    
    # Create empty feature class for the change points in memory and assign sr.
    arcpy.CreateFeatureclass_management(r"in_memory", "change_points", "POINT")
    arcpy.DefineProjection_management(r"in_memory\change_points", sr)
    # Add all desired fields and their original type from the input features.
    for field in keep_fields: 
        arcpy.AddField_management(
            r"in_memory\change_points", field, type_dict[field]
        )
        
    # Search for startpoints in the Polylines and insert them and their 
    # relevant attributes to change_points in memory. 
    search_cursor = arcpy.da.SearchCursor(
           pp_ot, ['SHAPE@'] + keep_fields
        )
    insert_cursor = arcpy.da.InsertCursor(
            r"in_memory\change_points", ['SHAPE@'] + keep_fields
        )
    for row in search_cursor:   
        point = row[0].firstPoint 
        insert_cursor.insertRow([point]+list(row[1:]))   
    
    # Save out the change points to the defined path.
    arcpy.CopyFeatures_management(r"in_memory\change_points", pp_change_points)
    
    # Give short information on number of created change points
    point_nr = arcpy.GetCount_management(pp_change_points)[0]
    print('Number of created change points:', point_nr)
    
    # Clean up. 
    arcpy.Delete_management(r"in_memory\change_points")

   
     

def populate_cutpoints(
    pp_in_points, pp_in_lines, pp_frame, pp_pop_points, 
    merge_fields, distance, statistic
):
    """Take a point layer pp_in_points, suck the attributes from pp_in_lines as
    defined in statistic in a distance as defined in distance, for those
    attributes defined in merge_field, and insert them to a new file in
    pp_pop_points. Also add the ID from the underlying frame based on the 
    in_seg_ID and proximity, and save out to pp_pop_points.
    
    """
    
    # Create string of for the spatial join defining all information that 
    # should be adopted from nearby polylines to the cutpoints.
    join_string = ''
    for merge_f in merge_fields:
        join_string = join_string + merge_f + ' ' + statistic + ';'
    join_string = join_string[:-1]
    
    # join features from closest line to startpoints. Automatically looses
    # features that have no line in the defined range.
    
    # Because of a bug in the arcpy JoinFeatures function, can't do the 
    # following in memory. 
    # Take information from pp_in_lines, and join them to pp_in_points, as
    # defined in the  join_string, based on distance
    arcpy.gapro.JoinFeatures(
        pp_in_points, pp_in_lines, pp_pop_points, "JOIN_ONE_TO_ONE", "NEAR", 
        distance, '', None, None, join_string, ''
    )
    
    # Fix field names. Now all start with statistic + '_' and are shortened. 
    new_fields = []
    fields = arcpy.ListFields(pp_pop_points)
    for field in fields:
        new_fields.append(field.name)
    new_fields = new_fields[len(new_fields)-len(merge_fields):]
    
    # This can only be done in memory, so copy to memory.
    arcpy.CopyFeatures_management(pp_pop_points, r"in_memory\pop_cutpoints")
    # Fix field names.
    for i, field in enumerate(new_fields):
        if merge_fields[i][:len(statistic)+2]!=field[len(statistic)+1:9]:
            print('Missaligned Field names!')
        else:
            arcpy.AlterField_management(r"in_memory\pop_cutpoints", field,
                                        merge_fields[i], merge_fields[i])
    # Clean up and replace by file with fixed names.
    for filename in glob.glob(pp_pop_points[:-4]+"*"):
        os.remove(filename)
    arcpy.CopyFeatures_management(r"in_memory\pop_cutpoints", pp_pop_points)
    
    # The final ID must still be transferred. 
    pp_temp = pp_pop_points[:-4] + 'temp' + '.shp'
    # Such final ID from pp_frame into pp_pop_points, and save to pp_temp.
    arcpy.gapro.JoinFeatures(
        pp_pop_points, pp_frame, pp_temp, "JOIN_ONE_TO_ONE", "NEAR", 
        distance, '', None, None, "ID " + statistic, 
        '$join["in_seg_ID"] == $target["in_seg_ID"]'
    )
    
    # Load the new file, with ID to memory. 
    arcpy.Delete_management(r"in_memory\pop_cutpoints")
    arcpy.CopyFeatures_management(pp_temp,r"in_memory\pop_cutpoints")
    # Fix field name 
    arcpy.AlterField_management(
        r"in_memory\pop_cutpoints", statistic + '_ID', 'ID', 'ID'
    )
    # Delete all before created files from the system and replace them with the
    # new fiel with ID and fixed name. 
    for filename in glob.glob(pp_pop_points[:-4]+"*"):
        os.remove(filename)
    arcpy.CopyFeatures_management(r"in_memory\pop_cutpoints", pp_pop_points)
    # Cleanup
    arcpy.Delete_management(r"in_memory\pop_cutpoints")

















