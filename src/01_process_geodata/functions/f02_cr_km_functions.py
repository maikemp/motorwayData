"""

"""
import re
import arcpy



def make_id_and_remove_short_sec(path, pdist):
    
    arcpy.AddField_management(path, 'ID', 'TEXT')
    
    delete_count = 0
    lost_length = 0
    with arcpy.da.UpdateCursor(
        path, ['SHAPE@LENGTH', 'ref', 'seg_ID','in_seg_ID', 'ID']
    ) as cursor:
        for row in cursor:
            # Delete too short sections. 
            if round(float(row[0]),2) != pdist:
                cursor.deleteRow()
                delete_count +=1
                lost_length += row[0]/1000
            # Make individual ID
            else:
                row[4] =  row[1][2:] + '0'*(3-len(row[2]))+row[2] +\
                            '0'*(4-len(row[3]))+row[3]
                cursor.updateRow(row)
                # Check if ID has a correct format.
                if len(row[4]) != 10:
                    print('Caution! Diverting ID length:',len(row[4]))
    
    print(str(delete_count), 'segments deleted, whereby', str(lost_length), 
      'km were lost')


def cut_sections(pp_in, pp_out, pp_out_points, spacial_ref, pnt_dist=1000):
    """Take Polylines in pp_in, cut them into pieces of length pnt_dict, 
    as defined in spacial_ref, that has to be the same as in the input file,
    and save out the cut features in pp_out, and the respective cut points in 
    pp_out_points. 
    
    """
    
    # Create containers for the cut lines and the break pionts, assign
    # spacial referenfe and add the required field names.             
    arcpy.CreateFeatureclass_management(r"in_memory", "break_points", "POINT")
    arcpy.CreateFeatureclass_management(r"in_memory", "cut_lines", "POLYLINE")
    sr = arcpy.SpatialReference(spacial_ref)
    arcpy.DefineProjection_management(r"in_memory\break_points", sr)
    arcpy.DefineProjection_management(r"in_memory\cut_lines", sr)
    # Container for the in segment ID.
    arcpy.AddField_management(r"in_memory\break_points", 'in_seg_ID', 'TEXT')
    arcpy.AddField_management(r"in_memory\cut_lines", 'in_seg_ID', 'TEXT')
    
    # Get list of field names that should be preserved from the input path and
    # add them to the cut_lines feature class.
    descr = arcpy.Describe(pp_in)
    fieldnames = [field.name for field in descr.fields \
                  if field.name not in ['FID','Shape']]
    fieldnames.index('seg_ID')
    for field in fieldnames:
        arcpy.AddField_management(r"in_memory\cut_lines", field, 'TEXT')
    
    
    search_cursor = arcpy.da.SearchCursor(pp_in, ['SHAPE@'] + fieldnames)  
    insert_cursor = arcpy.da.InsertCursor(
        r"in_memory\break_points", ['SHAPE@'] + ['in_seg_ID']
    )
    insert_cursor_output = arcpy.da.InsertCursor(
        r"in_memory\cut_lines", ['SHAPE@'] + fieldnames + ['in_seg_ID']
    )      
    
    # Search in the in input path. 
    for row in search_cursor: 
        # Take only one element from the search cursor and search only in it.
        arcpy.Select_analysis(
                pp_in, r"in_memory\row_i", 
                "\"seg_ID\"='{}'".format(str(row[fieldnames.index('seg_ID')+1]))
            )
        with arcpy.da.SearchCursor(
            r"in_memory\row_i", ['SHAPE@']+fieldnames
        ) as temp_cursor:
            for temp_row in temp_cursor:
                # Make point feature class only for this row and make it an
                # input cursor.
                arcpy.CreateFeatureclass_management(
                    r"in_memory", "temp_points", "POINT"
                )
                temp_insert_c = arcpy.da.InsertCursor(
                                    r"in_memory\temp_points", 'SHAPE@'
                                )
                i = 0
                point_dict = {}
                # Set a point all pnt_dist meters and save to the temporary 
                # point feature class and the break_point features.
                for dist in range(0, int(temp_row[0].length), pnt_dist):  
                    point = temp_row[0].positionAlongLine(dist).firstPoint  
                    insert_cursor.insertRow(tuple((point,str(i))))
                    temp_insert_c.insertRow([point])
                    
                    # Save shortended and rounded coordinates in a dictionary
                    # with a counter thats 0 for first cut point etc..
                    point_dict.update(
                        {(round(point.X/100), round(point.Y/100)): str(i)}
                    )
                    
                    i += 1
                # Split the row at the temp points
                arcpy.SplitLineAtPoint_management(
                    r"in_memory\row_i", r"in_memory\temp_points", 
                    r"in_memory\cut_row","0.05 Meters"
                )
                # In the cur roow, search for the start point and insert
                # the consecutive counter for elements in segment.
                with arcpy.da.SearchCursor(
                    r"in_memory\cut_row",['SHAPE@']+fieldnames
                ) as s_cursor:
                    for s_row in s_cursor:
                        start_point = (
                            round(s_row[0].firstPoint.X/100), 
                            round(s_row[0].firstPoint.Y/100)
                        )
                        # With the in_sec counter, insert row to output class.
                        insert_cursor_output.insertRow(
                            s_row + tuple((point_dict[start_point],))
                        )
    # Save out cut lines and break points.
    arcpy.CopyFeatures_management(r"in_memory\cut_lines", pp_out)       
    arcpy.CopyFeatures_management(r"in_memory\break_points", pp_out_points)
    
    # Clean up. 
    mefi = [r"in_memory\break_points",r"in_memory\cut_lines",
            r"in_memory\row_i", r"in_memory\temp_points", r"in_memory\cut_row"]
    for memory_file in mefi:
        arcpy.Delete_management(memory_file)



def del_short_and_assign_id_of_opposing_road(min_length, path, buffer_width):
    """Delete all segments in path that are shorter than min_length, and 
    insert a par_seg_ID that points to the seg_ID of the parallel running road,
    in the marging of buffer_width. Also, assign one_direct, that gives a 1 to
    one road and a 0 to all parallel running roads.
    
    """
    
    # Keep track of deleted elements and their respective length
    delete_count = 0
    lost_length = 0
    # Initiate dictionary with unique ID for each segment and their length.
    length_dict = {}
    # Prepare making 'FID' permanent as 'seg_ID'.
    arcpy.AddField_management(path, 'seg_ID', "TEXT")
    
    with arcpy.da.UpdateCursor(path,['SHAPE@LENGTH','FID','seg_ID']) as cursor:    
        for row in cursor:
            # Assign FID as permanent ID.
            row[2] = row[1]
            cursor.updateRow(row)
            # Delete too short segments.
            if int(row[0])/1000 < min_length:
                cursor.deleteRow()
                delete_count +=1
                lost_length += row[0]/1000
            else:
                # Update length_dict with remaining segments
                length_dict.update({row[2]: row[0]/1000})

    print(
        'Deleted because shorter than',min_length,'km:',delete_count,'Segments'
    )
    print('Lost length:', round(lost_length,2), 'km.')
    
    # Create paral_ID and that points to the ID of the parallel running 
    # motorway, and assing one_direct that is 1 for one road, and 0 for all
    # its parallel running roads.
    return_dict = _point_to_parallel_mtw(path, buffer_width, length_dict)
    
    return(return_dict, length_dict)
    

def _point_to_parallel_mtw(path, buffer_width, length_dict):
    """Insert a par_seg_ID that points to the ID of the parallel running road,
    in the marging of buffer_width. Also, assign one_direct, that gives a 1 to
    one road and a 0 to all parallel running roads. 
    
    """
    
    # Create a buffer layer of buffer_width.
    arcpy.Buffer_analysis(path, r"in_memory\buffer", buffer_width)
    
    buffer_cursor = arcpy.da.SearchCursor(
            r"in_memory\buffer", ["seg_ID", "SHAPE@LENGTH"]
        )
    buffer_dict = {}
    
    for rowbuff in buffer_cursor:
        # Make a layer that consists only of the one buffer.
        arcpy.Select_analysis(
            r"in_memory\buffer", r"in_memory\row_l", 
            "\"seg_ID\"='{}'".format(str(rowbuff[0]))
        )
        temp_dict = {}
        # Select only the lines that are within the selected buffer.
        arcpy.Clip_analysis(path, r"in_memory\row_l", r"in_memory\clip") 
        with arcpy.da.SearchCursor(
                r"in_memory\clip", ["seg_ID", "SHAPE@LENGTH"]
        ) as line_cursor:
            for rowline in line_cursor:
                if rowline[1]>0:
                    temp_dict.update({rowline[0]: rowline[1]/1000})
                    # print(rowline[0], rowline[1]/1000)
        temp_dict.pop(rowbuff[0])
        buffer_dict.update({rowbuff[0]:temp_dict})
        line_cursor.reset()
        arcpy.Delete_management(r"in_memory\clip")
        arcpy.Delete_management(r"in_memory\row_l")
    buffer_cursor.reset()
    
    # Add field that points to the parallel ID.
    arcpy.AddField_management(path, 'par_seg_ID', "TEXT")
    parallel_dict = {}
    with arcpy.da.UpdateCursor(path, ["seg_ID","par_seg_ID"]) as cursor:
        for row in cursor:
            row[1] = max(buffer_dict[row[0]], key=buffer_dict[row[0]].get)
            cursor.updateRow(row)
            parallel_dict.update({row[0]:row[1]})
    
    arcpy.AddField_management(path, 'one_direct', "TEXT")
    track_dict = {}
    return_dict = {}
    with arcpy.da.UpdateCursor(path, ["seg_ID","par_seg_ID","one_direct"]) as cursor:
        for row in cursor:
            if (row[0] not in list(track_dict.values()) and 
                list(parallel_dict.values()).count(row[1])==1):
                track_dict.update({row[0]: row[1]})
                row[2] = 0
                cursor.updateRow(row)
            else:
                row[2] = 1
                cursor.updateRow(row)
            return_dict.update({row[0]:row[1]})
        
    # Clean memory.
    mefi = [r"in_memory\buffer", r"in_memory\row_l", r"in_memory\clip"]
    for memory_file in mefi:
        arcpy.Delete_management(memory_file)
    
    return(return_dict)
    


def dissolve_lines_no_cutbacks(pp_in, pp_out, dissolve_field='ref', angle=45):
    """Take the polyline shapefile from pp_in, change 'ref' field to uniform 3 
    digit form, remove all cutbacks that have an angle smaller than defined in 
    angle, dissolve by dissolve_field and save out to pp_out.
    
    """
    # By changing ref = 'A 1' to ref = 'A 001' etc., the ordering will be 
    # right in the dissolved data.
    cursor = arcpy.da.UpdateCursor(pp_in, ["ref"])
    for row in cursor:
        row[0] = row[0].split(' ')[0] + ' ' + \
                '0'*(5-len(row[0]))+row[0].split(' ')[1]
        cursor.updateRow(row)
    
    arcpy.management.Dissolve(
        pp_in, r"in_memory\dissolved_no_cutbacks", 
        dissolve_field=dissolve_field, unsplit_lines='UNSPLIT_LINES'
    )
    # Remove cutbacks from the before created dissolved layer.
    arcpy.CheckOutExtension("Foundation")
    arcpy.topographic.RemoveCutbackVertices(
        r"in_memory\dissolved_no_cutbacks", angle, "ALL", "SKIP_COINCIDENT"
    )
    # Create points on vertices.
    arcpy.management.FeatureVerticesToPoints(
        r"in_memory\dissolved_no_cutbacks", r"in_memory\points_no_cutbacks", 
        "ALL"
    )
    
    # Create standard dissolved feature layer.
    arcpy.management.Dissolve(
        pp_in, r"in_memory\dissolved_cutbacks", 
        dissolve_field=dissolve_field, unsplit_lines='UNSPLIT_LINES'
    )
    # Create points on vertices.
    arcpy.management.FeatureVerticesToPoints(
        r"in_memory\dissolved_cutbacks",r"in_memory\points_cutbacks", "ALL"
    )
    
    # Keep only the points that are on the  cutbacks.
    arcpy.analysis.SymDiff(
        r"in_memory\points_cutbacks",r"in_memory\points_no_cutbacks",
        r"in_memory\points_only_cutbacks", "ALL", None
    )
    # Create a small buffer around these point.
    arcpy.analysis.Buffer(
        r"in_memory\points_only_cutbacks", r"in_memory\cutback_buffer", 
        "1 Meters","FULL", "ROUND", "NONE", None, "PLANAR"
    )
    # Cutout the roads inside the buffer.
    arcpy.Clip_analysis(
        r"in_memory\dissolved_cutbacks",
        r"in_memory\cutback_buffer", r"in_memory\buffer_roads"
    )
    # Remove the  roads in the buffer from the input file.
    arcpy.analysis.SymDiff(
        pp_in, r"in_memory\buffer_roads", r"in_memory\clean_roads", "ALL", None
    )
    # Dissolve the cleaned input file and save out to pp_out.
    arcpy.management.Dissolve(
            r"in_memory\clean_roads", pp_out, 
            dissolve_field=dissolve_field, unsplit_lines='UNSPLIT_LINES'
        )
    # Delete all the  files in memory.
    mefi = [r"in_memory\dissolved_no_cutbacks",r"in_memory\points_no_cutbacks",
            r"in_memory\dissolved_cutbacks", r"in_memory\points_cutbacks",
            r"in_memory\points_only_cutbacks", r"in_memory\cutback_buffer",
            r"in_memory\buffer_roads", r"in_memory\clean_roads"]
    for memory_file in mefi:
        arcpy.Delete_management(memory_file)



def create_length_by_ref_dict(date, field, path):
    """Measures the length of the Polylines grouped by field in km and outputs 
    a dictionary that contains the unique field values and the length of the 
    objects for which this field value is assigned.
    
    """
    arcpy.management.Dissolve(
        path, r"in_memory\features"+date.replace("-","")+"_dissolved", 
        dissolve_field=field
    )
    cursor = arcpy.da.SearchCursor(
        r"in_memory\features"+date.replace("-","")+"_dissolved", 
        ["SHAPE@LENGTH",field]
    )
    
    length_dict = {}
    
    for row in cursor:
        # print('length of ', row[1], ' : ', str(int(row[0])/1000), ' km')
        if field == 'ref':
            length_dict.update({
                row[1].split(' ')[0]+' '+'0'*(5-len(row[1]))+row[1].split(' ')[1]:
                    int(row[0])/1000
                })    
        else:
            length_dict.update({row[1]: int(row[0])/1000})
            
    cursor.reset()
    
    # For the field 'ref' make a more convenient ordering of the output dict.
    if field == 'ref':    
        length_dict_sorted = {}
    
        for key in sorted(length_dict.keys()):
            length_dict_sorted.update({key: length_dict[key]})
        for key in list(length_dict_sorted.keys()):
            mini_dict = {key: length_dict_sorted[key]}
            length_dict_sorted.pop(key)
            length_dict_sorted.update({re.sub(' 0*', ' ', key):mini_dict[key]})
        length_dict = length_dict_sorted
    
    arcpy.Delete_management(r"in_memory\features"+date.replace("-","")+"_dissolved")
    
    return(length_dict)



