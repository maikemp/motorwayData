"""
Contains all functions required for the processing of the overpass turbo data.

"""
import re
import arcpy


def rename_fields(geojson, rename_dict):
    '''Take the geojson file and change all keys in the dicitonaries
    in the list in geojson['features'] by replacing the string-sequences
    as defined in the keys of rename_dict by their corresponding values. 
    '''
    
    new_names = list(geojson['features'][0]['properties'].keys())
    for k, v in rename_dict.items():
        new_names = [x.replace(k, v) for x in new_names]
        
    field_renames = dict(zip(
            list(geojson['features'][0]['properties'].keys()), 
            new_names))
    
    for k, v in list(field_renames.items()):
        if k==v:
            field_renames.pop(k)
    
    for f in geojson['features']:
        for k, v in field_renames.items():
            f['properties'].update({v: f['properties'][k]})
            f['properties'].pop(k)
            
    return(geojson)


def impute_miss_from_neighbor_line(
        path, field, miss_val='missing', delete_miss=False, delete_isol=False
    ):
    """For a polyline, replace the miss_val of a field from a linked polyline's 
    field entry, if this entry is not equal to miss_val. Repeat until the 
    number of missing values does not change anymore. When delete = True, all 
    the remaining features with missings on this field are deleted. If 
    delete_isol = True, all features that do not share any start or end point
    with any other feature are deleted.
    
    """
    
    # Initialize counter list for the number of missing values.    
    count_miss_list = [0]
    count_miss = int(arcpy.GetCount_management(path)[0])
    count_miss_list.append(count_miss)
    
    # Initialize dictionary of end/start points as keys and field as value.
    point_dict = {}
    
    # Repeat operations until the missing values do not change anymore.
    while (count_miss_list[-1]!=count_miss_list[-2]):
        
        # For all features that are not missing, update the dictionary with the 
        # start/end points as values and field as value.
        cursor = arcpy.da.SearchCursor(
            path, ["OID@","SHAPE@",field],"\"{}\"<>'{}'".format(field,miss_val)
        )
        count_miss = int(arcpy.GetCount_management(path)[0])

        for row in cursor:
            
            count_miss = count_miss - 1
            
            startpt = row[1].firstPoint
            endpt = row[1].lastPoint
        
            point_dict.update({(startpt.X, startpt.Y): row[2]})
            point_dict.update({(endpt.X, endpt.Y): row[2]})
            
        cursor.reset()
        
        # For all missing values search in point_dict for the start and end
        # point and impute the field value from it.
        cursor = arcpy.da.UpdateCursor(
            path, ["OID@","SHAPE@",field], "\"{}\"='{}'".format(field,miss_val)
        )
    
        for row in cursor:
            
            startpt = row[1].firstPoint
            endpt = row[1].lastPoint
            
            startpoint = (startpt.X, startpt.Y)
            endpoint = (endpt.X, endpt.Y)
            
            
            if startpoint in list(point_dict.keys()):
                row[2] = point_dict[startpoint]
                count_miss = count_miss - 1
                cursor.updateRow(row)
            
            elif endpoint in list(point_dict.keys()):
                row[2] = point_dict[endpoint]
                count_miss = count_miss - 1
                cursor.updateRow(row)
                
        cursor.reset()
        
        # Update count_miss_list by the now counted value.
        count_miss_list.append(count_miss)
    
    print('Number of remaining missings in ref after imputing from neighbors:',
          str(count_miss_list[-1])
    )
    
    if delete_miss == True:    
        _delete_missings(path, field, miss_val)

    if delete_isol == True:    
        _delete_isolated(path)



def _delete_isolated(path):
    """Delete all features in path that do not share an end or start point with
    any other feature.
    
    """
    cursor = arcpy.da.UpdateCursor(path, ["F_id","SHAPE@"])        
    count_delete = 0
    point_dict = {}
    
    for row in cursor:
        
        startpt = row[1].firstPoint
        endpt = row[1].lastPoint
        
        startpoint = (startpt.X, startpt.Y)
        endpoint = (endpt.X, endpt.Y) 
        
        point_dict.update({str(row[0]) + '_start': startpoint})
        point_dict.update({str(row[0]) + '_end': endpoint})
    
    cursor.reset()
    
    for row in cursor:
        use_dict = point_dict.copy()
        use_dict.pop(str(row[0]) + '_start')
        use_dict.pop(str(row[0]) + '_end')
    
        if (point_dict[str(row[0])+'_start'] not in list(use_dict.values()) and
            point_dict[str(row[0]) + '_end'] not in list(use_dict.values())):
            cursor.deleteRow()
            count_delete += 1
    
    cursor.reset()
    print('Number of deleted isolated features: ', str(count_delete))



def _delete_missings(path, field, miss_val):
    """Delete all features in path that still have miss_val in field. 
    
    """  
    cursor = arcpy.da.UpdateCursor(
            path, [field], "\"{}\"='{}'".format(field,miss_val)
        )
        
    # Keep track of number of deleted values for safety reasons.
    count_delete = 0
    
    # Delete all the features that still have missing values.
    for row in cursor:
        cursor.deleteRow()
        count_delete = count_delete + 1
    
    cursor.reset()
    print('Number of deleted features with missings in ref:',str(count_delete)) 




def replace_if_missing(geojson, replace_dict, miss_val = "missing"):
    """For all features, if target feature defined in the key of the 
    replace_dict has value miss_val, insert the value from the  donor field
    defined in the value of the replace_dict.
    
    """
    
    for field in list(replace_dict.keys()):
        for feature in geojson['features']:
            if feature['properties'][field] == miss_val:
                feature['properties'][field] = \
                    feature['properties'][replace_dict[field]]
    
    return(geojson)


def change_features(geojson, change_dict):
    """For all features, if the filed defined in dict key rdoes not match the 
    first regular expression in the value list, but the second, then change 
    the field entry according to  the expressions defined in the third value 
    list entry.
    
    """
    for field in list(change_dict.keys()):
        for feature in geojson['features']:
            for change in change_dict[field]:
                if   (  re.match(
                        change[0],feature['properties'][field]
                        ) == None
                    and re.match(
                        change[1],feature['properties'][field]
                        ) != None
                ):    
                    
                    feature['properties'][field] = re.sub(
                        change[2][0],change[2][1], 
                        feature['properties'][field]
                    )
                
    return(geojson)


def del_features(geojson, del_dict, miss_val = 'missing'):
    """Take a geojson file, and a delete all features for which any field 
    defined in the expression_string_dict keys does not respond to the regular
    expression that is defined in the corresponding regular expression defined
    in the corresponding dictionary value.
    
    """
    remove_list = []
    for field in list(del_dict.keys()):
        for feature in geojson['features']:
            if (re.match(del_dict[field],feature['properties'][field]) == None
                and feature['properties'][field]!=miss_val):
                remove_list.append(feature)
    
    for feature in remove_list:
        geojson['features'].remove(feature)
    
    return(geojson)


def cr_prop_list(geojson, min_num, add_list=[], add_sw=[], 
                 excl_list=[], excl_sw=[]):
    """Return an in_list of all properties, that should be included in the 
    geojson for each feature, and an out_list with all property names that 
    should not be included.
    Take property names from the geojson style input data. Keep all property 
    names that are assigned at least min_num times. In addition keep properties
    in add_list and those that start with expressions in add_sw, but remove 
    properties in excl_list and those that start with expressions in excl_sw.
    
    """
    
    all_prop = []
    
    # Create a list of all properties that are ever assigned in the geojson
    for feature in geojson['features']:
        for key in feature['properties'].keys():
            if key not in all_prop:
                all_prop.append(key)
    
    # Create a dictionary with all properties as keys and 0 as values.
    prop_count = dict(zip(all_prop, [0]*len(all_prop)))
    
    # Count how often each key is assigned:
    for feature in geojson['features']:
        for key in feature['properties'].keys():
            prop_count[key]=prop_count[key] + 1
    
    # Add stuff that starts with expression to respective lists.
    for y in add_sw:
        add_list.extend([x for x in all_prop if x.startswith(y)])
    
    for y in excl_sw:
        excl_list.extend([x for x in all_prop if x.startswith(y)])
    
    in_list = [x for x in all_prop if prop_count[x]>min_num]
    in_list = list(set(in_list) - set(excl_list))
    in_list.extend(add_list)
    # Remove duplicates.
    in_list = list(set(in_list))
    
    out_list = list(set(all_prop) - set(in_list))
    
    return(in_list, out_list, prop_count)



def cr_uniform_geojson(geojson, in_list, out_list, miss_val = 'missing'):
    """Return a geojson style dictionary that is based on the input geojson  
    and contains every property in in_list and assigns miss_val if it was not 
    defined before. Delete all properties in out_list from all features. 
    
    """
    
    all_prop = list(set(in_list).union(out_list))
    
    for feature in geojson['features']:
        for prop in all_prop:
            if prop in in_list:
                if prop not in feature['properties'].keys():
                    feature['properties'].update({prop: miss_val})
            else:
                if prop in feature['properties'].keys():
                    feature['properties'].pop(prop)   
        
    return(geojson)












