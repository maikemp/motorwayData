"""
"""


import arcpy
import time
import numpy as np

def track_time_start():
    j=0
    print(j)
    start = time.time()
    return(j, start)


def track_time(j, info_freq, start):
    # Put in every loop to  keep track of the time. 
    j+=1
    if j in list(range(0,60000,info_freq)):
        print(j)
        print(round(time.time()-start,1)) 
    return(j)


def euclidean(x,y):
    return(np.sqrt((x[0] - y[0])**2 + (x[1] - y[1])**2))
    

def make_fieldlist(path, excl_fields_list = [], add_fields_dict = {}, 
                   make_type_dict = True, first_fields = [], last_fields=[]):
    """Make a list of all fields in feature class is path, except for the 
    fields in fexcl_fields_list and adding the fieldds in add_fields_dict. If 
    type_dict = True, also return a dictionary with field name and field type.
    
    """    
    
    # Intitate list of field names
    field_names = list(add_fields_dict.keys())
    
    # Extract all field names from path.
    fs = arcpy.ListFields(path)
    field_names = field_names + [f.name for f in fs]
#    for field in arcpy.ListFields(path):
#        field_names.append(field.name)
    # Remove unwanted fieldnames from it.
    fields = list(set(field_names) - set(excl_fields_list))
    
    # Bring first_fields to the front in the given ordering
    for ffs in first_fields:
        fields.pop(fields.index(ffs))
    fields = first_fields + fields
        
    # Bring last_fieldsd to the back of the list
    for lfs in last_fields:
        fields.pop(fields.index(lfs))
        fields.append(lfs)
        
    if make_type_dict:
        # Map wanted types to the given types 
        type_mapping = {
            'Integer':'FLOAT','String':'TEXT', 'SmallInteger':'SHORT', 
            'OID':'OID', 'Geometry':'Geometry', 'Double':'Double', 
            'Single':'Single'
            }
        # Initiate dictionary with names and field type.
        type_dict = {k:v for (k,v) in add_fields_dict.items()}
        
        # Add all wanted fields with name and type to dictionary.
        for field in fs:
            if field.name not in excl_fields_list:
                type_dict.update({field.name: type_mapping[field.type]})
        return(fields, type_dict)
    else:
        return(fields)











