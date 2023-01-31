"""
"""
import os
import arcpy
import pandas as pd
import numpy as np


#pp_bisstra_line = pp_sk
#pp_bisstra_nk = pp_nk

def get_locations(pp_bisstra_line, pp_bisstra_nk, pp_frame, pp_temp, 
                  pp_temp_NK, how='junction', dissolve_city_states=True,
                  dissolve_short_A=10000): 
    '''A location is defined as a section of consecutive or parallel segments
    that lie between two motorway intersections or triangles.
    
    * dissolve_short_segments: Either set to False or give a numeric value that
        gives the length of motorways (A..) in Meters that up to which it is 
        considered as one location regardless of borders etc. Give length into
        one direction, will be doubled by function to measure both directions.
    '''
    
    nk_string = r"in_memory\NK__"
    # Get what is needed from the bisstra shapefiles.     
    arcpy.management.MakeFeatureLayer(
            pp_bisstra_line, r"in_memory\A", 
            "Str_Klasse = 'A' And Sk_Subtyp_ = 'Abschnitt'", None
        )
    
    if how == 'junction':
        arcpy.management.MakeFeatureLayer(
            pp_bisstra_nk, nk_string, 
            "NK_Knotenp = 'Autobahndreieck' Or NK_Knotenp = 'Autobahnkreuz'", 
            None
        )
    elif how == 'border':    
        _get_isolated_border_points(pp_bisstra_nk, pp_frame, nk_string,
                                    dissolve_city_states=dissolve_city_states,
                                    dissolve_short_A=dissolve_short_A)    
    
    # Raise error if bug leads to empty class (if so change save location name)
    if int(arcpy.management.GetCount(nk_string)[0])<50:
        raise ValueError("Something went wrong saving location cutpoints.")
        
    arcpy.management.CopyFeatures(nk_string, pp_temp_NK)

    # Dissolve lines by 
    arcpy.management.Dissolve(r"in_memory\A", r"in_memory\A_dissolved", 
                              "Str_Kennun", None, "MULTI_PART", "DISSOLVE_LINES")
    arcpy.analysis.Buffer(r"in_memory\A_dissolved", r"in_memory\A_dissolved_buff", 
                          "50 Meters", "FULL", "FLAT", "NONE", None, "PLANAR")
    arcpy.analysis.Buffer(nk_string, r"in_memory\NK_Buffer", "200 Meters", 
                          "FULL", "ROUND", "NONE", None, "PLANAR")
    arcpy.analysis.Erase(r"in_memory\A_dissolved_buff", r"in_memory\NK_Buffer", 
                         r"in_memory\A_diss_buff_cutout", None)
    arcpy.management.MultipartToSinglepart(r"in_memory\A_diss_buff_cutout", 
                                           r"in_memory\locations")
    # Add a location id. 
    arcpy.AddField_management(r"in_memory\locations", "location", "LONG")
    i = 1
    with arcpy.da.UpdateCursor(r"in_memory\locations", ["location"]) as uc:
        for r in uc:
            r[0] = i
            i+=1
            uc.updateRow(r)
    # Spatially join location_id from locations to r"in_memory\A" to get df  
    # with BISStra Sk_Kennung and unique locations location_id. 
    arcpy.analysis.SpatialJoin(
            r"in_memory\A",r"in_memory\locations", r"in_memory\A_locationJoin", 
            "JOIN_ONE_TO_ONE", "KEEP_ALL", match_option="HAVE_THEIR_CENTER_IN", 
            search_radius="100 Meters"
        ) 
    
    # Make a copy to pp_temp.
    arcpy.management.CopyFeatures(r"in_memory\A_locationJoin", pp_temp)
    # Convert to df. 
    df_bi = pd.DataFrame.spatial.from_featureclass(
            r"in_memory\A_locationJoin"
        )[['Sk_Kennung', 'location']].drop_duplicates()
    
    return(df_bi)


def _get_isolated_border_points(
        pp_bisstra_nk, pp_frame, save_out_string=r"in_memory\NK_",
        dissolve_city_states=True, dissolve_short_A=False):
    
    arcpy.management.MakeFeatureLayer(
        pp_bisstra_nk, r"in_memory\NK", 
        "NK_Knotenp = 'Landesgrenze' And NK_GeoStrK LIKE '%A%'", 
        None
    )     
    # Join information from pp_frame to border points to be able to consider
    # borders on different motorways separately (get ref=A... variable). 
    arcpy.analysis.SpatialJoin(r"in_memory\NK", pp_frame, r"in_memory\NK_j",
                               "JOIN_ONE_TO_ONE", "KEEP_ALL",  
                               match_option="CLOSEST", 
                               search_radius="40 Meters")
    # Buffer border points.
    arcpy.analysis.Buffer(r"in_memory\NK_j", r"in_memory\NK_Buffer", 
                          "15 Kilometers", "FULL", "ROUND", "NONE")
    # Add variable NK_plain to Buffer to identify borders without direction.
    df_temp = pd.DataFrame.spatial.from_featureclass(r"in_memory\NK_Buffer")   
    df_temp['NK_plain'] = df_temp['NK_Name'].apply(lambda x: 
        '/'.join(sorted(x.replace('(','').replace(')','').split('/'))))
    df_temp.spatial.to_featureclass(location=r"in_memory\NK_Buffer")
    # Dissolve border point buffer by border (which two states?) and motorway.
    arcpy.management.Dissolve(r"in_memory\NK_Buffer", r"in_memory\NK_diss",
                              "ref;NK_plain", "Nr_ COUNT", "SINGLE_PART", 
                              "DISSOLVE_LINES")
    df_buff = pd.DataFrame.spatial.from_featureclass(r"in_memory\NK_diss")
    df_NK = pd.DataFrame.spatial.from_featureclass(r"in_memory\NK_j")
    df_NK['NK_plain'] = df_NK['NK_Name'].apply(lambda x: 
        '/'.join(sorted(x.replace('(','').replace(')','').split('/'))))

    # Create a dictionary that for each buffer, contain a df of all points  
    # belonging to that point.
    df_dict = {}
    for i in range(len(df_buff)):
        row_i = df_buff.iloc[i]
        df_temp = pd.DataFrame()
        for j in range(len(df_NK)):
            row_j = df_NK.iloc[j]
            if ((row_i.SHAPE.contains(row_j.SHAPE)) &  
                (row_i.nk_plain==row_j.NK_plain) & (row_i.ref==row_j.ref)):
                df_temp = df_temp.append(row_j)
        df_dict.update({row_i.OID: df_temp})
    df_NK_new_FID = []
    for k in list(df_dict.keys()):
        v = df_dict[k].copy()
        if len(v) > 1:
            if len(v) % 2 == 0:
                # Drop if road crosses same border twice in short distance.
                df_dict.pop(k)
            else:
               # If larger zero and uneven (= actually crossing border in the 
               # longer run): keep arbitrary point.
               df_dict[k] = v[int((len(v)-1)/2):int((len(v)+1)/2)]
               df_NK_new_FID = df_NK_new_FID+[df_dict[k].FID.values[0]]
        else:
               # If only one point in the buffer, just keep it.
               df_NK_new_FID = df_NK_new_FID+[df_dict[k].FID.values[0]]
    # Keep only before selected border points. 
    df_NK_new = df_NK[df_NK['FID'].isin(df_NK_new_FID)]
    
    if dissolve_city_states == True:
        # Dissolve HH into SH,HB into NI,BE into BB by dropping border points.
        df_NK_new = df_NK_new[
                df_NK_new['NK_plain'].isin(['HH/SH', 'HB/NI', 'BB/BE'])==False
            ]
    if dissolve_short_A != False:
        max_len = dissolve_short_A*2
        
        arcpy.management.Dissolve(r"in_memory\A", r"in_memory\A_dissolved", 
                            "Str_Kennun", None, "MULTI_PART", "DISSOLVE_LINES")
        
        # If motorway is very short, drop all NK points on them to keep as one. 
        df_t = pd.DataFrame.spatial.from_featureclass(r"in_memory\A_dissolved")  
        df_t['length'] = df_t['SHAPE'].apply(lambda x: x.length)
        keep_whole = [x[0]+' '+x[1:] for x in df_t[df_t['length'] < max_len].Str_Kennun]
        df_NK_new = df_NK_new[df_NK_new['ref'].isin(keep_whole)==False]
        
    df_NK_new.spatial.to_featureclass(location=save_out_string)  
    
    print('Saved results out to:', save_out_string)
    return()




















