"""
"""

import arcpy
import pandas as pd

from src.process_geodata.functions.f13_zeb_functions import compass_distance
from src.process_numdata.functions.f02_collapse_functions import long_to_wide


def sort_select_and_reshape(df):
    
    df.sort_values(['ID','FMEAS'], inplace=True)
    
    df['id'] = df['ID']
    df['num'] = df.sort_values(['FMEAS']).groupby(['id']).cumcount()
    df.set_index(['num','id'], inplace=True)
    
    sharp_turn_value = abs(df['angle_diff']).quantile(0.99)
    
    df_wide = long_to_wide(df, stubnames=['angle_diff'])
    
    angle_cols = [x for x in df_wide.columns if x.startswith('angle_diff')]
    
    return(df_wide, angle_cols, sharp_turn_value)


def get_and_join_vertices(pp_frame_split, pp_join_pts):
    
    # Make Start and End vertices.
    arcpy.management.FeatureVerticesToPoints(pp_frame_split, 
                                             r"in_memory\start_pts", "START")
    arcpy.management.FeatureVerticesToPoints(pp_frame_split, 
                                             r"in_memory\end_pts", "END")
    
    for fc in [r"in_memory\start_pts", r"in_memory\end_pts"]:
        arcpy.management.AddField(fc, "pt_type", "TEXT")
        with arcpy.da.UpdateCursor(fc, ["pt_type"]) as cursor:
            for r in cursor:
                r[0] = fc[10:]
                cursor.updateRow(r)
    
    # Join Start and endpoints together. 
    arcpy.analysis.SpatialJoin(r"in_memory\start_pts", r"in_memory\end_pts", 
                               pp_join_pts, "JOIN_ONE_TO_ONE", "KEEP_ALL",
                               'fid_ "fid_" true true false 10 Long 0 10,First,#,st_pts,fid_,-1,-1;ref "ref" true true false 10 Text 0 0,First,#,st_pts,ref,0,10;one_direct "one_direct" true true false 2 Text 0 0,First,#,st_pts,one_direct,0,2;compass_a "compass_a" true true false 19 Double 0 0,First,#,st_pts,compass_a,-1,-1;ORIG_FID "ORIG_FID" true true false 10 Long 0 10,First,#,st_pts,ORIG_FID,-1,-1;pt_type "pt_type" true true false 254 Text 0 0,First,#,st_pts,pt_type,0,254;fid1 "fid_" true true false 10 Long 0 10,First,#,end_pts,fid_,-1,-1;ref_1 "ref" true true false 10 Text 0 0,First,#,end_pts,ref,0,10;one_direct_1 "one_direct" true true false 2 Text 0 0,First,#,end_pts,one_direct,0,2;compass_a_1 "compass_a" true true false 19 Double 0 0,First,#,end_pts,compass_a,-1,-1;ORIG_FID_1 "ORIG_FID" true true false 10 Long 0 10,First,#,end_pts,ORIG_FID,-1,-1;pt_type_1 "pt_type" true true false 254 Text 0 0,First,#,end_pts,pt_type,0,254', 
                               "CLOSEST", "0.1 Meters", '')
    
    sdf_joinpts = pd.DataFrame.spatial.from_featureclass(pp_join_pts)
    # Subtract endpoints from startpoint angle to get how much angle has changed 
    # from first to second segment. Negative values are left turn, positive right. 
    sdf_joinpts['angle_diff'] = sdf_joinpts[['compass_a', 'compass_a_1']].apply(
                compass_distance, args=(False,), axis=1
            )
    # Keep only relevant columns. 
    sdf_joinpts = sdf_joinpts[["OID", "ref", "one_direct", "angle_diff", "SHAPE"]]
    sdf_joinpts.spatial.to_featureclass(pp_join_pts)
    
    print("Points with ldm difference saved in  in: " , pp_join_pts)

    return()


def segment_frame_and_get_ldm(pp_frame, pp_frame_split, split_length): 
    
    arcpy.management.Dissolve(pp_frame, r"in_memory\frame_diss", "one_direct;ref", 
                              None, "SINGLE_PART", "UNSPLIT_LINES")
    
    arcpy.management.GeneratePointsAlongLines(r"in_memory\frame_diss", 
                                              r"in_memory\pts", "DISTANCE", 
                                              split_length, None, None)
    arcpy.management.SplitLineAtPoint(r"in_memory\frame_diss", r"in_memory\pts", 
                                      pp_frame_split, "0.01 Meters")
    # Assign new permanent id. 
    arcpy.management.AddField(pp_frame_split, "newid", "LONG")
    
    with arcpy.da.UpdateCursor(pp_frame_split,["newid"]) as cursor:
        for i, r in enumerate(cursor):
            r[0] = i
            cursor.updateRow(r)
    
    # Get linear directional mean
    arcpy.stats.DirectionalMean(pp_frame_split, r"in_memory\dir_mean",
                                "DIRECTION", "newid")
    # Join linear directional mean to to frame_split. 
    sdf_dm = pd.DataFrame.spatial.from_featureclass(r"in_memory\dir_mean").set_index('newid')
    sdf_fsplt = pd.DataFrame.spatial.from_featureclass(pp_frame_split).set_index('newid')
    sdf_fsplt = sdf_fsplt.join(sdf_dm[['CompassA']], how='inner')
    sdf_fsplt.spatial.to_featureclass(pp_frame_split)

    print("Segmented frame with linear directional mean in: ", pp_frame_split)

    return()



