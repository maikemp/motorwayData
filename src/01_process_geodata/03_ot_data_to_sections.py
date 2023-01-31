# -*- coding: utf-8 -*-
"""
Runtime: approx. 5 hours
"""

import os
import arcpy

import src.process_geodata.functions.f03_ot_to_sec_functions as ot_to_sec
import src.process_geodata.functions.f00_general_functions as general
from project_paths import get_project_paths, get_params

arcpy.env.overwriteOutput = True


pp = get_project_paths()
param = get_params(path=os.path.join(pp['lib'], 'parameter.json'))

# Length of the segment.
pdist = param['section_length_m']

# Tolerance for which buffers to keep, dep on length of lines inside of it.
length_div_tolerance = 50
# How often to be informed by the long-running functions.
info_freq = 10000
# Spacial Reference
sr = arcpy.SpatialReference(param['spacial_reference'])

pop_cutpoint_stat = 'ANY'
pop_cutpoint_distance = '3 Meters'
buffer_width = "4 Meters"

dates = param['dates']
# dates = list(set(param['dates']) - set([param["reference_date"]]))

#############################  PATHS OVERVIEW  ################################
# paths:                     .    # in_memory\        .     # Differ by dates #
#-----------------------------------------------------------------------------#
# pp_frame                   .in  #                   .     #    No           #
# pp_cut_points              .in  #                   .     #    No           #
#                            .    # buffer            .temp #    No           #
# pp_ot                      .in  #                   .     #    Yes          #
# pp_change_points           .out #                   .     #    Yes          #
#                            .    # change_points_pop .temp #    Yes          #
# pp_cut_points_pop         ..out # cut_points_pop    .temp #    Yes          #
# pp_change_points_pop_dist ..out #                   .     #    Yes          # 
#-----------------------------------------------------------------------------#

# In paths.
pp_frame = os.path.join(pp['data_out_km'], 'frame.shp')
pp_cut_points = os.path.join(pp['data_out_km'], 'cutpoints.shp')

# the following must later be put in a loop, so everything ist done separately
# for each date.

###########
# date = '2019-07-01'
#
#pp_cut_points_pop = os.path.join(
#        pp['data_out_otsecpt'], 'cut_points_populated'+date+'.shp'
#    )
#pp_change_points_pop_dist = os.path.join(
#        pp['data_out_otsecpt'], 'change_points_dist'+date+'.shp'
#    )
#
#arcpy.CopyFeatures_management(pp_cut_points_pop, r"in_memory\cut_points_pop")
#arcpy.CopyFeatures_management(pp_change_points_pop_dist, r"in_memory\change_points_pop_dist")

###########


for date in dates:
    print("Start with loop for date:", date)
    # In path
    pp_ot = os.path.join(pp['data_out_otshp'], 'motorways_DE_'+date+'Polyline.shp')
    # Temporary out path, due to bug in arcpy.
    pp_change_points = os.path.join(pp['data_temp'], 'change_points'+date+'.shp')
    
    # Out path final for change points.
    pp_change_points_pop_dist = os.path.join(
            pp['data_out_otsecpt'], 'change_points_dist'+date+'.shp'
        )
    # Out path final for cut points. 
    pp_cut_points_pop = os.path.join(pp['data_out_otsecpt'], 
                                     'cut_points_populated'+date+'.shp')
    
    # Save all field names and the respective type in a list and a dictionary.
    keep_fields, type_dict = general.make_fieldlist(
                                 pp_ot, 
                                 ['F_id','ref','FID', 'Shape', 'OBJECTID']
                             )
       
    # Dissolve pp_ot data to only look at relevanatly different segments.  
    old_seg_nr = arcpy.GetCount_management(pp_ot)[0]
    print('Number of segments in input file:', old_seg_nr)
    arcpy.management.Dissolve(
            pp_ot, r"in_memory\ot", dissolve_field=keep_fields, 
            unsplit_lines='UNSPLIT_LINES'
        )
    # Give short information on relevant changes.
    new_seg_nr = arcpy.GetCount_management(r"in_memory\ot")[0]
    print('Number of relevantly differing segments in input:', new_seg_nr)
    
    # Make a point layer with all start points from pp_ot, containing the 
    # attributes, from the polylines in pp_ot and save out to 
    # pp_change_points (temp)
    ot_to_sec.ot_polylines_to_change_points(
            r"in_memory\ot", pp_change_points, keep_fields, type_dict, sr
        )

    
    # Take the before created cutpoints from pp_cut_points, that defined the regular
    # cutpoints for the frame, and populate them with the attributes from the 
    # closest Polyline, and the ID from the frame. 
    # Needs to take a pathbecause of a  bug... save out to pp_cut_points_pop.  
    # Produces some trash, maybe try to stop it by printing somewhere...
    ot_to_sec.populate_cutpoints(
        pp_cut_points, pp_ot, pp_frame, pp_cut_points_pop, 
        keep_fields, pop_cutpoint_distance, pop_cutpoint_stat
    )
    
    # Create a flat buffer around the frame and insert frame attributes like  
    # the ID to the points inside of this buffer. Save results in 
    # r"in_memory\change_points_pop". Keep the buffer in r"in_memory\buffer".
    ### Potential for improvements: widen buffer and if in two buffers, take 
    # the one, for which the point is closer to the respective frame line. 
    ot_to_sec.line_attr_to_close_points(
        pp_frame, pp_change_points,  r"in_memory\change_points_pop", 
        buffer_width, True,  r"in_memory\buffer"
    )
    
    arcpy.CopyFeatures_management(pp_cut_points_pop, r"in_memory\cut_points_pop") 
    
    # In addition to ID, add a numeric ID for all feature classes.
    # Produces a lot of trash, maybe try to stop it by printing somewhere...
#    for featclass in [r"in_memory\change_points_pop",  
#                      r"in_memory\cut_points_pop",    r"in_memory\buffer"]:
#        string_to_numeric_attr(featclass, 'ID', 'ID_n')
        
    ####### ab hier wirds arbeitsintensiv
    # Expected runtime: 1.5 hours on desktop pc...
    # Add length of the repective segment length to changepoints, and length 
    # from segment start to first changepoint to startpoints, and total segment
    # length to startpoints (ideal value = pdist)
    id_list = ot_to_sec.add_segm_length_to_points(
                    pp_buffer = r"in_memory\buffer", 
                    pp_ot = r"in_memory\ot", 
                    pp_changepts = r"in_memory\change_points_pop", 
                    pp_cutpts = r"in_memory\cut_points_pop", 
                    pdist = pdist, 
                    length_div_tolerance = length_div_tolerance, 
                    info_freq = info_freq
                )
    ##########
    # Still some mistakes... starttoend < 470 m oder so häufig durch fehler...
    # Untersuchen. Gelegentlich auch sehr kleine Werte (<10m)

    with open(
        os.path.join(pp['data_temp'], 'id_list' + date + '.txt'), "w"
    ) as f:
        f.write(str(id_list))
    # Save out to final path.
    arcpy.CopyFeatures_management(r"in_memory\change_points_pop", pp_change_points_pop_dist) 
    # Overwrite stuff in pp_cut_points_pop
    arcpy.CopyFeatures_management(r"in_memory\cut_points_pop", pp_cut_points_pop) 
    
    # Cleanup
#    for featclass in [r"in_memory\change_points_pop",  
#                      r"in_memory\cut_points_pop",    r"in_memory\buffer"]:
#        arcpy.Delete_management(featclass)


### Probleme: Immernoch viele Fehler 
        # wenn eine Linie direkt von einer anderen abgeht 
            # Cut point ist gleichzeitig change point für anderen buffer
            # in dem Fall fehlt auch der startpoint für einen eigentlich
            # in takten Buffer, nur weil da noch eine Linie draus abgeht...
        # oder wenn die ot Polyline relativ weit neben der frame line liegt 
            # aber (größtenteils) gerade noch so im buffer ist...
        # manchmal fehlen sogar cut points obwohl alles perfekt aussieht...
    




from datetime import datetime
print(datetime.now().strftime("%H:%M:%S"))


















