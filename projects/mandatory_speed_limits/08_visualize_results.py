"""

Runtime: approx. 1 Minute.
"""

import os 
import json
import pandas as pd

from project_paths import get_project_paths
import src.projects.mandatory_speed_limits.functions.visualize_functions as visual
from src.process_numdata.functions.concat_functions import load_csv_and_fix_id


pp = get_project_paths(laptop=False)



p_in = pp['data_out_a_grf']
p_out = pp['fin_data']

s = "" # "", "r_", "r_sl_", "r_ps_"

exps = ["b"] # ["a", "b", "ba"]
limits = ["maxspeed_100",  "maxspeed_120", "maxspeed_130"]
outcomes = ['total_rate', 'fatal_rate', 'severe_rate', 'light_rate']

in_p_me = os.path.join(pp['fin_data'], 'df_max_ext_imputed.csv')
df_me = load_csv_and_fix_id(in_p_me)

in_p = os.path.join(pp['fin_data'], 'analysis_df.csv')
df = load_csv_and_fix_id(in_p)
crash_counts = ['total', 'fatal', 'severely_injured', 'lightly_injured']
df = df.join(df_me[crash_counts])

cates_d = json.load(open(os.path.join(pp['lib'], 'cates_cutoffs.json')))


# Make a dictionary of rates and of cates and derive CI bounds and semi 
# elasticities. 
# Also add info to the df.
cd, df = visual.make_data_dicts(exps, limits, outcomes, df, "cates", 
                                cates_d, p_in, s)
exps_2 = [x for x in exps if not x in ["ab", "ba"]]
ra, df = visual.make_data_dicts(exps_2, limits, outcomes, df, "rate_q", 
                                    cates_d, p_in, s)


# Make graphics for each setup and each CATE.
# Colorful color list:
# color_list = ['grey', 'indianred', 'goldenrod', 'skyblue']
color_list=['dimgray', 'black', 'darkgray', 'silver']
# Plot ATEs.
visual.plot_ates(["ate", "ate_ov"], cd, s, exps, outcomes, limits, pp,
                 color_list)

# Graphics of CATEs
# Makes no sense to look into fatal rate --> too much uncertainty already!
visual.plot_cates({k:v for k, v in cates_d.items() if k in ['AADT','HT_share','ramps']}, 
                  cd, s, exps, outcomes = [x for x in outcomes if x!="fatal_rate"],
                  limits=limits, pp=pp, color_list=color_list)


