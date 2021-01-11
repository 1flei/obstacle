import numpy as np
import matplotlib.pylab as plt
import json
import random
from plot_obst_json_util import plot_data, gt_vessel0, gt_vessel1, gt_vessel2, plot_chosen_area, plot_data_obs_pnts, plot_data_dist, plot_data_and_distfig, plotERP, get_latest_json_filename



filename = get_latest_json_filename(resname= 'vessel1', sigma=0.05, delta=1.5, tau=500, folder='res')

print(filename)
# filename = 'res/vessel1_[10-10_22_31_sigma=0.1_delta=2_tau=300].json'
with open(filename, 'r') as f:
    datas = json.load(f)['obst']


for i, data in enumerate(datas):
    plot_data_obs_pnts(data, is_plot_density=False, is_plot_succ=True, color='b')
    
# plot_chosen_area(gt_vessel0, color='m')
plot_chosen_area(gt_vessel1, color='r')
plot_chosen_area(gt_vessel2, color='r')
plt.show()




