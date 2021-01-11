import numpy as np
import matplotlib.pylab as plt
import json
import random
from plot_obst_json_util import plot_data, gt_vessel0, gt_vessel1, gt_vessel2, plot_chosen_area, plot_data_obs_pnts, plot_data_dist, plot_data_and_distfig, plotERP

# filename = 'vessel2.json'
filename = 'res/vessel0_[10-07_14_24_sigma=0.1_delta=2_tau=100].json'
with open(filename, 'r') as f:
    datas = json.load(f)['obst']



# for i, data in enumerate(datas):
#     print(i)
#     if plot_data(data, z_threshold=1.68, is_plot_density=True, is_plot_dist=True):
#         # plot_chosen_area(gt)
#         pass

# plotERP()
# plt.show()

for i, data in enumerate(datas):
    plot_data_obs_pnts(data, is_plot_density=True, is_plot_succ=True, color='b')
    
plot_chosen_area(gt_vessel0, color='m')
plt.show()




