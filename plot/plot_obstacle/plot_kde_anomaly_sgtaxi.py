import numpy as np
import matplotlib.pylab as plt
import json
import random
from plot_obst_json_util import plot_data_obs_pnts_heatmap, plot_data_obs_pnts, getERPAtTime, plotERP, plotERPArea
from matplotlib.lines import Line2D

plt.rcParams.update({'font.size': 12})

# filename = 'taxi_morning_obst.json'
# filename = 'taxi_afternoon_obst.json'
# filename = 'res/taxi_morning_obst_delta2_tau2.326.json'
# filename = 'res/taxi_afternoon_obst.json'
filename = 'res/sgtaxi_morning_[10-12_14_56_threshold1.0].json'
# filename = 'obst_test_subt_density_sgtaxi.json'
with open(filename, 'r') as f:
    datas = json.load(f)

    

def plot_kde_anomaly(data_pnt, color='r', label=''):
    obs_pnts = np.array(data['flat_spatial_array'])
    plt.plot(obs_pnts[::2], obs_pnts[1::2], color=color, ls='-', alpha=0.3, label=label)


erp_morning =  getERPAtTime(t='08:00:01')
erp_afternoon =  getERPAtTime(t='18:00:01')
print('erp_morning=', len(erp_morning), '   erp_afternoon=', len(erp_afternoon))

plotERP(cm='c*', t='08:00:01')

# star = Line2D([0], [0], marker='*', color='w', label='Scatter',
#                     markerfacecolor='orange', markersize=20)
# plotERPArea(cm='m-', t='08:00:01')
# plotERP(cm='m*', t='08:00:01')
# plot_data_obs_pnts_heatmap(datas, sigma=0.0005)
for i, data in enumerate(datas):
    plot_kde_anomaly(data, label='Anomaly Trajectory' if i==0 else '')
    # for obs_json in data:
    #     plot_data(obs_json, is_plot_density=False, is_plot_dist=False, is_plot_q=False, qalpha=0.5)
        # plot_data_and_distfig_sgtaxi(obs_json, is_plot_dist=False)

plt.xlim(103.5, 104.3)
plt.ylim(1, 1.8)

plt.xticks([103.5, 103.7, 103.9, 104.1, 104.3])
plt.yticks([1., 1.2, 1.4, 1.6, 1.8])

plt.legend(fontsize=16)
plt.savefig('kde_anomaly.png', bbox_inches='tight', loc='best')
plt.show()




