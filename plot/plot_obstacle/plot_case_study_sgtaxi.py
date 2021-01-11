import numpy as np
import matplotlib.pylab as plt
import json
import random
from plot_obst_json_util import plot_data_obs_pnts_heatmap, plot_data_obs_pnts, getERPAtTime, plotERP, plotERPArea, get_latest_json_filename
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

plt.rcParams.update({'font.size': 12})

def plot_morning():
    fig = plt.figure(figsize=(7, 4))

    filename = get_latest_json_filename(resname='sgtaxi_morning', sigma=0.02, delta=2, tau=50, folder='res')
    with open(filename, 'r') as f:
        datas = json.load(f)['obst']


    erp_morning =  getERPAtTime(t='08:00:01')
    erp_afternoon =  getERPAtTime(t='18:00:01')
    print('erp_morning=', len(erp_morning), '   erp_afternoon=', len(erp_afternoon))

    plotERP(cm='r*', t='08:00:01')
    # plotERPArea(cm='m-', t='08:00:01')
    # plotERP(cm='m*', t='08:00:01')
    # plot_data_obs_pnts_heatmap(datas, sigma=0.0005)
    for i, data in enumerate(datas):
        plot_data_obs_pnts(data, is_plot_density=False, is_plot_succ=True, alpha=0.05, color='b', ls='-', tail_color='royalblue')
        # for obs_json in data:
        #     plot_data(obs_json, is_plot_density=False, is_plot_dist=False, is_plot_q=False, qalpha=0.5)
            # plot_data_and_distfig_sgtaxi(obs_json, is_plot_dist=False)


    plt.xlim(103.75, 104)
    plt.ylim(1.26, 1.4)

    plt.xticks([103.75, 103.8, 103.85, 103.9, 103.95, 104])
    # plt.yticks([1.26, 1.4])

    #rectA
    rectA = Rectangle((103.827, 1.283),0.03,0.021,linewidth=2,edgecolor='k',facecolor='none', zorder=4)
    plt.gca().add_patch(rectA)
    plt.text(103.82, 1.3, 'A', fontsize=16)

    #rectB
    rectB = Rectangle((103.852, 1.371),0.03,0.021,linewidth=2,edgecolor='k',facecolor='none', zorder=4)
    plt.gca().add_patch(rectB)
    plt.text(103.852+0.012, 1.363, 'B', fontsize=16)

    
    rect =  Patch(facecolor='w', edgecolor='b',
                         label='Operating Region')
    star = Line2D([0], [0], marker='*', color='w', label='Scatter',
                        markerfacecolor='r', markersize=20)

    line = Line2D([0], [0], marker='', color='b', ls='-', alpha=0.3)

    plt.legend([line, rect, star], ['Trajectories in $\mathcal{C}$', 'Obstacle Regions' ,'ERP Gantries'], 
            ncol=1, loc='lower right', fontsize=14)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=1.0)
    plt.savefig('case_study_taxi_morning.pdf')
    plt.savefig('case_study_taxi_morning.png')
    plt.show()


plot_morning()

# plot_afternoon()