import numpy as np
import matplotlib.pylab as plt
import json
import random
from plot_obst_json_util import plot_data_obs_pnts_heatmap, plot_data_obs_pnts, getERPAtTime, plotERP, plotERPArea, get_latest_json_filename
from plot_obst_json_util import gt_vessel0, gt_vessel1, gt_vessel2, plot_chosen_area
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

plt.rcParams.update({'font.size': 12})

def plot_vessel():
    fig = plt.figure(figsize=(5, 3.5))

    filename = get_latest_json_filename(resname= 'vessel1', sigma=0.05, delta=1.5, tau=500, folder='res')
    # filename = get_latest_json_filename(resname='sgtaxi_morning', sigma=0.02, delta=2, tau=50, folder='res')
    with open(filename, 'r') as f:
        datas = json.load(f)['obst']


    # erp_morning =  getERPAtTime(t='08:00:01')
    # erp_afternoon =  getERPAtTime(t='18:00:01')
    # print('erp_morning=', len(erp_morning), '   erp_afternoon=', len(erp_afternoon))

    gt = np.array([np.mean(gt_vessel1, axis=0), np.mean(gt_vessel2, axis=0)] )
    plt.plot(gt[:, 0], gt[:, 1], 'r*' , markersize=20, label='Operating Area')
    for i, data in enumerate(datas):
        plot_data_obs_pnts(data, is_plot_density=False, is_plot_succ=True, alpha=0.05, color='b', ls='-', tail_color='royalblue')


    #rectA
    # rectA = Rectangle((103.827, 1.283),0.03,0.021,linewidth=2,edgecolor='k',facecolor='none', zorder=4)
    # plt.gca().add_patch(rectA)
    # plt.text(103.827+0.003, 1.31, 'A', fontsize=16)
    
    rect =  Patch(facecolor='w', edgecolor='b',
                         label='Operating Region')
    star = Line2D([0], [0], marker='*', color='w', label='Scatter',
                        markerfacecolor='r', markersize=20)

    line = Line2D([0], [0], marker='', color='b', ls='-', alpha=0.3)

    plt.legend([line, rect, star], [r'Trajectories in $\mathcal{C}$', 'Obstacle Regions' ,'Operating Area'], 
            ncol=1, fontsize=14)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=1.0)
    # plt.savefig('case_study_vessel_12.pdf')
    plt.savefig('case_study_vessel_12.png')
    plt.show()

def plot_vessel_zoom_in():
    fig = plt.figure(figsize=(5, 3.5))

    filename = get_latest_json_filename(resname= 'vessel1', sigma=0.05, delta=1.5, tau=500, folder='res')
    # filename = get_latest_json_filename(resname='sgtaxi_morning', sigma=0.02, delta=2, tau=50, folder='res')
    with open(filename, 'r') as f:
        datas = json.load(f)['obst']


    # erp_morning =  getERPAtTime(t='08:00:01')
    # erp_afternoon =  getERPAtTime(t='18:00:01')
    # print('erp_morning=', len(erp_morning), '   erp_afternoon=', len(erp_afternoon))

    gt = np.array([np.mean(gt_vessel1, axis=0), np.mean(gt_vessel2, axis=0)] )
    plt.plot(gt[:, 0], gt[:, 1], 'r*' , markersize=20, label='Operating Area')
    for i, data in enumerate(datas):
        plot_data_obs_pnts(data, is_plot_density=False, is_plot_succ=True, alpha=0.05, color='b', ls='-', tail_color='royalblue')


    plt.xlim(103.5, 105)
    plt.ylim(1.1, 1.5)
    
    rect =  Patch(facecolor='w', edgecolor='b',
                         label='Operating Region')
    star = Line2D([0], [0], marker='*', color='w', label='Scatter',
                        markerfacecolor='r', markersize=20)

    line = Line2D([0], [0], marker='', color='b', ls='-', alpha=0.3)

    plt.legend([line, rect, star], [r'Trajectories in $\mathcal{C}$', 'Obstacle Regions' ,'Operating Area'], 
            ncol=1, fontsize=14)

    plt.xticks([103.5, 104, 104.5, 105])
    plt.yticks([1.1, 1.2, 1.3, 1.4, 1.5])

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=1.0)
    # plt.savefig('case_study_vessel_12.pdf')
    plt.savefig('case_study_vessel_12_zoom_in.png')
    plt.show()

plot_vessel_zoom_in()
plot_vessel()