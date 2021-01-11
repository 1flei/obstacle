import numpy as np
import matplotlib.pylab as plt
import json
import random
from plot_obst_json_util import plot_data_obs_pnts_heatmap, plot_data_obs_pnts, getERPAtTime, plotERP, plotERPArea, get_latest_json_filename

plt.rcParams.update({'font.size': 16})

def plot_morning():

    # filename = 'taxi_morning_obst.json'
    # filename = 'taxi_afternoon_obst.json'
    # filename = 'res/taxi_morning_obst_delta2_tau2.326.json'
    # filename = 'res/taxi_afternoon_obst.json'
    filename = 'res/sgtaxi_morning_[10-10_12_58_sigma0.03_delta1.5_tau100].json'
    # filename = 'obst_test_subt_density_sgtaxi.json'


    filename = get_latest_json_filename(resname='sgtaxi_morning', sigma=0.02, delta=2, tau=50, folder='res')
    with open(filename, 'r') as f:
        datas = json.load(f)['obst']


    erp_morning =  getERPAtTime(t='08:00:01')
    erp_afternoon =  getERPAtTime(t='18:00:01')
    print('erp_morning=', len(erp_morning), '   erp_afternoon=', len(erp_afternoon))

    plotERP(cm='m*', t='08:00:01')

    plotERPArea(cm='m-', t='08:00:01')
    # plotERP(cm='m*', t='08:00:01')
    # plot_data_obs_pnts_heatmap(datas, sigma=0.0005)
    for i, data in enumerate(datas):
        plot_data_obs_pnts(data, is_plot_density=False, is_plot_succ=True)
        # for obs_json in data:
        #     plot_data(obs_json, is_plot_density=False, is_plot_dist=False, is_plot_q=False, qalpha=0.5)
            # plot_data_and_distfig_sgtaxi(obs_json, is_plot_dist=False)

    plt.show()

def plot_afternoon():
    filename = get_latest_json_filename(resname='sgtaxi_afternoon', sigma=0.01, delta=1.5, tau=50, folder='res')
    with open(filename, 'r') as f:
        datas = json.load(f)['obst']


    erp_morning =  getERPAtTime(t='08:00:01')
    erp_afternoon =  getERPAtTime(t='18:00:01')
    print('erp_morning=', len(erp_morning), '   erp_afternoon=', len(erp_afternoon))

    plotERP(cm='m*', t='18:00:01')

    plotERPArea(cm='m-', t='18:00:01')
    # plotERP(cm='m*', t='08:00:01')
    # plot_data_obs_pnts_heatmap(datas, sigma=0.0005)
    for i, data in enumerate(datas):
        plot_data_obs_pnts(data, is_plot_density=False, is_plot_succ=True)
        # for obs_json in data:
        #     plot_data(obs_json, is_plot_density=False, is_plot_dist=False, is_plot_q=False, qalpha=0.5)
            # plot_data_and_distfig_sgtaxi(obs_json, is_plot_dist=False)

    plt.show()


plot_morning()

plot_afternoon()