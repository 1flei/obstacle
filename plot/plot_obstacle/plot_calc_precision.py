import numpy as np
import matplotlib.pylab as plt
import json
import random
from plot_obst_json_util import gt_vessel0, gt_vessel1, gt_vessel2, plot_data, plot_chosen_area, plot_data_obs_pnts, \
    plot_data_dist, plot_data_and_distfig_sgtaxi, plotERP, plotERPArea, polygon_point_mindist, readERPArea



plt.rcParams.update({'font.size': 16})

dist_threshold = 1e-4

erp_areas = readERPArea()

def min_dist2_erp(qs_json, erp_areas=erp_areas):
    qs = np.array([data['obs_pnts'][1:3] for data in qs_json])
    min_dist2= 1e9
    for q in qs:
        for area in erp_areas:
            # print(polygon_point_mindist(q, area))
            min_dist2 = np.minimum(polygon_point_mindist(q, area), min_dist2)
    return min_dist2

def calc_precision(datas, dist_threshold=dist_threshold, t = '08:00:01'):
    if datas is None:
        return 0.
    erp_areas = readERPArea(t = t)
    fits_cnt = 0
    for i, data in enumerate(datas):
        # plot_data_obs_pnts(data, is_plot_density=False)
        # for obs_json in data:
        #     plot_data(obs_json, is_plot_density=False, is_plot_dist=False, is_plot_q=False)
            # plot_data_and_distfig_sgtaxi(obs_json, is_plot_dist=False)
        dist2 = min_dist2_erp(data, erp_areas)
        # print(dist2, fits_cnt)
        if dist2 < dist_threshold:
            fits_cnt+= 1
    # plotERP(cm='m*', t=t)
    # plotERPArea(cm='m-', t=t)
    # plt.show()

    return 1.* fits_cnt / len(datas)

def calc_recall(datas, dist_threshold=dist_threshold, t='08:00:01'):
    if datas is None:
        return 0.
    
    erp_areas = readERPArea(t = t)
    fits_cnt = 0
    for erp in erp_areas:
        min_dist2= 1e9
        for data in datas:
            for obs_pnts in data:
                # print(obs_pnts, erp)
                q = np.array(obs_pnts['obs_pnts'][1:3] )
                min_dist2 = np.minimum(polygon_point_mindist(q, erp), min_dist2)
        if min_dist2 < dist_threshold:
            fits_cnt+= 1
    return fits_cnt*1. / len(erp_areas)


def min_dist2_vessel(qs_json, gt=gt_vessel0):
    qs = np.array([data['obs_pnts'][1:3] for data in qs_json])
    min_dist2= 1e9
    for q in qs:
            # print(polygon_point_mindist(q, area))
        min_dist2 = np.minimum(polygon_point_mindist(q, gt), min_dist2)
    return min_dist2

def calc_precision_vessel(datas, dist_threshold=dist_threshold, gt=gt_vessel0):
    if datas is None:
        return 0.
    gt_ = np.zeros_like(gt)
    gt_[:, 0] = gt[:, 1]
    gt_[:, 1] = gt[:, 0]
    fits_cnt = 0
    for i, data in enumerate(datas):
        # plot_data_obs_pnts(data, is_plot_density=False)
        # for obs_json in data:
        #     plot_data(obs_json, is_plot_density=False, is_plot_dist=False, is_plot_q=False)
            # plot_data_and_distfig_sgtaxi(obs_json, is_plot_dist=False)
        dist2 = min_dist2_vessel(data, gt_)
        # print(dist2, fits_cnt)
        if dist2 < dist_threshold:
            fits_cnt+= 1

    # plot_chosen_area(gt)
    # plt.show()  

    return 1.* fits_cnt / len(datas)

def parameter_setting_delta():
    dist_threshold= 1e-4
    filenames = ['taxi_afternoon_obst_delta%s.json'%(delta) for delta in ['0.5', '1', '1.5', '2', '2.5', '3', '3.5', '4']]
    print(filenames)
    for filename in filenames:
        with open(filename, 'r') as f:
            datas = json.load(f)['obst']
        precision_afternoon = calc_precision(datas, dist_threshold=dist_threshold, t='18:00:01')
        print(filename, precision_afternoon)

    filenames = ['taxi_morning_obst_delta%s.json'%(delta) for delta in ['0.5', '1', '1.5', '2', '2.5', '3', '3.5', '4']]
    print(filenames)
    for filename in filenames:
        with open(filename, 'r') as f:
            datas = json.load(f)['obst']
        precision_afternoon = calc_precision(datas, dist_threshold=dist_threshold, t='08:00:01')
        print(filename, precision_afternoon)

    
    dist_threshold_vessel = 0.01
    filenames = ['vessel0_obst_delta%s.json'%(delta) for delta in ['0.5', '1', '1.5', '2', '2.5', '3', '3.5', '4']]
    print(filenames)
    for filename in filenames:
        with open(filename, 'r') as f:
            datas = json.load(f)['obst']
        precision_vessel0 = calc_precision_vessel(datas, dist_threshold=dist_threshold_vessel, gt=gt_vessel0)
        print(filename, precision_vessel0)
        
    filenames = ['vessel1_obst_delta%s.json'%(delta) for delta in ['0.5', '1', '1.5', '2', '2.5', '3', '3.5', '4']]
    print(filenames)
    for filename in filenames:
        with open(filename, 'r') as f:
            datas = json.load(f)['obst']
        precision_vessel1 = calc_precision_vessel(datas, dist_threshold=dist_threshold_vessel, gt=gt_vessel1)
        print(filename, precision_vessel1)
        
    filenames = ['vessel2_obst_delta%s.json'%(delta) for delta in ['0.5', '1', '1.5', '2', '2.5', '3', '3.5', '4']]
    print(filenames)
    for filename in filenames:
        with open(filename, 'r') as f:
            datas = json.load(f)['obst']
        precision_vessel2 = calc_precision_vessel(datas, dist_threshold=dist_threshold_vessel, gt=gt_vessel2)
        print(filename, precision_vessel2)

def calc_f1score(precision, recall):
    return 2* (precision*recall)/(precision+recall+1e-9)

def quantitive_analysis():
    taus = ['1.282', '1.645', '1.96', '2.326', '2.576']
    # taus = ['1.645', '1.96', '2.326', '2.576']
    deltas = ['0.5', '1', '1.5', '2', '2.5', '3', '3.5', '4']

    # suffix = '_naive'
    suffix = ''

    best_setting = {}
    for tau in taus:
        for delta in deltas:
            print('tau=%s, delta=%s'%(tau, delta))

            filename = 'res/taxi_morning_obst_delta%s_tau%s%s.json'%(delta, tau, suffix)
            with open(filename, 'r') as f:
                datas_json = json.load(f)
                datas = datas_json['obst']
                query_time = datas_json['running_time']
            precision_morning = calc_precision(datas, dist_threshold=dist_threshold, t='08:00:01')
            recall_morning = calc_recall(datas, dist_threshold=dist_threshold, t='08:00:01')
            f1score = calc_f1score(precision_morning, recall_morning)

            if 'taxi_morning' not in best_setting or f1score > best_setting['taxi_morning'][2]:
                best_setting['taxi_morning'] = (precision_morning, recall_morning, f1score, query_time, tau, delta)
            print(precision_morning, recall_morning, f1score)

            filename = 'res/taxi_afternoon_obst_delta%s_tau%s%s.json'%(delta, tau, suffix)
            with open(filename, 'r') as f:
                datas_json = json.load(f)
                datas = datas_json['obst']
                query_time = datas_json['running_time']
            precision_afternoon = calc_precision(datas, dist_threshold=dist_threshold, t='18:00:01')
            recall_afternoon = calc_recall(datas, dist_threshold=dist_threshold, t='18:00:01')
            f1score = calc_f1score(precision_afternoon, recall_afternoon)

            if 'taxi_afternoon' not in best_setting or f1score > best_setting['taxi_afternoon'][2]:
                best_setting['taxi_afternoon'] = (precision_afternoon, recall_afternoon, f1score, query_time, tau, delta)
            print(precision_afternoon, recall_afternoon, f1score)

            dist_threshold_vessel = 0.01
            filename = 'res/vessel0_obst_delta%s_tau%s%s.json'%(delta, tau, suffix)
            with open(filename, 'r') as f:
                datas_json = json.load(f)
                datas = datas_json['obst']
                query_time = datas_json['running_time']
            precision_vessel0 = calc_precision_vessel(datas, dist_threshold=dist_threshold_vessel, gt=gt_vessel0)
            recall_vessel0 = 1. if precision_vessel0 >0 else 0.
            f1score = calc_f1score(precision_vessel0, recall_vessel0)

            
            if 'vessel0' not in best_setting or f1score > best_setting['vessel0'][2]:
                best_setting['vessel0'] = (precision_vessel0, recall_vessel0, f1score, query_time, tau, delta)
            print(precision_vessel0, recall_vessel0, f1score)


            filename = 'res/vessel1_obst_delta%s_tau%s%s.json'%(delta, tau, suffix)
            with open(filename, 'r') as f:
                datas_json = json.load(f)
                datas = datas_json['obst']
                query_time = datas_json['running_time']
            precision_vessel1 = calc_precision_vessel(datas, dist_threshold=dist_threshold_vessel, gt=gt_vessel1)
            recall_vessel1 = 1. if precision_vessel1 >0 else 0.
            f1score = calc_f1score(precision_vessel1, recall_vessel1)

            
            if 'vessel1' not in best_setting or f1score > best_setting['vessel1'][2]:
                best_setting['vessel1'] = (precision_vessel1, recall_vessel1, f1score, query_time, tau, delta)
            print(precision_vessel1, recall_vessel1, f1score)


            filename = 'res/vessel2_obst_delta%s_tau%s%s.json'%(delta, tau, suffix)
            with open(filename, 'r') as f:
                datas_json = json.load(f)
                datas = datas_json['obst']
                query_time = datas_json['running_time']
            precision_vessel2 = calc_precision_vessel(datas, dist_threshold=dist_threshold_vessel, gt=gt_vessel2)
            recall_vessel2 = 1. if precision_vessel2 >0 else 0.
            f1score = calc_f1score(precision_vessel2, recall_vessel2)
            # print(precision_vessel0, precision_vessel1, precision_vessel2)


            if 'vessel2' not in best_setting or f1score > best_setting['vessel2'][2]:
                best_setting['vessel2'] = (precision_vessel2, recall_vessel2, f1score, query_time, tau, delta)
            print(precision_vessel2, recall_vessel2, f1score)

    print(best_setting)

if __name__ == '__main__':
    quantitive_analysis()