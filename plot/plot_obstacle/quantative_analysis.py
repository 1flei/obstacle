import numpy as np
import matplotlib.pylab as plt
import json
import random
from plot_obst_json_util import plot_data_obs_pnts_heatmap, plot_data_obs_pnts, getERPAtTime, plotERP, plotERPArea, get_latest_json_filename, \
    readERPArea, polygon_point_mindist, gt_vessel0, gt_vessel1, gt_vessel2



def min_dist2_poly(qs_json, erp_areas, stop_threshold=0):
    qs = np.array([data['obs_pnts'][1:3] for data in qs_json])
    min_dist2= 1e9
    for q in qs:
        for area in erp_areas:
            # print(q, area, polygon_point_mindist(q, area))
            min_dist2 = np.minimum(polygon_point_mindist(q, area), min_dist2)
            if min_dist2 < stop_threshold:
                return min_dist2
    return min_dist2

def calc_precision(datas, dist_threshold, t = '08:00:01'):
    if datas is None:
        return 0.
    erp_areas = readERPArea(t = t)
    fits_cnt = 0
    for i, data in enumerate(datas):
        # for obs_json in data:
        #     plot_data(obs_json, is_plot_density=False, is_plot_dist=False, is_plot_q=False)
            # plot_data_and_distfig_sgtaxi(obs_json, is_plot_dist=False)
        dist2 = min_dist2_poly(data, erp_areas, stop_threshold= dist_threshold)
        # print(dist2, fits_cnt)
        if dist2 < dist_threshold:
            fits_cnt+= 1

    return 1.* fits_cnt / len(datas)

def calc_recall(datas, dist_threshold, t='08:00:01'):
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

def calc_precision_vessel(datas, dist_threshold, gt=[gt_vessel0]):
    if datas is None:
        return 0.
    fits_cnt = 0
    for i, data in enumerate(datas):
        # plot_data_obs_pnts(data, is_plot_density=False)
        # for obs_json in data:
        #     plot_data(obs_json, is_plot_density=False, is_plot_dist=False, is_plot_q=False)
            # plot_data_and_distfig_sgtaxi(obs_json, is_plot_dist=False)
        dist2 = min_dist2_poly(data, gt)
        # print(dist2, fits_cnt)
        if dist2 < dist_threshold:
            fits_cnt+= 1

    return 1.* fits_cnt / len(datas)

def calc_recall_vessel(datas, dist_threshold, gt=[gt_vessel0]):
    if datas is None:
        return 0.
    fits_cnt = 0
    for erp in gt:
        min_dist2= 1e9
        for data in datas:
            for obs_pnts in data:
                # print(obs_pnts, erp)
                q = np.array(obs_pnts['obs_pnts'][1:3] )
                min_dist2 = np.minimum(polygon_point_mindist(q, erp), min_dist2)
        if min_dist2 < dist_threshold:
            fits_cnt+= 1
    return fits_cnt*1. / len(gt)

def calc_f1score(precision, recall):
    return 2* (precision*recall)/(precision+recall+1e-9)


#-----
def test_sgtaxi(
    name = 'sgtaxi_morning', 
    taus = ['50', '100', '150'], 
    deltas = ['1.5', '2', '2.5'], 
    sigmas = ['0.01', '0.02', '0.03'], 
    erp_time = '08:00:01', 
    dist_threshold = 1e-4
    ):

    best_score = 0
    for sigma in sigmas:
        for tau in taus:
            for delta in deltas:
                filename = get_latest_json_filename(name, sigma=sigma, delta=delta, tau=tau)
                # print(filename)
                with open(filename, 'r') as f:
                    datas_json = json.load(f)
                    datas = datas_json['obst']
                    query_time = datas_json['running_time']
                precision_morning = calc_precision(datas, dist_threshold=dist_threshold, t=erp_time)
                recall_morning = calc_recall(datas, dist_threshold=dist_threshold, t=erp_time)
                f1score = calc_f1score(precision_morning, recall_morning)

                if f1score > best_score:
                    best_tau = tau
                    best_delta = delta
                    best_sigma = sigma
                    best_score = f1score
                    best_precision = precision_morning
                    best_recall = recall_morning
                    best_query_time = query_time

                # if 'taxi_morning' not in best_setting or f1score > best_setting['taxi_morning'][2]:
                #     best_setting['taxi_morning'] = (precision_morning, recall_morning, f1score, query_time, tau, delta)
                print("delta={}, tau={}, sigma={}".format(delta, sigma, tau))
                print(precision_morning, recall_morning, f1score)
    print('{}: best_score={}, best_precision={}, best_recall={}, best_tau={}, best_delta={}, best_sigma={}, best_query_time={}'.format(
        name, best_score, best_precision, best_recall, best_tau, best_delta, best_sigma, best_query_time))
    return best_score, best_tau, best_delta, best_sigma


#-----
def test_vessel(
    name = 'vessel0', 
    taus = ['100', '200', '300'], 
    deltas = ['1.5', '2', '2.5'], 
    sigmas = ['0.05', '0.1', '0.15'], 
    dist_threshold = 0.02, 
    gt = [gt_vessel0]
    ):

    best_score = -1
    for sigma in sigmas:
        for tau in taus:
            for delta in deltas:
                # print("delta={}, tau={}, sigma={}".format(delta, tau, sigma))
                filename = get_latest_json_filename(name, sigma=sigma, delta=delta, tau=tau)
                print(filename)
                with open(filename, 'r') as f:
                    datas_json = json.load(f)
                    datas = datas_json['obst']
                    query_time = datas_json['running_time']
                precision = calc_precision_vessel(datas, dist_threshold=dist_threshold, gt=gt)
                recall = calc_recall_vessel(datas, dist_threshold=dist_threshold, gt=gt)
                f1score = calc_f1score(precision, recall)

                if f1score > best_score:
                    best_tau = tau
                    best_delta = delta
                    best_sigma = sigma
                    best_score = f1score
                    best_precision = precision
                    best_recall = recall
                    best_query_time = query_time

                # if 'taxi_morning' not in best_setting or f1score > best_setting['taxi_morning'][2]:
                #     best_setting['taxi_morning'] = (precision_morning, recall_morning, f1score, query_time, tau, delta)
                print(precision, recall, f1score)
    print('{}: best_score={}, best_precision={}, best_recall={}, best_tau={}, best_delta={}, best_sigma={}, best_query_time={}'.format(
            name, best_score, best_precision, best_recall, best_tau, best_delta, best_sigma, best_query_time))
    return best_score, best_tau, best_delta, best_sigma


if __name__ == '__main__':
    deltas = ['1.5', '2', '2.5']

    taus = ['50', '100', '150']
    sigmas = ['0.01', '0.02', '0.03']

    taus_vessel = [500, 1000, 1500]
    sigmas_vessel = ['0.05', '0.1', '0.15']

    # test_sgtaxi(name = 'sgtaxi_morning', taus=taus, deltas=deltas, sigmas=sigmas)
    # test_sgtaxi(name = 'sgtaxi_afternoon', taus=taus, deltas=deltas, sigmas=sigmas, erp_time='18:00:01')
    # test_vessel(name = 'vessel1', taus=[500, 1000, 1500], deltas=[1.5, 2], sigmas=[0.05], gt=[gt_vessel1, gt_vessel2])
    # test_vessel(name = 'vessel1', taus=taus_vessel, deltas=deltas, sigmas=sigmas_vessel, gt=[gt_vessel1, gt_vessel2])
    test_vessel(name = 'vessel0', taus=taus_vessel, deltas=deltas, sigmas=sigmas_vessel, gt=[gt_vessel0])
    # test_sgtaxi(name = 'vessel1', taus=taus, deltas=deltas, sigmas=sigmas)