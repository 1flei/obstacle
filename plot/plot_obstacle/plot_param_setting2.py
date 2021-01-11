import numpy as np
import matplotlib.pylab as plt
import json
import random
from plot_obst_json_util import gt_vessel0, gt_vessel1, gt_vessel2, plot_data, plot_chosen_area, plot_data_obs_pnts, \
    plot_data_dist, plot_data_and_distfig_sgtaxi, plotERP, plotERPArea, polygon_point_mindist, readERPArea, get_latest_json_filename
from quantative_analysis import *


plt.rcParams.update({'font.size': 16})

def plot_single_figure(xs, ys, xlabel, ylabel):

    deltas_number = [float(x) for x in deltas]
    for f1score, marker, color, label in zip(f1scores.T, markers, colors, labels):
        plt.plot(deltas_number, f1score, color=color, marker=marker, label=label, markerfacecolor='none', 
        markersize=10)
    
    plt.xlim(0.5, 4.)
    plt.ylim(0, 1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(ncol=2, fontsize=15, loc='upper center')
    # plt.title('F1-score over various delta')
    plt.savefig('f1score_delta.pdf')
    plt.show()

    for querytime, marker, color, label in zip(querytimes.T, markers, colors, labels):
        plt.semilogy(deltas_number, querytime, color=color, marker=marker, label=label, markerfacecolor='none', 
        markersize=10)
    
    plt.xlim(0.5, 4.)
    plt.ylim(1, 3000)
    plt.legend(ncol=2, fontsize=15, loc='upper center')
    plt.xlabel('$\delta$')
    plt.ylabel('Query Time (s)')
    # plt.title('F1-score over various delta')
    plt.savefig('querytime_delta.pdf')
    plt.show()

taus = [1.5, 1.75, 2, 2.25, 2.5]
sigmas = [0.01, 0.015, 0.02, 0.025, 0.03]
deltas = [50, 75, 100, 125, 150]
markers = ['o','s','^','p','*']
colors = ['r','g','b','c','m','y']
labels = ['Morning ERP', 'Afternoon ERP']
# taus = [1.5, 2, 2.5]
# sigmas = [0.01, 0.02, 0.03]
# deltas = [50, 100, 150]

def parameter_setting_delta(sigma=0.01, tau=2, deltas=deltas):
    settings = ['sgtaxi_morning', 'sgtaxi_afternoon']
    erp_times = ['08:00:01', '18:00:01']
    dist_threshold = 1e-4

    xyzs = []
    for setting, erp_time in zip(settings, erp_times):

        xyz = []
        for delta in deltas:
            filename = get_latest_json_filename(setting, sigma=sigma, delta=tau, tau=delta)
            print(filename)
            with open(filename, 'r') as f:
                datas_json = json.load(f)
                datas = datas_json['obst']
                query_time = datas_json['running_time']
            precision_morning = calc_precision(datas, dist_threshold=dist_threshold, t=erp_time)
            recall_morning = calc_recall(datas, dist_threshold=dist_threshold, t=erp_time)
            f1score = calc_f1score(precision_morning, recall_morning)

            print(delta, f1score, query_time)
            xyz += [[delta, f1score, query_time]]
        xyzs += [xyz]
    xyzs = np.array(xyzs)
    print(xyzs)

    #plot x vs. F1-score
    fig = plt.figure(figsize=(5, 3.5))
    for i, xyz, color, marker, label in zip(range(10), xyzs, colors, markers, labels):
        plt.plot(xyz[:, 0], xyz[:, 1], color=color, marker=marker, label=label, markerfacecolor='none', 
            markersize=10)
    
    plt.xlim(np.min(deltas), np.max(deltas))
    plt.ylim(0, 0.8)
    plt.xticks(deltas)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8])
    plt.xlabel('$\delta$')
    plt.ylabel('F1-score')
    plt.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.1)
    plt.legend()
    plt.savefig('param_delta_f1score.pdf')
    plt.savefig('param_delta_f1score.png')
    # plt.show()

    #plot x vs. query-time
    fig = plt.figure(figsize=(5, 3.5))
    for i, xyz, color, marker, label in zip(range(10), xyzs, colors, markers, labels):
        plt.plot(xyz[:, 0], xyz[:, 2], color=color, marker=marker, label=label, markerfacecolor='none', 
            markersize=10)
    
    plt.xlim(np.min(deltas), np.max(deltas))
    plt.xticks(deltas)
    plt.ylim(0, 80)
    plt.yticks([0, 40, 80])
    plt.xlabel('$\delta$')
    plt.ylabel('Query Time (s)')
    plt.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.1)
    plt.legend()
    plt.savefig('param_delta_qt.pdf')
    plt.savefig('param_delta_qt.png')
    # plt.show()

def parameter_setting_tau(sigma=0.01, tau=taus, delta=50):
    settings = ['sgtaxi_morning', 'sgtaxi_afternoon']
    erp_times = ['08:00:01', '18:00:01']
    dist_threshold = 1e-4

    xyzs = []
    for setting, erp_time in zip(settings, erp_times):

        xyz = []
        for tau in taus:
            filename = get_latest_json_filename(setting, sigma=sigma, delta=tau, tau=delta)
            print(filename)
            with open(filename, 'r') as f:
                datas_json = json.load(f)
                datas = datas_json['obst']
                query_time = datas_json['running_time']
            precision_morning = calc_precision(datas, dist_threshold=dist_threshold, t=erp_time)
            recall_morning = calc_recall(datas, dist_threshold=dist_threshold, t=erp_time)
            f1score = calc_f1score(precision_morning, recall_morning)

            xyz += [[tau, f1score, query_time]]
        xyzs += [xyz]
    xyzs = np.array(xyzs)

    #plot x vs. F1-score
    fig = plt.figure(figsize=(5, 3.5))
    for i, xyz, color, marker, label in zip(range(10), xyzs, colors, markers, labels):
        plt.plot(xyz[:, 0], xyz[:, 1], color=color, marker=marker, label=label, markerfacecolor='none', 
            markersize=10)
    
    plt.xlim(np.min(taus), np.max(taus))
    plt.ylim(0, 0.8)
    plt.xticks(taus)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8])
    plt.xlabel(r'$\tau$')
    plt.ylabel('F1-score')
    plt.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.1)
    plt.legend()
    plt.savefig('param_tau_f1score.pdf')
    plt.savefig('param_tau_f1score.png')
    # plt.show()

    #plot x vs. query-time
    fig = plt.figure(figsize=(5, 3.5))
    for i, xyz, color, marker, label in zip(range(10), xyzs, colors, markers, labels):
        plt.plot(xyz[:, 0], xyz[:, 2], color=color, marker=marker, label=label, markerfacecolor='none', 
            markersize=10)
    
    plt.xlim(np.min(taus), np.max(taus))
    plt.xticks(taus)
    plt.ylim(0, 80)
    plt.yticks([0, 40, 80])
    plt.xlabel(r'$\tau$')
    plt.ylabel('Query Time (s)')
    plt.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.1)
    plt.legend()
    plt.savefig('param_tau_qt.pdf')
    plt.savefig('param_tau_qt.png')
    # plt.show()

def parameter_setting_sigma(sigma=sigmas, tau=2, delta=50):
    settings = ['sgtaxi_morning', 'sgtaxi_afternoon']
    erp_times = ['08:00:01', '18:00:01']
    dist_threshold = 1e-4

    xyzs = []
    for setting, erp_time in zip(settings, erp_times):

        xyz = []
        for sigma in sigmas:
            filename = get_latest_json_filename(setting, sigma=sigma, delta=tau, tau=delta)
            print(filename)
            with open(filename, 'r') as f:
                datas_json = json.load(f)
                datas = datas_json['obst']
                query_time = datas_json['running_time']
            precision_morning = calc_precision(datas, dist_threshold=dist_threshold, t=erp_time)
            recall_morning = calc_recall(datas, dist_threshold=dist_threshold, t=erp_time)
            f1score = calc_f1score(precision_morning, recall_morning)

            xyz += [[sigma, f1score, query_time]]
        xyzs += [xyz]
    xyzs = np.array(xyzs)

    #plot x vs. F1-score
    fig = plt.figure(figsize=(5, 3.5))
    for i, xyz, color, marker, label in zip(range(10), xyzs, colors, markers, labels):
        plt.plot(xyz[:, 0], xyz[:, 1], color=color, marker=marker, label=label, markerfacecolor='none', 
            markersize=10)
    
    plt.xlim(np.min(sigmas), np.max(sigmas))
    plt.ylim(0, 0.8)
    plt.xticks(sigmas)
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8])
    plt.xlabel(r'$\sigma$')
    plt.ylabel('F1-score')
    plt.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.1)
    plt.legend()
    plt.savefig('param_sigma_f1score.pdf')
    plt.savefig('param_sigma_f1score.png')
    # plt.show()

    #plot x vs. query-time
    fig = plt.figure(figsize=(5, 3.5))
    for i, xyz, color, marker, label in zip(range(10), xyzs, colors, markers, labels):
        plt.plot(xyz[:, 0], xyz[:, 2], color=color, marker=marker, label=label, markerfacecolor='none', 
            markersize=10)
    
    plt.xlim(np.min(sigmas), np.max(sigmas))
    plt.xticks(sigmas)
    plt.ylim(0, 150)
    plt.yticks([0, 50, 100, 150])
    plt.xlabel(r'$\sigma$')
    plt.ylabel('Query Time (s)')
    plt.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.1)
    plt.legend()
    plt.savefig('param_sigma_qt.pdf')
    plt.savefig('param_sigma_qt.png')
    # plt.show()


if __name__ == '__main__':
    parameter_setting_delta()
    parameter_setting_tau() 
    parameter_setting_sigma() 
    plt.show()