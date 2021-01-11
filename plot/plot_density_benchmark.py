import os
import re
import numpy as np
import matplotlib.pylab as plt
import matplotlib
import json
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

from scipy.spatial import ConvexHull
from itertools import chain, count
from scipy.interpolate import interp1d

from collections import defaultdict
from plot_utility import *

import glob, os

possible_datasets = ['sgtaxi']

dataset_labels_map = {
    'sgtaxi': 'sgtaxi', 
}

possible_methods = ['lsh_density_pie_cm']

datasets = ['sgtaxi']
# datasets = ['vessel']
methods = ['lsh_density_pie_cm', 'lsh_density_pie_ht', 'lsh_density_bks']

method_colors = ['red', 'green', 'blue', 'purple', 'darkgoldenrod', 'c', 'y']
method_markers = ['o', '^', 's', 'd', '*', 'p', 'x']

result_folder_path = 'density_test_res'

def get_latest_res(file_prefix):
    pat = re.compile(r'%s_\[.*\]'%(file_prefix))
    
    latestFile = None
    for file in os.listdir(result_folder_path):
        if pat.match(file):
            if latestFile is None:
                latestFile = file
            elif latestFile < file:
                latestFile = file
    if latestFile is None:
        return None
    return os.path.join(result_folder_path, latestFile)

def get_file_prefix(datasetName, methodName):
    return '%s_%s'%(datasetName, methodName)


#res is a file with multiple jsons
def parse_res(filename):
    with open(filename, 'r') as f:
        data_str = f.read()
        cnt = 0
        bracket_beg = 0
        for curloc, ch in enumerate(data_str):
            if ch=='{':
                if cnt ==0:
                    bracket_beg = curloc
                cnt += 1
            elif ch=='}':
                cnt -= 1
                if cnt == 0:
                    json_str = data_str[bracket_beg:curloc+1]
                    data = json.loads(json_str)
                    yield data
                    
def getratio(res):
    return res[1]   
def getrecall(res):
    return res[3]   
def gettime(res):
    return res[2]
def getindexingtime(res):
    return res[4]
def getindexsize(res):
    return res[5]
def get_l(res):
    return int(res[0][1]['L'])
def get_c(res):
    return int(res[0][0])
def get_time(res):
    return float(res[1][2])
def get_recall(res):
    return float(res[1][3])
def get_p(res):
    return int(res[0][1]['p'])

#fine the best parameter setting at given recall level
def best_config_at_recall_level(xys, settings, ress, recallThreshold):
    besttime = 1e9
    bestsetting = settings[0]
    bestres = ress[0]
    
    for xy, setting, res in zip(xys, settings, ress):
        recall = xy[1]
        time = xy[0]
        
        if recall > recallThreshold and time < besttime:
            besttime = time
            bestsetting = setting
            bestres = res
    return bestsetting, bestres

def config_more_than_recall_level(xys, settings, ress, recallThreshold):
    ret_settings = []
    ret_ress = []
    ret_times = []
    
    for xy, setting, res in zip(xys, settings, ress):
        recall = xy[1]
        time = xy[0]
        
        if recall > recallThreshold:
#             print(recall, time, setting, res)
            ret_times += [time]
            ret_settings += [setting]
            ret_ress += [res]
    return ret_times, ret_settings, ret_ress

def lower_bound_curve(xys, extra_pnts = [[0, 0]]):
#     print(xys, extra_pnts)
#     xys = np.append(xys, extra_pnts, axis=0)
    
    eps = np.random.normal(size=xys.shape) * 1e-4
    xys += eps
#     print(xys)
#     xys = np.array(sorted(xys, key=lambda x:x[1]) )
#     print(xys)
    hull = ConvexHull(xys)
#     print(hull.vertices)
    
    hull_vs = xys[hull.vertices]
    
    v1s = []
    maxv0 = [-1, -1]
    for v0, v1 in zip(hull_vs, chain(hull_vs[1:], hull_vs[:1])):
#         print(v0, v1)
        if v0[1] > v1[1] and v0[0] > v1[0]:
#             plt.semilogy([v0[1], v1[1]], [v0[0], v1[0]],  'k-')
            v1s = np.append(v1s, v1, axis=-1)
            if v0[1] > maxv0[1]:
                maxv0 = v0
                
    
    vs = np.array(np.append(maxv0, v1s)).reshape(-1, 2)
#     print(vs)
#     plt.semilogy(vs[:, 1], vs[:, 0], 'k-')
        
    f = interp1d(vs[:, 1], vs[:, 0])
    
    minx = np.min(vs[:, 1])+1e-6
    maxx = np.max(vs[:, 1])-1e-6
    x = np.arange(minx, maxx, 1)
    y = list(map(f, x))
#     print(x, y)
#     plt.semilogy(x, y, 'k-')
    return x, y

def lower_bound_curve2(xs, ys):
#     print('xs, ys', xs, ys)
    xys = np.zeros(shape=(len(xs), 2))
    xys[:, 0] = xs
    xys[:, 1] = ys
    
    if len(xs)>2 and xs[-1]>0:
        hull = ConvexHull(xys)
        
        hull_vs = xys[hull.vertices]
        ret_vs = []
        
#         print("hull_vs: ", hull_vs)
        
        pflg = False
        for v0, v1, v2 in zip(chain(hull_vs[-1:], hull_vs[:-1]), hull_vs, chain(hull_vs[1:], hull_vs[:1])):
    #         print(v0, v1)
            if v0[0] < v1[0]:
    #             plt.semilogy([v0[1], v1[1]], [v0[0], v1[0]],  'k-')
                ret_vs = np.append(ret_vs, v1, axis=-1)
            elif v1[0] < v2[0]:
                ret_vs = np.append(ret_vs, v1, axis=-1)
        ret_vs = ret_vs.reshape((-1, 2))
        ret_vs = np.array(sorted(ret_vs, key=lambda x:x[0]) )
        return ret_vs
    return xys

    
def plot_indexing_phase(datasets, methods, method_labels, fig_width=0.55+3.333*len(datasets), fig_height=6.2):
    plt_helper = PlotHelper(plt, fig_width, fig_height)
    plt_helper.plot_subplots_adjust()
    n_datasets = len(datasets)
    
    recall_level = 50
#     recall_level = 85
    
#     y_lim_dataset = {
#         ('Mnist784', ):(0.1, 10),
#     }
#     y_ticks_dataset = {
#         'Mnist784':(0.1, 10),
#     }
    
    for di, (dataset, dataset_label) in enumerate(zip(datasets, dataset_labels)):
        ax_size = plt.subplot(2, n_datasets, di+1)
        plt.xlabel('Index size (GB)')
        
        plt.title(dataset_label)
        ax_time = plt.subplot(2, n_datasets, n_datasets+di+1)
#         plt.xlim(0, 100)
        plt.xlabel('Indexing time (s)')
            
#             if dataset in y_lim_dataset:
#                 ax_size.set_ylim(y_lim_dataset[dataset])
#                 ax_time.set_ylim(y_lim_dataset[dataset])
#             if dataset in x_lim_size_dataset:
#                 ax_size.set_xlim(x_lim_size_dataset[dataset])
#             if dataset in x_lim_time_dataset:
#                 ax_time.set_xlim(x_lim_time_dataset[dataset])
                
            
#         plt.ylim(ymin=0)
        
        miny = 1e9
        maxy = -1e9
        
        for method_idx, method, method_label, method_color, method_marker in zip(count(), methods, method_labels, method_colors, method_markers):
            filename_prefix = get_file_prefix(dataset, method, distance)
            filename = get_latest_res(filename_prefix)
            if filename is None:
                continue
            print('-------------', filename, '------------')
            xys = []
            settings = []
            ress = []
            
            index_timesize_dict = defaultdict(list)
            for setting, res in parse_res(filename):
                qtime = gettime(res)
                qrecall = getrecall(res)
                index_time = getindexingtime(res)
                index_size = getindexsize(res)
                
                index_timesize_dict[(index_time, index_size)] += [[qrecall, qtime]]

#             print(xys)

            index_times, index_sizes, qtimes_at_50recall = [], [], []
            for (index_time, index_size), qrecall_times in index_timesize_dict.items():
#                 print(index_size, qrecall_times)
                qrecall_times = np.array(qrecall_times+[[0, 0]])
                qrecalls = qrecall_times[:, 0]
                qtimes = qrecall_times[:, 1]
                
                if np.max(qrecalls) > recall_level:                    
                    
#                     print(qrecalls, qtimes)
                    f = interp1d(qrecalls, qtimes)
                    time_at_50recall = f(recall_level)
                
#                     print('iit', index_time, index_size, time_at_50recall)
                    index_times += [index_time]
                    index_sizes += [index_size]
                    qtimes_at_50recall += [time_at_50recall]
            
            index_times = np.array(index_times)
            index_sizes = np.array(index_sizes)
            qtimes_at_50recall = np.array(qtimes_at_50recall)
            
            
            index_size_qtimes = lower_bound_curve2(index_sizes/1e9, qtimes_at_50recall)
            if len(index_size_qtimes)>0:
                print('min_qtime=', np.min(index_size_qtimes[:, 1]))
                
                ax_size.semilogy(index_size_qtimes[:, 0], index_size_qtimes[:, 1], '-', color=method_color, marker=method_marker, label=method_label if di==0 else "", 
                            markerfacecolor='none', markersize=10)
                
                
                index_time_qtimes = lower_bound_curve2(index_times, qtimes_at_50recall)
                print(method, index_time_qtimes)
                ax_time.semilogy(index_time_qtimes[:, 0], index_time_qtimes[:, 1], '-', color=method_color, marker=method_marker, label="", 
                            markerfacecolor='none', markersize=10, zorder=len(methods)-method_idx)
                
                miny = min(miny, np.min(index_time_qtimes[:, 1]) )
                maxy = max(maxy, np.max(index_time_qtimes[:, 1]) ) 
                miny = min(miny, np.min(index_size_qtimes[:, 1]) )
                maxy = max(maxy, np.max(index_size_qtimes[:, 1]) ) 
                
        plt_helper.set_y_axis_close(ax_time, miny, maxy)
        plt_helper.set_y_axis_close(ax_size, miny, maxy)
        if di==0:
            ax_size.set_ylabel('Query time (ms)')
            ax_time.set_ylabel('Query time (ms)')
                
#             print('xys_=', xys_)
#             if len(xys_) > 1:
#                 query_time, indexing_time = lower_bound_curve(xys_)
#             else:
#                 query_time, indexing_time = xys_[:, 0], xys_[:, 1]
#             
# #             ax.plot(xys[:, 1], xys[:, 0], '.', color=method_color, marker=method_marker, label=method_label)
# #             ax.semilogy(xys[:, 1], xys[:, 0], '.', color=method_color, marker=method_marker, label=method_label)
#             ax.semilogy(indexing_time, query_time, '-', color=method_color, marker=method_marker, label=method_label, markevery=10, 
#                         markerfacecolor='none', markersize=10)
    
#     plt.figlegend(fontsize=16, bbox_to_anchor=(0.07,0.85,0.75,0.2), loc="center",
#                 mode="expand", borderaxespad=0, ncol=len(methods))
    plt_helper.plot_fig_legend(ncol=len(methods))    
    plt_helper.plot_and_save('saved_figure/indexing_time_recall_%s'%distance)
    
    
                    
def plot_methods(datasets, methods, fig_width = 3.*len(datasets), fig_height = 2.7 + 0.8):
    # plt_helper = PlotHelper(plt, fig_width, fig_height)
    # plt_helper.plot_subplots_adjust()
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.set_xlabel('avg_query_time')
    ax.set_ylabel('memory_usage')
    ax.set_zlabel('log10(MSE)')

    # ax.set_xlim(0, 0.001)
    # ax.set_ylim(0, 5e9)

    # ax.zaxis.set_scale('log')

    for di, dataset in enumerate(datasets):
        for mi, (method, method_color, method_marker) in enumerate(zip(methods, method_colors, method_markers)):
            file_prefix = get_file_prefix(dataset, method)
            filename = get_latest_res(file_prefix)
            if filename is None:
                continue
            print(filename)

            xyzs = []
            for res_json in parse_res(filename):
                xyzs += [[res_json['avg_query_time'], res_json['memory_usage'], res_json['mse']]]

                # t = ''
                # if 'n_repeat' in res_json:
                #     t += str(res_json['n_repeat']) + ' '
                # if 'ht_size' in res_json:
                #     t += str(res_json['ht_size']) + ' '
                # if 'n_cm_repeat' in res_json:
                #     t += str(res_json['n_cm_repeat']) + ' '
                # if 'bks_size' in res_json:
                #     t += str(res_json['bks_size']) + ' '
                # ax.text(res_json['avg_query_time'], res_json['memory_usage'], np.log10(res_json['mse']), t)

            xyzs = np.array(xyzs)
            
            # ax.plot_trisurf(xyzs[:, 0], xyzs[:, 1], xyzs[:, 2], cmap='hot')
            ax.scatter(xyzs[:, 0], xyzs[:, 1], np.log10(xyzs[:, 2] ), marker=method_marker, color=method_color, label=method, depthshade=False)

    plt.legend()
    plt.show()
    
if __name__ == '__main__':
    # datasets = ['sgtaxi']
    datasets = ['vessel']
    plot_methods(datasets, methods)
    # plot_methods_recall_ratio(datasets, methods, method_labels)
    # plot_indexing_phase(datasets, methods, method_labels)
#     plot_indexing_time(datasets, methods, method_labels)