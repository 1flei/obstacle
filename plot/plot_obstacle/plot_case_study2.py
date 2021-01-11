import numpy as np
import matplotlib.pylab as plt
import json
import random
from plot_obst_json_util import plot_data, gt_vessel1, gt_vessel2, plot_chosen_area, plot_data_obs_pnts, plot_data_dist, plot_data_and_distfig_sgtaxi, plotERP
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

plt.rcParams.update({'font.size': 12})

def plot_erp_case_study():

    fig = plt.figure(figsize=(7, 4.5))
    filename = 'res/taxi_morning_obst_delta2_tau1.645.json'
    # filename = 'obst_test_subt_density_sgtaxi.json'
    with open(filename, 'r') as f:
        data_taxi_morning = json.load(f)['obst']

    def plotLegendTaxi():    
        class ObjB(object):
            pass
        class ObjR(object):
            pass
        class ObjG(object):
            pass

        class data_handler(object):
            def __init__(self, color):
                self.color = color
            def legend_artist(self, legend, orig_handle, fontsize, handlebox):
                scale = fontsize / 22
                x0, y0 = handlebox.xdescent, handlebox.ydescent
                width, height = handlebox.width, handlebox.height
                patch_line = mpatches.Rectangle([x0, y0+height/2], width, height/10,
                        color=self.color, transform=handlebox.get_transform())
                patch_circ = mpatches.Circle([x0 + width, y0 + height/2], height/2 * scale,
                        color=self.color, transform=handlebox.get_transform())
                patch_circ_small = mpatches.Circle([x0, y0 + height/2], height/5 * scale,
                        color=self.color, transform=handlebox.get_transform())
                patch_circ_small2 = mpatches.Circle([x0 + width/2, y0 + height/2], height/5 * scale,
                        color=self.color, transform=handlebox.get_transform())

                handlebox.add_artist(patch_line)
                handlebox.add_artist(patch_circ)
                handlebox.add_artist(patch_circ_small)
                handlebox.add_artist(patch_circ_small2)
                return patch_line

        class data_handler2(object):
            def __init__(self, color):
                self.color = color
            def legend_artist(self, legend, orig_handle, fontsize, handlebox):
                scale = fontsize / 22
                x0, y0 = handlebox.xdescent, handlebox.ydescent
                width, height = handlebox.width, handlebox.height
                patch_line = mpatches.Rectangle([x0, y0+height/2], width, height/10,
                        color=self.color, transform=handlebox.get_transform())
                patch_circ = mpatches.Circle([x0 + width, y0 + height/2], height/2 * scale,
                        color=self.color, transform=handlebox.get_transform())

                handlebox.add_artist(patch_line)
                handlebox.add_artist(patch_circ)
                return patch_line

        star = Line2D([0], [0], marker='*', color='w', label='Scatter',
                            markerfacecolor='orange', markersize=20)

        plt.legend([ObjR(), ObjG(), ObjB(), star], ['$t_{(w)}\in N_{\mathcal{T}}(q)\in \mathcal{C}$', '$N_{\mathcal{T}}(t)$' ,'$N_{\mathcal{Q}}(t)$', 'ERP Gantries'], handler_map={ObjR: data_handler2('r'), ObjB: data_handler('b'), ObjG: data_handler('g')})

    for data in data_taxi_morning:
        for obs_json in data:
            plot_data(obs_json, is_plot_density=False, is_plot_dist=False, is_plot_q=False, qalpha=0.5, is_plot_obst_traj=False)


    plotERP(cm='*', markersize=15, color='orange', zorder=3)
    plt.xlim(103.75, 104.)
    plt.ylim(1.26, 1.4)
    plotLegendTaxi()
    # plt.title('Taxi:morning')

    #rectA
    rectA = Rectangle((103.777, 1.34263),0.03,0.021,linewidth=2,edgecolor='k',facecolor='none', zorder=4)
    plt.gca().add_patch(rectA)
    plt.text(103.777+0.012, 1.364, 'A', fontsize=16)

    #rectB
    rectB = Rectangle((103.862, 1.371),0.03,0.021,linewidth=2,edgecolor='k',facecolor='none', zorder=4)
    plt.gca().add_patch(rectB)
    plt.text(103.862+0.012, 1.363, 'B', fontsize=16)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=1.0)
    plt.savefig('case_study_taxi_morning.pdf')
    plt.savefig('case_study_taxi_morning.png')
    plt.show()

def plot_vessel_case_study():
    filename = 'res/vessel1_obst_delta1.5_tau1.96.json'
    # filename = 'obst_test_subt_density_sgtaxi.json'
    with open(filename, 'r') as f:
        data_vessel1 = json.load(f)['obst']

    def plotLegendVessel():    
        class ObjB(object):
            pass
        class ObjR(object):
            pass
        class ObjG(object):
            pass

        class data_handler(object):
            def __init__(self, color):
                self.color = color
            def legend_artist(self, legend, orig_handle, fontsize, handlebox):
                scale = fontsize / 22
                x0, y0 = handlebox.xdescent, handlebox.ydescent
                width, height = handlebox.width, handlebox.height
                patch_line = mpatches.Rectangle([x0, y0+height/2], width, height/10,
                        color=self.color, transform=handlebox.get_transform())
                patch_circ = mpatches.Circle([x0 + width, y0 + height/2], height/2 * scale,
                        color=self.color, transform=handlebox.get_transform())
                patch_circ_small = mpatches.Circle([x0, y0 + height/2], height/5 * scale,
                        color=self.color, transform=handlebox.get_transform())
                patch_circ_small2 = mpatches.Circle([x0 + width/2, y0 + height/2], height/5 * scale,
                        color=self.color, transform=handlebox.get_transform())

                handlebox.add_artist(patch_line)
                handlebox.add_artist(patch_circ)
                handlebox.add_artist(patch_circ_small)
                handlebox.add_artist(patch_circ_small2)
                return patch_line

        class data_handler2(object):
            def __init__(self, color):
                self.color = color
            def legend_artist(self, legend, orig_handle, fontsize, handlebox):
                scale = fontsize / 22
                x0, y0 = handlebox.xdescent, handlebox.ydescent
                width, height = handlebox.width, handlebox.height
                patch_line = mpatches.Rectangle([x0, y0+height/2], width, height/10,
                        color=self.color, transform=handlebox.get_transform())
                patch_circ = mpatches.Circle([x0 + width, y0 + height/2], height/2 * scale,
                        color=self.color, transform=handlebox.get_transform())

                handlebox.add_artist(patch_line)
                handlebox.add_artist(patch_circ)
                return patch_line

        star =  Patch(facecolor='w', edgecolor='orange',
                            label='Operating Region')

        plt.legend([ObjR(), ObjG(), ObjB(), star], ['$t_{(w)}\in N_{\mathcal{T}}(q)\in \mathcal{C}$', '$N_{\mathcal{T}}(t)$' ,'$N_{\mathcal{Q}}(t)$', 'Operating Region'], handler_map={ObjR: data_handler2('r'), ObjB: data_handler('b'), ObjG: data_handler('g')}, 
                ncol=1)

    fig = plt.figure(figsize=(7, 4.5))
    for data_json in data_vessel1:
        for obs_json in data_json:
            plot_data(obs_json, is_plot_density=False, is_plot_dist=False, is_plot_q=False, qalpha=0.5, is_plot_obst_traj=False)
    plot_chosen_area(gt_vessel1, color='orange', zorder=4, linewidth=2)
    plt.xlim(103.8, 103.95)
    plt.ylim(1.15, 1.25)
    plotLegendVessel()
    # plt.title('Vessel:THORCO CLOUD')
    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=1.0)
    plt.savefig('case_study_vessel1.pdf')
    plt.savefig('case_study_vessel1.png')
    plt.show()

plot_erp_case_study()
# plot_vessel_case_study()