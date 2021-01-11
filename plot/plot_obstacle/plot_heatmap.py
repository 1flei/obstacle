import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv
import json 
from math import floor
from plot_obst_json_util import getERPAtTime

mydpi = 100

def plot_tata():
    # filename = 'tata_json.json'
    # output_filename = 'tata_json.png'
    filename = 'tata_json_aug.json'
    output_filename = 'tata_json_aug.png'
    with open(filename, 'r') as f:
        tata_ref = json.load(f)

    # print(tata_ref)

    grid_width = 0.0002
    arr_width = int((103.71-103.42)/grid_width)+2
    arr_height = int((1.2-1.0)/grid_width)+2
    print(arr_width, arr_height)

    def grid_fx(x_grid_idx):
        x = (x_grid_idx - floor(103.42/grid_width+0.5) )
        return x
    def grid_fy(y_grid_idx):
        y = arr_height -1 - (y_grid_idx - floor(1.0/grid_width+0.5) )
        return y

    def get_arr(j):
        arr = np.zeros(shape=(arr_height, arr_width))
        for items in j:
            #since y=0 is the top pix
            x, y = grid_fx(items[0][0]), grid_fy(items[0][1])
            cnt = items[1]
            arr[y, x] = cnt
        return arr

    arr = get_arr(tata_ref)
    max_cnt = np.max(arr)

    cnt97 = np.quantile(arr, 0.97)
    print(arr, 'max_cnt=', max_cnt, 'cnt@95=', cnt97)
    arr[arr>cnt97] = cnt97


    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0.005, right=0.995, bottom=0, top=1)
#     fig, ax = plt.subplots(facecolor='white')
    fig.patch.set_facecolor('white')
#     m = plotBaseMap(mindata, maxdata)
    fig.set_size_inches(arr_width/mydpi, arr_height/mydpi)

#     ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])

    pcm = ax.imshow(arr, cmap='hot')

    chosen_area = np.array([(1.113583, 103.5263),(1.06417, 103.6067),(1.0482, 103.596616),(1.09683, 103.516183),(1.113583, 103.5263)])
    
    def plot_chosen_area():
        data = np.floor(chosen_area/grid_width + 0.5)
        data[:, 0] = grid_fy(data[:, 0])
        data[:, 1] = grid_fx(data[:, 1])
        print(data)

        for p, pp in zip(data, data[1:]):
            plt.plot([p[1], pp[1]], [p[0], pp[0]], 'b-')

    plot_chosen_area()
    plt.savefig(output_filename, dpi=mydpi)
    plt.show()

def plot_sunken_vessel0():
    grid_width = 0.0002
    # // [104.10, 104.25] * [1.25, 1.35]
    x_max = 104.25
    x_min = 104.10
    y_max = 1.35
    y_min = 1.25
    chosen_area = np.array([
        (1.286, 104.183), 
        (1.2807, 104.18317), 
        (1.28, 104.17683), 
        (1.2855, 104.17617), 
        (1.286, 104.183)
    ])
    
    arr_width = int((x_max-x_min)/grid_width)+2
    arr_height = int((y_max-y_min)/grid_width)+2

    filename = 'sunken_vessel0.json'
    output_filename = 'sunken_vessel0.png'
    # filename = 'sunken_vessel0_query.json'
    # output_filename = 'sunken_vessel0_query.png'
    with open(filename, 'r') as f:
        tata_ref = json.load(f)

    # print(tata_ref)
    # [104.15, 104.20] * [1.27, 1.30]

    print(arr_width, arr_height)

    def grid_fx(x_grid_idx):
        x = (x_grid_idx - floor(x_min/grid_width+0.5) )
        return x
    def grid_fy(y_grid_idx):
        y = arr_height -1 - (y_grid_idx - floor(y_min/grid_width+0.5) )
        return y

    def get_arr(j):
        arr = np.zeros(shape=(arr_height, arr_width))
        for items in j:
            #since y=0 is the top pix
            x, y = grid_fx(items[0][0]), grid_fy(items[0][1])
            cnt = items[1]
            arr[y, x] = cnt
        return arr

    arr = get_arr(tata_ref)
    max_cnt = np.max(arr)

    cnt97 = np.quantile(arr, 0.97)
    print(arr, 'max_cnt=', max_cnt, 'cnt@95=', cnt97)
    arr[arr>cnt97] = cnt97


    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0.005, right=0.995, bottom=0, top=1)
#     fig, ax = plt.subplots(facecolor='white')
    fig.set_size_inches(arr_width/mydpi, arr_height/mydpi)
    fig.patch.set_facecolor('white')
#     m = plotBaseMap(mindata, maxdata)

#     ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])

    pcm = ax.imshow(arr, cmap='hot')

    def plot_chosen_area():
        data = np.floor(chosen_area/grid_width + 0.5)
        data[:, 0] = grid_fy(data[:, 0])
        data[:, 1] = grid_fx(data[:, 1])
        print(data)

        for p, pp in zip(data, data[1:]):
            plt.plot([p[1], pp[1]], [p[0], pp[0]], 'b-')

    plot_chosen_area()
    plt.savefig(output_filename, dpi=mydpi)
    plt.show()

def plot_heatmap(grid_width = 0.0002, x_max = 104.25, x_min = 104.10, 
    y_max = 1.35, y_min = 1.25, 
    chosen_area = np.array([
        (1.286, 104.183), 
        (1.2807, 104.18317), 
        (1.28, 104.17683), 
        (1.2855, 104.17617), 
        (1.286, 104.183)
    ]), 
    filenames = ['sunken_vessel0.json', 'sunken_vessel0_query.json'], 
    output_filenames = ['sunken_vessel0.png', 'sunken_vessel0_query.png'], 
    pratio = 0.97
):
    
    # // [104.10, 104.25] * [1.25, 1.35]
    for filename, output_filename in zip(filenames, output_filenames):
        arr_width = int((x_max-x_min)/grid_width)+2
        arr_height = int((y_max-y_min)/grid_width)+2
        with open(filename, 'r') as f:
            tata_ref = json.load(f)

        # print(tata_ref)
        # [104.15, 104.20] * [1.27, 1.30]

        print(arr_width, arr_height)

        def grid_fx(x_grid_idx):
            x = (x_grid_idx - floor(x_min/grid_width+0.5) )
            return x
        def grid_fy(y_grid_idx):
            y = arr_height -1 - (y_grid_idx - floor(y_min/grid_width+0.5) )
            return y

        def get_arr(j):
            arr = np.zeros(shape=(arr_height, arr_width))
            for items in j:
                #since y=0 is the top pix
                x, y = grid_fx(items[0][0]), grid_fy(items[0][1])
                cnt = items[1]
                if items[0][0]<0:
                    continue
                # print(items, y, x)
                arr[y, x] = cnt
            return arr

        arr = get_arr(tata_ref)
        max_cnt = np.max(arr)

        cnt97 = np.quantile(arr, pratio)
        print(arr, 'max_cnt=', max_cnt, 'cnt@95=', cnt97)
        arr[arr>cnt97] = cnt97


        fig, ax = plt.subplots()
        fig.subplots_adjust(left=0.005, right=0.995, bottom=0, top=1)
    #     fig, ax = plt.subplots(facecolor='white')
        fig.set_size_inches(arr_width/mydpi, arr_height/mydpi)
        fig.patch.set_facecolor('white')
    #     m = plotBaseMap(mindata, maxdata)

    #     ax.axis('off')
        ax.set_xticks([])
        ax.set_yticks([])
        pcm = ax.imshow(arr, cmap='Blues')
        plt.axis('off')

        def plot_chosen_area():
            data = np.floor(chosen_area/grid_width + 0.5)
            data[:, 0] = grid_fy(data[:, 0])
            data[:, 1] = grid_fx(data[:, 1])
            print(data)

            for p, pp in zip(data, data[1:]):
                plt.plot([p[1], pp[1]], [p[0], pp[0]], 'r-')

        plot_chosen_area()
        plt.savefig(output_filename, dpi=mydpi)
        # plt.show()
    

            


def plot_heatmap_taxi(grid_width = 0.0002, x_max = 104.25, x_min = 104.10, 
    y_max = 1.35, y_min = 1.25, 
    filenames = ['sunken_vessel0.json', 'sunken_vessel0_query.json'], 
    output_filenames = ['sunken_vessel0.png', 'sunken_vessel0_query.png'], 
    pratio = 0.97, 
    working_time = None
):
    
    # // [104.10, 104.25] * [1.25, 1.35]
    for filename, output_filename in zip(filenames, output_filenames):
        arr_width = int((x_max-x_min)/grid_width)+2
        arr_height = int((y_max-y_min)/grid_width)+2
        with open(filename, 'r') as f:
            tata_ref = json.load(f)

        # print(tata_ref)
        # [104.15, 104.20] * [1.27, 1.30]

        print(arr_width, arr_height)

        def grid_fx(x_grid_idx):
            x = (x_grid_idx - floor(x_min/grid_width+0.5) )
            return x
        def grid_fy(y_grid_idx):
            y = arr_height -1 - (y_grid_idx - floor(y_min/grid_width+0.5) )
            return y

        def get_arr(j):
            arr = np.zeros(shape=(arr_height, arr_width))
            for items in j:
                #since y=0 is the top pix
                x, y = grid_fx(items[0][0]), grid_fy(items[0][1])
                cnt = items[1]
                if items[0][0]<0:
                    continue
                arr[y, x] = cnt
            return arr

        arr = get_arr(tata_ref)
        max_cnt = np.max(arr)

        cnt97 = np.quantile(arr, pratio)
        print(arr, 'max_cnt=', max_cnt, 'cnt@95=', cnt97)
        arr[arr>cnt97] = cnt97


        fig, ax = plt.subplots()
        fig.subplots_adjust(left=0.005, right=0.995, bottom=0, top=1)
    #     fig, ax = plt.subplots(facecolor='white')
        fig.set_size_inches(arr_width/mydpi, arr_height/mydpi)
        fig.patch.set_facecolor('white')
    #     m = plotBaseMap(mindata, maxdata)

    #     ax.axis('off')
        ax.set_xticks([])
        ax.set_yticks([])
        pcm = ax.imshow(arr, cmap='Blues')
        plt.axis('off')

        def plotERP_heatmap(t='08:00:01', fname = 'erp2.json', cm='b*', alpha=1., markersize=15):
            erps = getERPAtTime(t, fname)
            erps = np.array(erps/grid_width + 0.5)
            plt.plot(grid_fx(erps[:, 0]), grid_fy(erps[:, 1]), cm, markersize=markersize, alpha=alpha, label='working ERPs')

        if working_time is not None:
            plotERP_heatmap(cm='r*', t=working_time, markersize=15)
        plt.savefig(output_filename, dpi=mydpi)
        # plt.show()



if __name__=='__main__':
    # plot_tata()
    # # HARITA BERLIAN 18
    # plot_heatmap(grid_width = 0.0002, x_max = 104.25, x_min = 104.10, 
    #     y_max = 1.35, y_min = 1.25, 
    #     chosen_area = np.array([
    #         (1.286, 104.183), 
    #         (1.2807, 104.18317), 
    #         (1.28, 104.17683), 
    #         (1.2855, 104.17617), 
    #         (1.286, 104.183)
    #     ]), 
    #     filenames = ['sunken_vessel_0_ref.json', 'sunken_vessel_0_q.json'], 
    #     output_filenames = ['sunken_vessel_0.png', 'sunken_vessel_0_query.png'], 
    #     pratio = 0.97
    # )

    # # THORCO CLOUD
    # plot_heatmap(grid_width = 0.0002, x_max = 104.00, x_min =103.80, 
    #     y_max = 1.25, y_min = 1.15, 
    #     chosen_area = np.array([
    #         (1.204300, 103.898650), 
    #         (1.208133, 103.901667), 
    #         (1.217517, 103.917367), 
    #         (1.213350, 103.917333), 
    #         (1.204500, 103.902600), 
    #         (1.204300, 103.898650), 
    #     ]), 
    #     filenames = ['sunken_vessel_1_may.json', 'sunken_vessel_1_aug.json'], 
    #     output_filenames = ['sunken_vessel_1_may.png', 'sunken_vessel_1_aug.png'], 
    #     pratio = 0.97
    # )

    # Cai Jun 3
    # plot_heatmap(grid_width = 0.0002, x_max = 104.55, x_min =104.35, 
    #     y_max = 1.5, y_min = 1.35, 
    #     chosen_area = np.array([
    #         (1.431917, 104.456483), 
    #         (1.437067, 104.459750), 
    #         (1.426917, 104.462750), 
    #         (1.431617, 104.456517), 
    #         (1.431917, 104.456483)
    #     ]), 
    #     filenames = ['sunken_vessel_2_jun.json', 'sunken_vessel_2_aug.json'], 
    #     output_filenames = ['sunken_vessel_2_jun.png', 'sunken_vessel_2_aug.png'], 
    #     pratio = 0.97
    # )

    # Putri Sea, 104.05, 104.25, 1.2, 1.35
    # plot_heatmap(grid_width = 0.0002, x_max = 104.25, x_min =104.05, 
    #     y_max = 1.35, y_min = 1.2, 
    #     chosen_area = np.array([
    #         (1.297500, 104.120067), 
    #         (1.297500, 104.120667), 
    #         (1.298500, 104.120667), 
    #         (1.298500, 104.120067), 
    #         (1.297500, 104.120067)
    #     ]), 
    #     filenames = ['sunken_vessel_3_may.json', 'sunken_vessel_3_aug.json'], 
    #     output_filenames = ['sunken_vessel_3_may.png', 'sunken_vessel_3_aug.png'], 
    #     pratio = 0.97
    # )

    # THORCO CLOUD's heatmap from interpolated data
    # plot_heatmap(grid_width = 0.0002, x_max = 104.00, x_min =103.80, 
    #     y_max = 1.25, y_min = 1.15, 
    #     chosen_area = np.array([
    #         (1.204300, 103.898650), 
    #         (1.208133, 103.901667), 
    #         (1.217517, 103.917367), 
    #         (1.213350, 103.917333), 
    #         (1.204500, 103.902600), 
    #         (1.204300, 103.898650), 
    #     ]), 
    #     filenames = ['interp_may1_hm.json', 'interp_aug_hm.json'], 
    #     output_filenames = ['interp_may1_hm.png', 'interp_aug_hm.png'], 
    #     pratio = 0.97
    # )

    # sgtaxi heatmap
    # plot_heatmap_taxi(grid_width = 0.0002, x_min =103.68, x_max = 104.00, 
    #     y_min = 1.23, y_max = 1.47, 
    #     filenames = ['res/sgtaxi_heatmap.json', 'res/sgtaxi_heatmap_afternoon.json'], 
    #     output_filenames = ['sgtaxi_heatmap_afternoon_ref.png', 'sgtaxi_heatmap_afternoon.png'], 
    #     pratio = 0.97, 
    #     working_time = '18:00:01'
    # )
    # sgtaxi heatmap
    plot_heatmap_taxi(grid_width = 0.0002, x_min =103.675, x_max = 104.00, 
        y_min = 1.23, y_max = 1.47, 
        filenames = ['res/sgtaxi_heatmap.json', 'res/sgtaxi_heatmap_afternoon.json'], 
        output_filenames = ['sgtaxi_heatmap_morning_ref.png', 'sgtaxi_heatmap_morning.png'], 
        pratio = 0.97, 
        working_time = '08:00:01'
    )