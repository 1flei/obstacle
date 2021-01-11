from run_config import *
from util import create_results_folder, list_expansion
import os
from time import gmtime, strftime




class SGTAXI_KDE_ANOMALY:
    def __init__(self, thresholds= [5, 10, 15, 20]):
        self.name = 'sgtaxi_morning'
        self.ref_name = 'sgtaxi_ref'
        self.query_name = 'sgtaxi_ref'
        self.thresholds = thresholds

        self.hnsw_max_objs = 300000

    
    def get_output_filename(self, threshold, curtime=''):
        return '../res/{}_[{}_threshold{}].json'.format(self.name, curtime, threshold)

    def for_params(self, curtime=None):
        if curtime is None:
            curtime = strftime("%m-%d_%H_%M", gmtime())
        for threshold in self.thresholds:
            yield r'-D{} -Q{}  --n_thread_building_index 1 --hnsw_max_objs {} -O{} --sigma 0.01 --threshold {}' \
                .format(self.ref_name, self.query_name, self.hnsw_max_objs, 
                self.get_output_filename(threshold, curtime=curtime), threshold)

class VESSEL0:
    def __init__(self, thresholds= [20, 40, 60, 80]):
        self.name = 'vessel0'
        self.ref_name = 'interp_may_jun'
        self.query_name = 'interp_may_jun'
        self.time_interval = 600
        self.thresholds = thresholds

        self.hnsw_max_objs = 6000000

    
    def get_output_filename(self, threshold, curtime=''):
        return '../res/{}_[{}_threshold{}].json'.format(self.name, curtime, threshold)


    def for_params(self, curtime=None):
        if curtime is None:
            curtime = strftime("%m-%d_%H_%M", gmtime())
        for threshold in self.thresholds:
            yield r'-D{} -Q{}  --time_interval {} --n_thread_building_index 1 --hnsw_max_objs {} -O{} --sigma 0.001 --threshold {}' \
                .format(self.ref_name, self.query_name, self.time_interval, self.hnsw_max_objs, 
                self.get_output_filename(threshold, curtime=curtime), threshold)




@list_expansion
def run_alg(dataset, curtime=None, binary_name='./kde_anomaly_test'):
    for param in dataset.for_params():
        cmd = '{} {} '.format(binary_name, param)
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    create_results_folder('../res')
    # run_alg([SGTAXI_KDE_ANOMALY(), VESSEL0()])
    run_alg([SGTAXI_KDE_ANOMALY(thresholds=[1.])])