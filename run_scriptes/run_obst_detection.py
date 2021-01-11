from run_config import *
from util import create_results_folder, list_expansion
import os
from time import gmtime, strftime






@list_expansion
def run_alg(dataset, curtime=None, binary_name='./obst_test'):
    for param in dataset.for_params():
        cmd = '{} {}'.format(binary_name, param)
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    create_results_folder('../res')
    # run_alg([SGTAXI_MORNING(deltas=[1.25, 1.75, 2.25]), SGTAXI_AFTERNOON(deltas=[1.25, 1.75, 2.25])])
    # run_alg([VESSEL0(sigmas = [0.1], deltas=[1.5], taus=[1000, 2000, 3000]), VESSEL12(sigmas = [0.1], deltas=[1.5], taus=[1000, 2000, 3000])])
    # run_alg([SGTAXI_MORNING(sigmas = [0.02], deltas=[2], taus=[100])])
    # run_alg([VESSEL0(sigmas=[0.05, 0.1, 0.15], taus=[500, 1000, 1500], deltas=[1.5, 2, 2.5])])
    # run_alg([VESSEL12(sigmas=[0.05, 0.1, 0.15], taus=[500, 1000, 1500], deltas=[1.5, 2, 2.5])])
    # run_alg([VESSEL0()])

    
    # run_alg([SGTAXI_MORNING(sigmas=[0.015, 0.025], taus=[50], deltas=[2.])])
    # run_alg([SGTAXI_MORNING(sigmas=[0.01], taus=[75, 125], deltas=[2.])])
    # run_alg([SGTAXI_MORNING(sigmas=[0.01], taus=[50], deltas=[1.75, 2.25])])
    
    # run_alg([SGTAXI_AFTERNOON(sigmas=[0.015, 0.025], taus=[50], deltas=[2.])])
    # run_alg([SGTAXI_AFTERNOON(sigmas=[0.01], taus=[75, 125], deltas=[2.])])
    # run_alg([SGTAXI_AFTERNOON(sigmas=[0.01], taus=[50], deltas=[1.75, 2.25])])

    run_alg([SGTAXI_MORNING(sigmas=[0.005, 0.01, 0.015, 0.02, 0.025, 0.03], taus=[50], deltas=[2.])])
    run_alg([SGTAXI_MORNING(sigmas=[0.01], taus=[25, 50, 75, 100, 125, 150], deltas=[2.])])
    run_alg([SGTAXI_MORNING(sigmas=[0.01], taus=[50], deltas=[1.25, 1.5, 1.75, 2., 2.25, 2.5])])
    
    run_alg([SGTAXI_AFTERNOON(sigmas=[0.005, 0.01, 0.015, 0.02, 0.025, 0.03], taus=[50], deltas=[2.])])
    run_alg([SGTAXI_AFTERNOON(sigmas=[0.01], taus=[25, 50, 75, 100, 125, 150], deltas=[2.])])
    run_alg([SGTAXI_AFTERNOON(sigmas=[0.01], taus=[50], deltas=[1.25, 1.5, 1.75, 2., 2.25, 2.5])])