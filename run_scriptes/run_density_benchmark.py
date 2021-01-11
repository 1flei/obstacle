from run_config import *
import os
from time import gmtime, strftime

#a decorator
def list_expansion(f):
    def g(*args, **kwargs):
        args_copy = list(args)
        for i, arg in enumerate(args):
            if type(arg) == list:
                for element in arg:
                    args_copy[i] = element
                    g(*args_copy, **kwargs)
                return 
        f(*args)
    return g


def get_dataset_path(dataset):
    return '../density_test_dataset/%s.ds'%(dataset.name,)
def get_query_path(dataset):
    return '../density_test_dataset/%s.q'%(dataset.name,)

def get_grount_truth_path(dataset, sigma=0.001):
    return '../density_test_dataset/%s_sigma=%.3f.gt'%(dataset.name, sigma)

def create_results_folder(path='../density_test_res'):
    try: 
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise


def get_output_filename(dataset, method, curtime=''):
    return '../density_test_res/%s_%s_[%s].json'%(dataset.name, method.name, curtime)

def get_dataset_params(ds):
    return '--sigma %f -D %s -Q %s'%(ds.sigma, 
        get_dataset_path(ds), 
        get_query_path(ds))


def get_common_params(ds, method, sigma=0.1, curtime=None):
    if curtime is None:
        curtime = strftime("%m-%d_%H_%M", gmtime())
    return '-D%s -Q%s --%s -G%s --sigma %f -O%s --force_build_index'%\
        (get_dataset_path(ds), get_query_path(ds), method.name, get_grount_truth_path(ds, sigma), sigma, get_output_filename(ds, method, curtime))


# ./density_benchmarker -D../density_test_dataset/sgtaxi.ds -Q../density_test_dataset/sgtaxi.q --lsh_density_pie_cm  -G../density_test_dataset/sgtaxi_sigma=0.001.gt --sigma 0.001\
#     --ht_size $ht_size --n_repeat $n_repeat --n_cm_repeat $n_cm_repeat  \
#     -O../res_json/sgtaxi_pie_cm_h${ht_size}_nr${n_repeat}_ncr${n_cm_repeat}.json 
@list_expansion
def run_alg(dataset, method, sigma=0.1, curtime=None, binary_name='./density_benchmarker'):
    common_params = get_common_params(dataset, method, sigma=sigma)

    for method_param in method.for_param():
        cmd = '%s %s %s'%(binary_name, common_params, method_param)
        print(cmd)
        os.system(cmd)

@list_expansion
def compute_ground_truth(dataset, sigma=0.1, binary_name='./density_benchmarker'):
    files = '-D {} -Q {} -G {}'.format(get_dataset_path(dataset), get_query_path(dataset), get_grount_truth_path(dataset, sigma))
    cmd = '{} --compute_ground_truth --sigma {} {}'.format(binary_name, sigma, files)
    print(cmd)
    os.system(cmd)


if __name__ == '__main__':
    create_results_folder()
    # run_alg([SGTAXI()], [LSH_DENSITY_PIE_CM()])
    # run_alg([SGTAXI()], [LSH_DENSITY_PIE_HT()])
    # run_alg([SGTAXI()], [LSH_DENSITY_PIE_HT(), LSH_DENSITY_BKS()])
    # compute_ground_truth([VESSEL_DENSITY_DATASET()])
    # run_alg([VESSEL_DENSITY_DATASET()], [LSH_DENSITY_PIE_CM(), LSH_DENSITY_LINEAR()])
    run_alg([VESSEL_DENSITY_DATASET()], [LSH_DENSITY_LINEAR()])