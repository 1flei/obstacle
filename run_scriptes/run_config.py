from time import gmtime, strftime

class SGTAXI_MORNING:
    def __init__(self, sigmas = [0.01, 0.02, 0.03, 0.04], deltas=[1.5, 2, 2.5], taus=[50, 100, 150]):
        self.sigmas = sigmas
        self.deltas = deltas
        self.taus = taus
        self.name = 'sgtaxi_morning'
        self.ref_name = 'sgtaxi_ref'
        self.query_name = 'sgtaxi_morning_q'

        self.hnsw_max_objs = 300000

    
    def get_output_filename(self, sigma, delta, tau, curtime=''):
        return '../res/{}_[{}_sigma{}_delta{}_tau{}].json'.format(self.name, curtime, sigma, delta, tau)

    def for_params(self, curtime=None):
        if curtime is None:
            curtime = strftime("%m-%d_%H_%M", gmtime())
        for sigma in self.sigmas:
            for delta in self.deltas:
                for tau in self.taus:
                    yield r'-D{} -Q{} --sigma {} --significance {} --min_support {}  --n_thread_building_index 1 --hnsw_max_objs {} -O{}' \
                        .format(self.ref_name, self.query_name, sigma, delta, tau, self.hnsw_max_objs, 
                        self.get_output_filename(sigma, delta, tau, curtime=curtime))

class SGTAXI_AFTERNOON:
    def __init__(self, sigmas = [0.01, 0.02, 0.03, 0.04], deltas=[1.5, 2, 2.5], taus=[50, 100, 150]):
        self.sigmas = sigmas
        self.deltas = deltas
        self.taus = taus
        self.name = 'sgtaxi_afternoon'
        self.ref_name = 'sgtaxi_ref'
        self.query_name = 'sgtaxi_afternoon_q'

        self.hnsw_max_objs = 300000

    
    def get_output_filename(self, sigma, delta, tau, curtime=''):
        return '../res/{}_[{}_sigma={}_delta={}_tau={}].json'.format(self.name, curtime, sigma, delta, tau)

    def for_params(self, curtime=None):
        if curtime is None:
            curtime = strftime("%m-%d_%H_%M", gmtime())
        for sigma in self.sigmas:
            for delta in self.deltas:
                for tau in self.taus:
                    yield r'-D{} -Q{} --sigma {} --significance {} --min_support {}  --n_thread_building_index 1 --hnsw_max_objs {} -O{}' \
                        .format(self.ref_name, self.query_name, sigma, delta, tau, self.hnsw_max_objs, 
                        self.get_output_filename(sigma, delta, tau, curtime=curtime))


class VESSEL0:
    def __init__(self, sigmas = [0.05, 0.1, 0.15, 0.2], deltas=[1.5, 2, 2.5], taus=[100, 200, 300]):
        self.sigmas = sigmas
        self.deltas = deltas
        self.taus = taus
        self.name = 'vessel0'
        self.ref_name = 'interp_may_jun'
        self.query_name = 'interp_aug'
        self.time_interval = 600

        self.hnsw_max_objs = 6000000

    
    def get_output_filename(self, sigma, delta, tau, curtime=''):
        return '../res/{}_[{}_sigma={}_delta={}_tau={}].json'.format(self.name, curtime, sigma, delta, tau)

    def for_params(self, curtime=None):
        if curtime is None:
            curtime = strftime("%m-%d_%H_%M", gmtime())
        for sigma in self.sigmas:
            for delta in self.deltas:
                for tau in self.taus:
                    yield r'-D{} -Q{}  --sigma {} --significance {} --min_support {} --time_interval {} --n_thread_building_index 1 --hnsw_max_objs {} -O{}' \
                        .format(self.ref_name, self.query_name, sigma, delta, tau, self.time_interval, self.hnsw_max_objs, 
                        self.get_output_filename(sigma, delta, tau, curtime=curtime))

class VESSEL12:
    def __init__(self, sigmas = [0.05, 0.1, 0.15, 0.2], deltas=[1.5, 2, 2.5], taus=[100, 200, 300]):
        self.sigmas = sigmas
        self.deltas = deltas
        self.taus = taus
        self.name = 'vessel1'
        self.ref_name = 'interp_aug'
        self.query_name = 'interp_may_jun'
        self.time_interval = 600

        self.hnsw_max_objs = 6000000
    
    def get_output_filename(self, sigma, delta, tau, curtime=''):
        return '../res/{}_[{}_sigma={}_delta={}_tau={}].json'.format(self.name, curtime, sigma, delta, tau)

    def for_params(self, curtime=None):
        if curtime is None:
            curtime = strftime("%m-%d_%H_%M", gmtime())
        for sigma in self.sigmas:
            for delta in self.deltas:
                for tau in self.taus:
                    yield r'-D{} -Q{}  --sigma {} --significance {} --min_support {} --time_interval {} --n_thread_building_index 1 --hnsw_max_objs {} -O{}' \
                        .format(self.ref_name, self.query_name, sigma, delta, tau, self.time_interval, self.hnsw_max_objs, 
                        self.get_output_filename(sigma, delta, tau, curtime=curtime))

# class VESSEL2:
#     def __init__(self, sigmas = [0.05, 0.1, 0.15, 0.2], deltas=[1.5, 2, 2.5], taus=[50, 100, 150]):
#         self.sigmas = sigmas
#         self.deltas = deltas
#         self.taus = taus
#         self.name = 'vessel2'
#         self.ref_name = 'interp_aug'
#         self.query_name = 'interp_aug_q'
#         self.time_interval = 600

    
#     def get_output_filename(self, sigma, delta, tau, curtime=''):
#         return '../res/{}_[{}_sigma={}_delta={}_tau={}].json'.format(self.name, curtime, sigma, delta, tau)

#     def for_params(self, curtime=None):
#         if curtime is None:
#             curtime = strftime("%m-%d_%H_%M", gmtime())
#         for sigma in self.sigmas:
#             for delta in self.deltas:
#                 for tau in self.taus:
#                     yield r'-D{} -Q{}  --sigma {} --significance {} --min_support {} --time_interval {} --n_thread_building_index 1 -O{}' \
#                         .format(self.ref_name, self.query_name, sigma, delta, tau, self.time_interval, 
#                         self.get_output_filename(sigma, delta, tau, curtime=curtime))


class VESSEL_DENSITY_DATASET:
    def __init__(self, sigma = 0.1):
        self.name = 'vessel'
        self.sigma = sigma


class LSH_DENSITY_PIE_CM:
    def __init__(self, ht_sizes=[4194301], n_repeats=[64], n_cm_repeats=[1]):
        self.name = 'lsh_density_pie'
        self.ht_sizes = ht_sizes
        self.n_repeats = n_repeats
        self.n_cm_repeats = n_cm_repeats

    def for_param(self):
        for ht_size in self.ht_sizes:
            for n_repeat in self.n_repeats:
                for n_cm_repeat in self.n_cm_repeats:
                    yield '--ht_size %d --n_repeat %d --n_cm_repeat %d'%(ht_size, n_repeat, n_cm_repeat)
                    
class LSH_DENSITY_LINEAR:
    def __init__(self, n_repeats=[64], n_cm_repeats=[1]):
        self.name = 'lsh_density_linear'
        self.n_repeats = n_repeats
        self.n_cm_repeats = n_cm_repeats

    def for_param(self):
        for n_repeat in self.n_repeats:
            for n_cm_repeat in self.n_cm_repeats:
                yield '--n_repeat %d --n_cm_repeat %d'%(n_repeat, n_cm_repeat)



                
