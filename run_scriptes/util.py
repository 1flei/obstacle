import os




def create_results_folder(path='../res'):
    try: 
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise


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

