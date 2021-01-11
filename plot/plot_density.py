import matplotlib.pylab as plt
import numpy as np
import json


def read_data(filename):
    with open(filename, 'r') as f:
        j = json.load(f)

        ret = []
        for pnt in j:
            ret += [[pnt['lsh_cnt'], pnt['brute_force_kde']]]

        return np.array(ret)





if __name__ == '__main__':
    data = read_data('res/vessel0_density_test.json')

    plt.plot(data[:, 0], data[:, 1], 'b.')
    plt.show()