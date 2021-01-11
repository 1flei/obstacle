from plot_collision_prob import *
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np



def plot_guassian_3d(c='r', m='x'): 
    # ax = plt.gca(projection='3d')

    xx = np.arange(0, 5, step=0.01)
    yy = np.arange(0, 5, step=0.01)
    x, y = np.meshgrid(xx, yy)
    z = gaussian(x) * gaussian(y)
    # plt.show()
    plt.contour(x, y, z, 5, colors='r')


def plog_collision_prob_l2_kl_3d(r, k=1, l=1, c='b', m='', label=None):
    # ax = plt.gca(projection='3d')

    xx = np.arange(0, 5, step=0.01)
    yy = np.arange(0, 5, step=0.01)
    x, y = np.meshgrid(xx, yy)
    px = 1-(1-collision_prob_l2(x, r)**k)**l
    py = 1-(1-collision_prob_l2(y, r)**k)**l
    z = px*py
    if label is None:
        label='r=%.2f_k=%d_l=%d'%(r, k, l)
    plt.contour(x, y, z, 5, colors='b')


def plog_collision_prob_diff(r, k=1, l=1, c='b', m='', label=None):
    fig = plt.figure()  
    ax = plt.gca(projection='3d')

    xx = np.arange(0, 5, step=0.01)
    yy = np.arange(0, 5, step=0.01)
    x, y = np.meshgrid(xx, yy)
    px = 1-(1-collision_prob_l2(x, r)**k)**l
    py = 1-(1-collision_prob_l2(y, r)**k)**l

    z = px*py - gaussian(x) * gaussian(y)

    if label is None:
        label='r=%.2f_k=%d_l=%d'%(r, k, l)
    # ax.set_zlim(-0.1, 0.1)
    ax.set_zlim(-0.05, 0.05)
    ax.set_zticks([-0.05, 0, 0.05])
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 5)

    ax.set_xlabel('$\Delta_1$')
    ax.set_ylabel('$\Delta_2$')
    ax.set_zlabel(r'$\kappa_d(\Delta_1)\kappa_d(\Delta_2)-p_{succ}$')
    ax.xaxis.labelpad=10
    ax.yaxis.labelpad=10
    ax.zaxis.labelpad=22
    ax.tick_params(axis='z', which='major', pad=8)

    fig.subplots_adjust(left=0.05)
    # fig.subplots_adjust(right=0.9)
    ax.plot_wireframe(x, y, z, color=c, rcount=10, ccount=10)

    # plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.1)

    ax.legend()

if __name__ == '__main__':
    plt.rcParams.update({'font.size': 22})
    # plot_guassian_3d()
    plog_collision_prob_diff(5.123, 8, 3, c='b', label='PIE-based approach')

    plt.savefig('kernel_pie_3d.png')
    plt.show()