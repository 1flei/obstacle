import matplotlib.pylab as plt
import numpy as np
import sympy as sp
import math
from scipy.special import erf
from itertools import product, cycle
from scipy.stats import norm

#plot the collision probability and the gaussian pdf


plt.rcParams.update({'font.size': 16})

def gaussian(x):
    return np.exp(-x**2*0.5)

def laplacian(x):
    return np.exp(-np.abs(x))


def phi(x):
    #cdf of normal distirbution
    return (1.0 + erf(x / np.sqrt(2.0))) / 2.0

# def collision_prob_l2(c, r):
#     x = r/c
#     return 1 - 2*phi(-x) - 2/((2*np.pi)**0.5*x) * (1-np.exp(-x**2/2))


def collision_prob_l2(d, r):
    return 1-2*norm.cdf(-r/d)-2/(2**0.5*np.pi**0.5*(r/d))*(1-np.exp(-(r/d)**2/2))


def plog_collision_prob_l2_kl(r, k=1, l=1, c='b', m='', label=None):
    x = np.arange(0, 5, step=0.001)
    y = collision_prob_l2(x, r)
    y = 1-(1-(y**k))**l
    if label is None:
        label='r=%.2f_k=%d_l=%d'%(r, k, l)
    plt.plot(x, y, label=label, color=c, marker=m, markevery=100)
    plt.xlim(0, 5)
    plt.ylim(0, 1)
    # plt.show()

def collision_prob_l1(c, r):
    return 2*(np.arctan(r/c)/np.pi) - 1/np.pi/(r/c) *np.log(1+(r/c)**2)

def plog_collision_prob_l1_kl(r, k=1, l=1, c='b', m='.'):
    x = np.arange(0, 5, step=0.001)
    y = collision_prob_l1(x, r)
    y = 1-(1-(y**k))**l
    plt.plot(x, y, label='r=%.2f_k=%d_l=%d'%(r, k, l), color=c, marker=m, markevery=100)
    plt.xlim(0, 5)
    plt.ylim(0, 1)
    # plt.show()
    

def plot_guassian(c='r', m='x'):
    x = np.arange(0, 5, step=0.001)
    y = gaussian(x)
    plt.plot(x, y, label='exp(-x^2/2)', color=c, marker=m, markevery=100)
    plt.xlim(0, 5)
    # plt.show()

def plot_laplacian(c='r', m='x'):
    x = np.arange(0, 5, step=0.001)
    y = laplacian(x)
    plt.plot(x, y, label='exp(-|x|)', color=c, marker=m, markevery=100)
    plt.xlim(0, 5)


def find_param(rs = np.arange(0.001, 10, 0.001), ks=list(range(1, 16)), ls=list(range(1, 32)) ):
    x = np.arange(0, 10, step=0.001)
    y = gaussian(x)

    ress = []
    
    for r in rs:
        y_hash = collision_prob_l2(x, r)

        for k in ks:
            for l in ls:
                y_ = 1-(1-(y_hash**k))**l

                err = np.sum((y-y_)**2)

                ress += [[r, k, l, err]]

    sorted_ress = sorted(ress, key=lambda x:x[3])
    print(sorted_ress[: 10])


def find_param_laplacian(rs = np.arange(0.001, 10, 0.001), ks=list(range(1, 16)), ls=list(range(1, 32)) ):
    x = np.arange(0, 10, step=0.001)
    y = laplacian(x)

    ress = []
    
    for r in rs:
        y_hash = collision_prob_l1(x, r)

        for k in ks:
            for l in ls:
                y_ = 1-(1-(y_hash**k))**l

                err = np.nansum((y-y_)**2)

                ress += [[r, k, l, err]]

    sorted_ress = sorted(ress, key=lambda x:x[3])
    print(sorted_ress[: 10])

#remarks
#l1-laplacion: 5.008, 4, 2
#l2-gaussian


if __name__ == '__main__':
    # find_param(rs = np.arange(0.001, 10, 0.001), ls=list(range(1, 5) ) )
    # find_param(rs = np.arange(0.001, 10, 0.01), ks=range(1, 3), ls=[1] )
    # find_param_laplacian(rs = np.arange(0.001, 10, 0.001), ls=list(range(1, 5)))

    # plot_laplacian()
    # plog_collision_prob_l1_kl(5.008, 4, 2, c='b')
    
    plt.rcParams.update({'font.size': 24})
    plot_guassian()
    plog_collision_prob_l2_kl(5.123, 8, 3, c='b', label='PIE Approach')
    plt.xlim(0, 5)
    plt.ylim(-0.0, 1.0)
    plt.xticks([0, 1, 2, 3, 4, 5])
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.])
    plt.xlabel(r'$\Delta$')
    plt.ylabel(r'$\kappa_d(\Delta)$')
    plt.legend(loc='upper right')
    plt.savefig('kernel_pie.png', bbox_inches='tight')
    # plot_guassian()
    # plog_collision_prob_l2_kl(1., 1, 1, c='b')

    # ds = np.arange(0, 5, step=0.001)
    # # ys = erf(1/ds/2**0.5)
    # ys = 1-2*norm.cdf(-1/ds)
    # ys2 = 2/(2**0.5*np.pi**0.5*(1/ds))*(1-np.exp(-(1/ds)**2/2))
    # plt.plot(ds, ys, 'g-')
    # plt.plot(ds, ys2, 'c-')
    # markers = ['o', 'v', 's', 'p']
    # colors = ['g', 'b', 'c', 'm', 'y']

    # cm = cycle(product(colors, markers) )
    # for r in [3, 4, 5, 6]:
    #     for k in [7]:
    #         for l in [3, 4, 5]:
    #             c, m = next(cm)
    #             plog_collision_prob_l2_kl(r=r, k=k, l=l, c=c, m=m)
    # plog_collision_prob_l2_kl(r=1, k=1, l=1, c='b')
    # plog_collision_prob_l2_kl(r=3, k=5, l=4, c='g')
    plt.show()