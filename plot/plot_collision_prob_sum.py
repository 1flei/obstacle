import matplotlib.pylab as plt
import numpy as np
import sympy as sp
import math
from scipy.special import erf
from itertools import product, cycle
from scipy.stats import norm
from sklearn.linear_model import LinearRegression, Ridge

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


def plog_collision_prob_l2_kl(r, k=1, l=1, c='b', m='.'):
    x = np.arange(0, 5, step=0.001)
    y = collision_prob_l2(x, r)
    y = 1-(1-(y**k))**l
    plt.plot(x, y, label='r=%.2f_k=%d_l=%d'%(r, k, l), color=c, marker=m, markevery=100)
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

def find_param(rs = [4., 8.], ks=list(range(1, 16)), m=3):
    x = np.arange(0.001, 10, step=0.01)
    y = gaussian(x)

    cur_config = [(0, 0) for i in range(m)]
    print(cur_config)

    best_plan = [None]
    def dfs(cur_lvl, min_kidx, best_plan):
        if cur_lvl>=m:
            # do the actual computation
            fss = np.array([((collision_prob_l2(x, rs[ridx])**ks[kidx])) for (ridx, kidx) in cur_config] )
            # fss /= np.sum(fss*(1-fss), axis = 1).reshape((-1, 1))**0.5
            # print('fss=', fss)

            alpha = 0.001
            reg_res = Ridge(fit_intercept=False, alpha=alpha).fit(fss.T, y)
            ws = reg_res.coef_
            cp = np.dot(ws.T, fss)
            error = np.sum((cp-y)**2) + np.sum(ws**2) * alpha 
            bias = np.sum((cp-y)**2)
            if best_plan[0] is None or best_plan[0][0] >= error:
                print(cur_config, ws, error, bias)
                best_plan[0] = (error, cur_config.copy(), ws)
            return 
        for ri in range(0, len(rs)):
            for ki in range(min_kidx, len(ks)):
                cur_config[cur_lvl] = (ri, ki)
                dfs(cur_lvl+1, ki, best_plan)
    dfs(0, 0, best_plan)

    return best_plan[0]

#remarks
#l1-laplacion: 5.008, 4, 2
#l2-gaussian

xs= np.arange(0.001, 10, 0.05)

def find_param_greedy(rs = [1,2,4], maxk=8, maxl=3, target=gaussian(xs)):
    x = np.arange(0, 10, step=0.01)
    y = gaussian(x)

    #generate fs

    ps = [collision_prob_l2(xs, r) for r in rs]
    pidxs = []
    def dfs(cur_lvl, lastidx, cur_pidx, pidxs):
        if cur_lvl > maxk:
            return 
        if cur_lvl >= 1:
            # print(cur_lvl, lastidx)
            # print(cur_pidx)
            pidxs += [cur_pidx]
        
        for i in range(lastidx, len(rs)):
            dfs(cur_lvl+1, i, cur_pidx+[i], pidxs)

    dfs(0, 0, [], pidxs)
    # print(fs)
    fs = np.ones(shape=(len(pidxs), len(xs)))
    for i, pidx in enumerate(pidxs):
        for pi in pidx:
            fs[i] = fs[i]*ps[pi]

    fs /= np.sum(fs*(1-fs), axis = 1).reshape((-1, 1))**0.5

    best_error = [1e9]
    best_res = None
    def dfs2(cur_lvl, lastidx, cur_fidxs):
        if cur_lvl >= maxl:
            #compute error
            fss = np.array([fs[i] for i in cur_fidxs])
            cur_pidxs = [pidxs[i] for i in cur_fidxs]
            # print(fss.shape)

            reg_res = Ridge(fit_intercept=False, alpha=1.).fit(fss.T, target)
            ws = reg_res.coef_
            cp = np.dot(ws.T, fss)
            error = np.sum((cp-target)**2) + np.sum(ws**2)

            if error < best_error[0]:
                best_error[0] = error
                best_res = (error, ws, cur_pidxs)
                print(error, ws, cur_pidxs)
            return 
        
        for i in range(lastidx, len(fs)):
            dfs2(cur_lvl+1, i, cur_fidxs+[i])

    dfs2(0, 0, [])

def plot_guassian_linear_combine(ws = np.array([-3.8232685518081553, -8.620702462685585, 13.441813962243032]), 
        ks=np.array([6, 7, 9]), rs=np.array([2.5, 4, 5]), c='r', m='o'):
    # [-3.8232685518081553, 2.501, 6, -8.620702462685585, 4.001, 7, 13.441813962243032, 5.001, 9, 0.0008193039957631349]

    x = np.arange(0, 10, step=0.001)
    # y = gaussian(x)

    y = np.zeros_like(x)

    for w, k, r in zip(ws, ks, rs):
        y += w*(collision_prob_l2(x, r) **k)
    plt.plot(x, y, label='params=(%s, %s, %s)'%(str(ws), str(ks), str(rs)), color=c, marker=m, markevery=100)

def plot_plan(plan, rs, ks, *args, **kwargs):
    print('plan0=', plan)
    x = np.arange(0.001, 10, step=0.01)

    fs = []
    error, rks, ws = plan
    for ridx, kidx in rks:
        r = rs[ridx]
        k = ks[kidx]
        fs += [collision_prob_l2(xs, r)**k]
    
    # print(fs)

    cp = np.dot(ws.T, fs)
    # print(ws, error)

    plt.plot(xs, cp, *args, **kwargs)

    return error, ws

if __name__ == '__main__':
    # find_param(rs = np.arange(0.001, 8, 0.1), ks=range(1, 10) )
    # find_param(rs=np.arange(1, 8, 0.5), maxk=4)

    plt.rcParams.update({'font.size': 17})
    # rs = [1, 2, 4, 8]
    rs = [4]
    ks = list(range(1, 16))

    # plan0 = (1.1202124643187696, [(0, 0), (0, 3), (0, 3), (0, 3), (0, 3), (0, 18), (0, 18)], np.array([-0.08297805,  0.39493347,  0.39493347,  0.39493347,  0.39493347,
    #    -0.27749366, -0.27749366]))   #k<=24
    # plan0= (1.1856431873768931, [(0, 0), (0, 3), (0, 3), (0, 3), (0, 3), (0, 14), (0, 14)], np.array([-0.09107084,  0.41073369,  0.41073369,  0.41073369,  0.41073369,
    #    -0.28826238, -0.28826238])) #maxk<=16
    plan0= (0.02347270587591245, [(0, 0), (0, 2), (0, 4), (0, 4), (0, 9), (0, 10), (0, 14)], np.array([ 0.03581134, -0.87357163,  1.97736006,  1.97736006, -1.65327229,
       -1.6066025 ,  1.1412126 ]))
    # plan0 = find_param(rs = rs, ks=ks, m=7)
    plot_guassian()
    plot_plan(plan0, rs = rs, ks=ks, color='b', label='linear_regression')
    plt.xlim(0, 5)
    plt.ylim(-0.05, 1.05)
    plt.legend()
    plt.savefig('kernel_linear.png', bbox_inches='tight')
    plt.show()

    # [0.35340151512900775, 5.001, 5, 0.6328616002703029, 5.001, 5, 0.23401174100234812, 4.501, 4, 1.6462740178400752]
    # [0.7242831100993499, 5.101000000000001, 5, 0.04772281072837359, 4.301, 4, 0.44470348690083483, 4.2010000000000005, 4, 1.640864830962128],

    # ws = [0.7242831100993499, 0.04772281072837359, 0.44470348690083483]
    # ks = [5, 4, 4]
    # rs = [5.1, 4.3, 4.2]
    # plot_guassian_linear_combine(c='b', m='o', ws= ws, ks=ks, rs=rs)
    # plot_guassian()

    # plt.show()