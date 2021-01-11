import matplotlib.pylab as plt
import numpy as np
import sympy as sp
import math
from scipy.special import erf
from itertools import product, cycle
from scipy.stats import norm
from sklearn.linear_model import LinearRegression, Ridge
from copy import copy, deepcopy

xs = np.arange(0.001, 10, 0.05)

def collision_prob_l2(d, r):
    return 1-2*norm.cdf(-r/d)-2/(2**0.5*np.pi**0.5*(r/d))*(1-np.exp(-(r/d)**2/2))

def gaussian(x):
    return np.exp(-x**2*0.5)


def gen_unit_prob(ds=xs, rs=(0.1, 8, 0.1)):
    height = len(rs)
    width = len(ds)

    unit_prob = np.zeros(shape=(height, width))

    for i, r in enumerate(rs):
        unit_prob[i] = collision_prob_l2(ds, r)

    return unit_prob

def gen_target(ds=xs):
    return gaussian(ds)


unit_probs = gen_unit_prob()
target = gen_target()

class Plan:
    __slots__ = ['ps', 'fs', 'ws', 'err', 'cost', 'ds', 'pvecs', 'fvecs']

    def __init__(self, ps=[], fs=[], ws=[], err=1e9, cost=0):
        self.ps = deepcopy(ps)
        self.fs = deepcopy(fs)
        self.ws = ws
        self.err = err
        self.cost = cost

    @classmethod
    def new_from(cls, another_plan):
        return Plan(ps=another_plan.ps, fs=another_plan.fs)

    def new_plan_by_append_new_p(self, p):
        new_plan = Plan(self.ps.copy(), self.fs.copy(), self.ws.copy(), self.err, self.cost)

        new_plan.ps += [p]
        if len(self.fs) > 0:
            new_plan.fs[-1] += [len(new_plan.ps)-1 ]
        else:
            new_plan.fs = [[len(new_plan.ps)-1 ]]
        return new_plan



    def get_ws(self, ds=xs, target= target):
        #suppose ps, fs aer given, compute ws
        self.pvecs = np.zeros(shape=(len(self.ps), len(ds)) )
        for i, pr in enumerate(self.ps):
            self.pvecs[i] = collision_prob_l2(ds, pr)

        self.fvecs = np.ones(shape=(len(self.fs), len(ds)))
        for i, fidxs in enumerate(self.fs):
            for idx in fidxs:
                self.fvecs[i] = self.fvecs[i] * self.pvecs[idx]

        # self.fvecs /= np.mean(self.fvecs*(1-self.fvecs), axis = 1)**0.5
        fvars = np.sum(self.fvecs*(1-self.fvecs), axis = 1)
        fvecs_normalized = self.fvecs / fvars.reshape((-1, 1))**0.5
        # self.fvecs /= fvars.reshape((-1, 1))**0.5
        # print(fvecs)

        reg_res = Ridge(fit_intercept=False).fit(fvecs_normalized.T, target)

        # error = reg_res.score(fvecs.T, target)
        ws = reg_res.coef_ 

        cp = np.dot(ws.T, fvecs_normalized)
        error = np.sum((cp-target)**2) + np.sum(ws**2)

        return error, ws/ fvars**0.5


def get_plan_l1(rs, target, max_m = 10):
    plans = [Plan() for i in range(max_m)]
    best_err_all = 1e9
    best_plan_all = None
    for i in range(1, max_m):
        best_err = 1e9
        best_plan = None
        for r in rs:
            new_plan = plans[i-1].new_plan_by_append_new_p(r)
            new_err, ws = new_plan.get_ws()
            print(new_plan.ps, new_plan.fs, new_err, ws)
            if new_err < best_err:
                best_err = new_err
                best_plan = Plan.new_from(new_plan)
        plans[i] = Plan.new_from(best_plan)

        if best_err < best_err_all:
            best_err_all = best_err
            best_plan_all = Plan.new_from(best_plan)
        else:
            break
    return plans

# def get_plans(unit_probs, target, l_threshold, m_threshold):

# print(unit_probs)

# plans = get_plan_l1(np.arange(1., 10, 0.1), target) 

def plot_plan(plan, *args, **kwargs):
    error, ws = plan.get_ws()

    cp = np.dot(ws.T, plan.fvecs)
    print(ws, error)

    plt.plot(xs, cp, *args, **kwargs)

    return error, ws


# plan_to_plot = Plan(ps=[1, 2, 3, 4, 5] ,fs=[[0, 0], [0, 1, 1], [1, 1, 2, 3, 3, 3]])
# plan_to_plot = Plan(ps=[1.7000000000000007, 2.5000000000000014, 4.100000000000003] ,fs=[[0, 1, 2]])
plan_to_plot = Plan(ps=np.arange(1, 8, 0.5), fs=[[0, 0, 0, 0], [7, 7, 1, 7], [7, 7, 7, 5], [7, 7, 6, 7], [7, 7, 7, 7], [7, 7, 7, 7], [7, 7, 7, 7], [7, 7, 7, 7]] )
plot_plan(plan_to_plot, 'b-')
plt.plot(xs, target, 'r-')
plt.show()