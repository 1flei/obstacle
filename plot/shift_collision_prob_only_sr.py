import matplotlib.pylab as plt
import numpy as np
import sympy as sp
import math
from scipy.special import erf, erfi, gammainc
from itertools import product, cycle
from scipy.stats import norm
from math import factorial
from sklearn.linear_model import LinearRegression, Ridge

#plot the collision probability and the gaussian pdf


def gaussian(x):
    return np.exp(-x**2*0.5)

def laplacian(x):
    return np.exp(-np.abs(x))


def shift_collsion_prob_l2(d, w, s):
    sqrt2over2 = 2**0.5 / 2
    pi = np.pi
    exp = np.exp
    integral0= lambda t:0.5*(1-s/w)*erf(sqrt2over2*t/d) - sqrt2over2*pi**(-0.5)*d/w*exp(-0.5*t**2/d**2)
    integral1= lambda t:0.5*(1+s/w)*erf(sqrt2over2*t/d) + sqrt2over2*pi**(-0.5)*d/w*exp(-0.5*t**2/d**2)

    return integral0(s) - integral0(s-w) + integral1(s+w) - integral1(s)



def fit_gaussian_using_sr(maxk=8, xs=np.arange(0.001, 5, 0.05)):
    target = gaussian(xs)
    print('gen fs')
    fvecs = np.array([shift_collsion_prob_l2(xs, w, s*w)**k for w in np.arange(0.1, 10, 0.02) for s in np.arange(0, 1, 0.02) for k in range(1, maxk)])

    print(fvecs.shape)

    fvars = np.sum(fvecs*(1-fvecs), axis = 1)
    fvecs_normalized = fvecs / fvars.reshape((-1, 1))**0.5
    # self.fvecs /= fvars.reshape((-1, 1))**0.5
    # print(fvecs)


    print('now fitting')
    reg_res = Ridge(fit_intercept=False).fit(fvecs_normalized.T, target)

    # error = reg_res.score(fvecs.T, target)
    ws = reg_res.coef_ 
    cp = np.dot(ws.T, fvecs_normalized)

    error = np.sum((cp-target)**2) + np.sum(ws**2)

    # ws_ = ws.copy()
    # ws_[ws_<0] = 0.
    # cp_ = np.dot(ws_.T, fvecs_normalized)

    plt.plot(xs, cp, 'b-')
    plt.plot(xs, target, 'r-')
    plt.show()

    return error, ws/ fvars**0.5
    
if __name__ == '__main__':
    error, w = fit_gaussian_using_sr()
    print(error, list(w))