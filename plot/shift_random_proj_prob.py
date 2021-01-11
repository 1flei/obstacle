import matplotlib.pylab as plt
import numpy as np
import sympy as sp
import math
from scipy.special import erf, erfi, gammainc
from itertools import product, cycle
from scipy.stats import norm
from math import factorial

#plot the collision probability and the gaussian pdf


def gaussian(x):
    return np.exp(-x**2*0.5)

def laplacian(x):
    return np.exp(-np.abs(x))


def phi(x):
    #cdf of normal distirbution
    #2*phi = 1+erf(x/sqrt2)
    return (1.0 + erf(x / np.sqrt(2.0))) / 2.0

# def collision_prob_l2(c, r):
#     x = r/c
#     return 1 - 2*phi(-x) - 2/((2*np.pi)**0.5*x) * (1-np.exp(-x**2/2))


def collision_prob_l2(d, r):
    return 1-2*norm.cdf(-r/d)-2/(2**0.5*np.pi**0.5*(r/d))*(1-np.exp(-(r/d)**2/2))


def shift_collsion_prob_l2(d, w, s):
    # t, delta, s = sp.symbols('t d, s')
    # p0 = (1-(s-t)/w)
    # p1 = (1-(t-s)/w)
    # t_ = t/d
    # f = sp.exp(t_**2 * -0.5)*(2*sp.pi)**-0.5
    # print('integral0=', sp.simplify(sp.integrate(p0*f/d, t)) )
    # print('integral1=', sp.simplify(sp.integrate(p1*f/d, t)) )

    # res1 = sp.integrate(p0*f, t, s-w, s)+sp.integrate(p1*f, t, s-w, s)

    sqrt2over2 = 2**0.5 / 2
    pi = np.pi
    exp = np.exp
    integral0= lambda t:0.5*(1-s/w)*erf(sqrt2over2*t/d) - sqrt2over2*pi**(-0.5)*d/w*exp(-0.5*t**2/d**2)
    integral1= lambda t:0.5*(1+s/w)*erf(sqrt2over2*t/d) + sqrt2over2*pi**(-0.5)*d/w*exp(-0.5*t**2/d**2)

    return integral0(s) - integral0(s-w) + integral1(s+w) - integral1(s)



def gaussian_sp_diff(d0):
    d = sp.symbols('d')
    f = sp.exp(-d**2*0.5)
    df = f.diff(d).limit(d, d0).simplify()
    print(df)

def collision_prob_l2_sp_diff(r, d0=0.):
    d = sp.symbols('d')
    f=  -sp.erf(-r/d /2**0.5)-2/(2**0.5*sp.pi**0.5*(r/d))*(1-sp.exp(-(r/d)**2/2))
    df = sp.N(sp.diff(f, d).limit(d, d0))
    print(df)

def shift_collsion_prob_l2_sp_diff(w, s, d0=0.):
    s0 = s
    d, t, s = sp.symbols('d, t, s')
    sqrt2over2 = sp.sqrt(2) / 2
    pi = sp.pi
    exp = sp.exp
    erf = sp.erf
    integral0= 0.5*(1-s/w)*erf(sqrt2over2*t/d) - sqrt2over2*pi**(-0.5)*d/w*exp(-0.5*t**2/d**2)
    integral1= 0.5*(1+s/w)*erf(sqrt2over2*t/d) + sqrt2over2*pi**(-0.5)*d/w*exp(-0.5*t**2/d**2)
    f = integral0.subs(t, s) - integral0.subs(t, s-w) + integral1.subs(t, s+w) - integral1.subs(t, s)
    print(f.simplify())
    f = f.subs(s, s0)
    df = sp.N(sp.diff(f, d).limit(d, d0))
    print(df)

def exp_collsion_prob(d, lmd):
    # t, d = sp.symbols('t, d')
    # ft = sp.exp(-t**2 /2)/sp.sqrt(2*sp.pi)
    # fd = sp.integrate(ft*(1-sp.exp(-lmd*t*d)), t).simplify()
    # print('fd=', fd)
    # sqrt, pi, erf, e = sp.sqrt, sp.pi, sp.erf, sp.E
    # fd_ = 2*pi*erf(sqrt(2)*t/2)/2 -pi*sp.exp(d^2/2)*erf(t/sqrt(2)+d/sqrt(2))
    e, sqrt = np.e, np.sqrt
    a= lmd
    fd_ = -(e**((a**2*d**2)/2)*(erf((a*d)/sqrt(2))-1))
    # sp.plot(fd_, ylim=(0, 1))
    return fd_

def plot_incomplete_gamma(n=2):
    xs = np.arange(0.001, 5, 0.01)

    fx = 2**(n*0.5-1) * gammainc((n+1)*0.5, xs**2/2)
    plt.plot(xs, fx)
    plt.show()

def coeff_incoplete_gamma(n=2):
    if n==0:
        return 1
    if n%2==1:
        return 0.
    n=n/2
    return (-1)**n* np.sqrt(np.pi) / 2**n / n /factorial(n) / 2

    
if __name__ == '__main__':
    # plot_incomplete_gamma(8)
    for i in range(0, 20):
        if i%2==0:
            print(r'{}*x^{}'.format(coeff_incoplete_gamma(i), i) )
    # d0 = 0.005
    # gaussian_sp_diff(d0=d0)
    # collision_prob_l2_sp_diff(r=1., d0=d0)
    # shift_collsion_prob_l2_sp_diff(w=1., s=0.02, d0=d0)

    # print('d=1, w=1')
    # shift_collsion_prob_l2(d=1, w=1, s0=0.1)
    # print('d=2, w=1')
    # shift_collsion_prob_l2(d=2, w=1, s0=0.1)
    # print('d=1, w=2')
    # shift_collsion_prob_l2(d=1, w=2, s0=0.1)


    s = 0.
    xs = np.arange(0.001, 5, 0.01)
    f = shift_collsion_prob_l2(xs, w=1, s=s)

    expcp = exp_collsion_prob(xs, 2.)
    plt.plot(xs, gaussian(xs), 'r-')
    # plt.plot(xs, expcp, 'g-')
    from plot_collision_prob import plog_collision_prob_l2_kl
    # plog_collision_prob_l2_kl(np.pi**0.5)
    plt.plot(xs, f, 'b-')
    plt.ylim(0, 1)
    plt.show()