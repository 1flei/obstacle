import numpy as np
import matplotlib.pylab as plt
from plot_collision_prob import plog_collision_prob_l2_kl, plot_guassian
from shift_random_proj_prob import shift_collsion_prob_l2

n = 100
d = 2

xs = np.random.normal(size=(n, d)) / d**0.5

q = np.random.normal(size=(1, d)) / d**0.5

dists2 = np.sum((xs-q)**2, axis=-1)

kd = np.exp(-0.5 * dists2)


def get_hash(xs, a, b, r, l, k):
    p = np.ceil((np.dot(xs, a) + b) / r).astype(np.int32)

    for i in range(l):
        pi = p[:, i*k: (i+1)*k]
        pih = [hash(tuple(pij)) for pij in pi]
        yield pih

def get_hash_q(xs, a, b, r, l, k, s):
    p = np.ceil((np.dot(xs, a) + s + b) / r).astype(np.int32)

    for i in range(l):
        pi = p[:, i*k: (i+1)*k]
        pih = [hash(tuple(pij)) for pij in pi]
        yield pih

def get_cnts(xs, q, ntrials = 128, r=5.12, k=8, l=3):
    cnts = np.zeros(shape=(n))
    ass = np.random.normal(size=(d, k*l*ntrials)) 
    bss = np.random.uniform(0, r, size=(1, k*l*ntrials))
    for t in range(ntrials):
        # a = np.random.normal(size=(d, k*l))
        # b = np.random.uniform(0, r, size=(1, k*l))
        a = ass[:, t*(k*l): (t+1)*k*l]
        b = bss[:, t*(k*l): (t+1)*k*l]

        hxs = list(get_hash(xs, a, b, r, l, k) )
        hqs = list(get_hash(q,  a, b, r, l, k) )

        hxqs = np.array([[hi==hq[0] for hi in hx] for hx, hq in zip(hxs, hqs)])

        hxq = np.logical_or.reduce(hxqs, axis=0)
        # res = [h0 or h1 or h2 for (h0, h1, h2) in zip(hxq0, hxq1, hxq2)]
        cnts += hxq

    # print(cnts)
    return cnts /ntrials

def get_cnts_async(xs, q, ntrials = 128, r=5.12, s=0.01, k=8, l=3):
    cnts = np.zeros(shape=(n))
    ass = np.random.normal(size=(d, k*l*ntrials)) 
    bss = np.random.uniform(0, r, size=(1, k*l*ntrials))
    for t in range(ntrials):
        # a = np.random.normal(size=(d, k*l))
        # b = np.random.uniform(0, r, size=(1, k*l))
        a = ass[:, t*(k*l): (t+1)*k*l]
        b = bss[:, t*(k*l): (t+1)*k*l]

        hxs = list(get_hash(xs, a, b, r, l, k) )
        hqs = list(get_hash_q(q,  a, b, r, l, k, s) )

        hxqs = np.array([[hi==hq[0] for hi in hx] for hx, hq in zip(hxs, hqs)])

        hxq = np.logical_or.reduce(hxqs, axis=0)
        # res = [h0 or h1 or h2 for (h0, h1, h2) in zip(hxq0, hxq1, hxq2)]
        cnts += hxq

    # print(cnts)
    return cnts /ntrials

# for r in [5.12]:
#     for k in [8]:
#         for l in [3]:
#             print('r=', r, 'k=', k, 'l=', l)
#             cnts = get_cnts(xs, q, r=r, k=k, l=l)
#             plog_collision_prob_l2_kl(r=r, k=k, l=l, c='b')

#             plt.plot(dists2**0.5, cnts, 'g.')
#             plot_guassian()
#             plt.show()


def plog_shift_collision_prob_l2_kl(r, s= 0., c='b', m='.'):
    x = np.arange(0, 5, step=0.001)
    y = shift_collsion_prob_l2(x, r, s)
    y = 1-(1-(y**k))**l
    plt.plot(x, y, label='r=%.2f_k=%d_l=%d'%(r, k, l), color=c, marker=m, markevery=100)
    plt.xlim(0, 5)
    plt.ylim(0, 1)
    # plt.show()

for r in [1]:
    for k in [1]:
        for l in [1]:
            print('r=', r, 'k=', k, 'l=', l)
            s0 = 2
            cnts = get_cnts_async(xs, q, r=r, k=k, l=l, s= s0)
            plog_shift_collision_prob_l2_kl(r=r, c='b', s=s0)

            plt.plot(dists2**0.5, cnts, 'g.')
            plot_guassian()
            plt.show()
