import numpy as np
import matplotlib.pylab as plt
import json
import random
from math import sqrt
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from itertools import chain
import re
import os

from scipy.spatial import Delaunay
import numpy as np

#HARITA BERLIAN
gt_vessel0 = np.array([
            [1.286, 104.183], 
            [1.2807, 104.18317], 
            [1.28, 104.17683], 
            [1.2855, 104.17617], 
            [1.286, 104.183]
        ])[:, ::-1]

#THORCO CLOUD
gt_vessel1 = np.array([[1.204300, 103.898650], 
            [1.208133, 103.901667], 
            [1.217517, 103.917367], 
            [1.213350, 103.917333], 
            [1.204500, 103.902600], 
            [1.204300, 103.898650]])[:, ::-1]

#Cai Jun 3
gt_vessel2 = np.array([
            [1.431917, 104.456483], 
            [1.437067, 104.459750], 
            [1.426917, 104.462750], 
            [1.431617, 104.456517], 
            [1.431917, 104.456483]
        ])[:, ::-1]


def pntdist2(ai, bi):
    return np.sum((ai[1:]-bi[1:])**2)

def dtw_flat(a, b, l=6, c=2):
    a = a.reshape((-1, 3))
    b = b.reshape((-1, 3))
    return dtw(a, b, l, c)

def traj_length(a):
    length = 1e-9
    for ai, aip1 in zip(a, a[1:]):
        length += pntdist2(ai, aip1)**0.5
    return length

def dtw_normalized_flat(a, b, l=6, c=2):
    a = a.reshape((-1, 3))
    b = b.reshape((-1, 3))
    dtw2 = dtw(a, b, l, c)
    lena = traj_length(a)
    lenb = traj_length(b)

    # print('dtw2={}, lena={}, lenb={}'.format(dtw2, lena, lenb) )

    return dtw2/lena/lenb

def get_z_score(dq, dr, dq_comover, dr_comover):
    dq_ncom = dq-dq_comover
    dr_ncom = dr-dr_comover

    if dr == 0:
        dr = 1e-9
    if dq == 0:
        dq = 1e-9

    p1_hat = dr_comover/dr
    p2_hat = dq_comover/dq

    p_hat = (dr*p1_hat + dq*p2_hat) / (dr+dq)
    sigmad = sqrt(p_hat*(1-p_hat)*(1/dr + 1/dq))
    if sigmad==0:
        return 0

    z = (p1_hat-p2_hat) / sigmad

    return z

def dtw(a, b, l=6, c=2):
    assert(len(a)==l and len(b)==l)

    dtwij = np.zeros(shape=(l+1, l+1))
    dtwij[:, :] = 1e9
    dtwij[0, 0] = 0
    for i in range(1, l+1):
        for j in range( max(1, i-c), min(l, i+c)+1 ):
            cost = pntdist2(a[i-1], b[j-1])
            dtwij[i, j] = cost + min(dtwij[i-1, j], dtwij[i, j-1], dtwij[i-1, j-1])
    # print('dtwij=', dtwij)

    # for ai in a:
    #     for bj in b:
    #         d = pntdist2(ai, bj)
    #         print('{}, '.format(d), end='')
    #     print()
    return dtwij[l, l]


def alpha_shape(points, alpha, only_outer=True):
    """
    Compute the alpha shape (concave hull) of a set of points.
    :param points: np.array of shape (n,2) points.
    :param alpha: alpha value.
    :param only_outer: boolean value to specify if we keep only the outer border
    or also inner edges.
    :return: set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
    the indices in the points array.
    """
    assert points.shape[0] > 3, "Need at least four points"

    def add_edge(edges, i, j):
        """
        Add an edge between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            assert (j, i) in edges, "Can't go twice over same directed edge right?"
            if only_outer:
                # if both neighboring triangles are in shape, it's not a boundary edge
                edges.remove((j, i))
            return
        edges.add((i, j))

    tri = Delaunay(points)
    edges = set()
    # Loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.vertices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Computing radius of triangle circumcircle
        # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)
    return edges


def plot_data(data, z_threshold=None, is_plot_dist=False, is_plot_density=True, dr_thres=None, dq_thres=None, 
        is_plot_q=True, qalpha=1., neighboralpha=0.2, is_plot_obst_traj=True):
    q = np.array(data['q'])
    t = np.array(data['t'])
    nqts = np.array(data['Nqt']) if 'Nqt' in data else np.array([])
    nrts = np.array(data['Nrt']) if 'Nrt' in data else np.array([])
    # print(nqts, nrts)

    q_succ = np.array(data['q_succ'])
    t_succ = np.array(data['t_succ'])
    nqt_succs = np.array(data['Nqt_succ']) if 'Nqt_succ' in data else np.array([])
    nrt_succs = np.array(data['Nrt_succ']) if 'Nrt_succ' in data else np.array([])

    densities = data['densities']

    if z_threshold is not None:
        dq, dr, dq_comover, dr_comover = tuple(densities)
        z_score = get_z_score(dq, dr, dq_comover, dr_comover)
        if z_score < z_threshold:
            return False
        print('z=', z_score)

    if dq_thres is not None:
        dq = densities[0]
        if dq < dq_thres:
            return False
    if dr_thres is not None:
        dr = densities[1]
        if dr < dr_thres:
            return False


    # 103.70, 104.10, 1.1, 1.3
    # if t[-2] < 103.70 or t[-2] > 104.1 or t[-1] < 1.1 or t[-1] > 1.3:
    #     return False

    for xi in nqts:
        plt.plot(xi[1::3], xi[2::3], 'b.-', alpha=neighboralpha)
        d = dtw_normalized_flat(t, xi)
        if is_plot_dist:
            plt.text(xi[-2], xi[-1], '%f'%d**0.5, alpha=0.2, color='b')
    for xi in nrts:
        plt.plot(xi[1::3], xi[2::3], 'g.-', alpha=neighboralpha)
        d = dtw_normalized_flat(t, xi)
        if is_plot_dist:
            plt.text(xi[-2], xi[-1], '%f'%d**0.5, alpha=0.2, color='b')

    for xi in nqt_succs:
        plt.plot(xi[-5::3], xi[-4::3], 'b-', alpha=neighboralpha)
        plt.plot(xi[-2], xi[-1], 'bo', alpha=neighboralpha)
        d = dtw_normalized_flat(t_succ, xi)
        if is_plot_dist:
            plt.text(xi[-2], xi[-1], '%f'%d**0.5, alpha=0.2, color='b')
    for xi in nrt_succs:
        plt.plot(xi[-5::3], xi[-4::3], 'g-', alpha=neighboralpha)
        plt.plot(xi[-2], xi[-1], 'go', alpha=neighboralpha)
        d = dtw_normalized_flat(t_succ, xi)
        if is_plot_dist:
            plt.text(xi[-2], xi[-1], '%f'%d**0.5, alpha=0.2, color='b')

    if is_plot_q:
        plt.plot(q[1::3], q[2::3], 'rx-', zorder=3, alpha=qalpha)
        plt.plot(q_succ[-5::3], q_succ[-4::3], 'rx-', zorder=3, alpha=qalpha)

    if is_plot_obst_traj:
        plt.plot(t[1::3], t[2::3], 'r.-', zorder=3, alpha=qalpha)
    plt.plot(t_succ[-5::3], t_succ[-4::3], 'r-', zorder=3, alpha=qalpha)
    plt.plot(t_succ[-2], t_succ[-1], 'ro', zorder=3, alpha=qalpha)

    if is_plot_density:
        plt.text(t[-2], t[-1], '%s'%densities, alpha=1.)
    return True

    
def plot_data_dist(data):
    q = np.array(data['q'])
    t = np.array(data['t'])
    nqts = np.array(data['Nqt']) if 'Nqt' in data else np.array([])
    nrts = np.array(data['Nrt']) if 'Nrt' in data else np.array([])
    # print(nqts, nrts)

    q_succ = np.array(data['q_succ'])
    t_succ = np.array(data['t_succ'])
    nqt_succs = np.array(data['Nqt_succ']) if 'Nqt_succ' in data else np.array([])
    nrt_succs = np.array(data['Nrt_succ']) if 'Nrt_succ' in data else np.array([])


    dists_nqts = np.array([dtw_normalized_flat(t, xi) for xi in nqts])
    dists_nrts = np.array([dtw_normalized_flat(t, xi) for xi in nrts])
    dists_nqtsucc = np.array([dtw_normalized_flat(t_succ, xi) for xi in nqt_succs])
    dists_nrtsucc = np.array([dtw_normalized_flat(t_succ, xi) for xi in nrt_succs])

    densities = data['densities']
    print(densities)
    plt.plot(dists_nqts**0.5, 'b-')
    plt.plot(dists_nqtsucc**0.5, 'b--')
    plt.plot(dists_nrts**0.5, 'r-')
    plt.plot(dists_nrtsucc**0.5, 'r--')
    return True


def plot_data_and_distfig(data, gt=gt_vessel1, **kwargs):
    # print(i)
    ax = plt.subplot(1, 2, 1)
    if not plot_data(data, **kwargs):
        return False
        
    fig = plt.gcf()
    fig.set_size_inches(12, 5)
    plot_chosen_area(gt)

    ax = plt.subplot(1, 2, 2)
    plot_data_dist(data)
    return True

def plot_data_and_distfig_sgtaxi(data, z_threshold=None, is_plot_dist=False, is_plot_density=True):
    plt.figure(figsize=(12, 5))
    # print(i)
    ax = plt.subplot(1, 2, 1)
    plotERP(cm='r*')
    plot_data(data, z_threshold, is_plot_dist, is_plot_density)

    ax = plt.subplot(1, 2, 2)
    plot_data_dist(data)

def plot_data_obs_pnts(datas, is_plot_density=True, is_plot_succ=True, color='r', alpha=0.1, ls='--', tail_color=None):
    obs_pnts = np.array([data['obs_pnts'] for data in datas])

    try:
        hull_edges = ConvexHull(obs_pnts[:, 1:]).simplices
        # hull_edges = alpha_shape(obs_pnts[:, 1:], alpha=0.25)
        for simplex in hull_edges:
            plt.plot(obs_pnts[simplex, 1], obs_pnts[simplex, 2], color = color, ls=ls, marker='.', lw=2)
        # plt.plot(obs_pnts[hull.vertices,1], obs_pnts[hull.vertices,2], color = color, ls='--', marker='.', lw=2)
    except:
        # plt.plot(obs_pnts[:, 1], obs_pnts[:, 2], color = color, ls='', marker='.')
        pass

    # plt.plot(obs_pnts[:, 1], obs_pnts[:, 2], color = color, ls='', marker='.')

    for data in datas:
        obs_pnt = data['obs_pnts']

        if is_plot_succ:
            t_succ = np.array(data['t_succ'])
            plt.plot(t_succ[-8::3], t_succ[-7::3], color=tail_color if tail_color is not None else color, ls='-', alpha=alpha)
    if is_plot_density:
        densities = datas[0]['densities']
        plt.text(obs_pnt[-2], obs_pnt[-1], '%s'%densities, alpha=alpha)
    
    return True

# def get_gaussian_filtered_2dheat(pnts, sigma = 0.01):
#     minx = np.min(pnts[:, 1])
#     maxx = np.max(pnts[:, 1])
#     miny = np.min(pnts[:, 2])
#     maxy = np.max(pnts[:, 2])

#     nbins_x = int(np.ceil((maxx-minx)/sigma ) )
#     nbins_y = int(np.ceil((maxy-miny)/sigma ) )

#     hist = np.zeros(shape=(nbins_x, nbins_y))
#     for pnt in pnts:
#         xx = int(np.floor((pnt[1] - minx)/sigma ) )
#         yy = int(np.floor((pnt[2] - miny)/sigma ) )
#         for sigmax in range(-5, 6):
#             for sigmay in range(-5, 6):
#                 pnt_sigma = pnt + np.array([0, sigmax*sigma, sigmay*sigma])

#                 dist2 = np.sum((pnt-pnt_sigma)**2) 

#                 hist[xx, yy] += np.exp(-dist2/2/sigma/sigma)

#     return hist

def plot_data_obs_pnts_heatmap(datass, is_plot_succ=True, color='r', sigma=0.01, vmax=50):
    obs_pnts = np.array([data['obs_pnts'] for datas in datass for data in datas])

    from scipy.stats.kde import gaussian_kde

    x = obs_pnts[:, 1]
    y = obs_pnts[:, 2]
    print(x.min(), x.max(), y.min(), y.max())
    k = gaussian_kde(obs_pnts[:, 1:].T)
    k.set_bandwidth(sigma*5)
    nbins_x = int(np.ceil((x.max()-x.min())/sigma ) )
    nbins_y = int(np.ceil((y.max()-y.min())/sigma ) )
    xi, yi = np.mgrid[x.min():x.max():(nbins_x)*1j, y.min():y.max():nbins_y*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap='OrRd', vmax=vmax)

    # hist = get_gaussian_filtered_2dheat(obs_pnts, sigma=sigma)
    # plt.hist2d(hist[:, 0], hist[:, 1], bins=hist.shape, density=False)
    # plt.imshow(hist, cmap='OrRd')

    for datas in datass:
        for data in datas:
            obs_pnt = data['obs_pnts']

            if is_plot_succ:
                t_succ = np.array(data['t_succ'])
                plt.plot(t_succ[-8::3], t_succ[-7::3], color=color, ls='-', alpha=0.1)

    # if is_plot_density:
    #     densities = datas[0]['densities']
    #     plt.text(obs_pnt[-2], obs_pnt[-1], '%s'%densities, alpha=0.2)
    
    return True



# def plot_chosen_area(data, **kwargs):
#     for p, pp in zip(data, data[1:]):
#         plt.plot([p[1], pp[1]], [p[0], pp[0]], ls='-', **kwargs)
def plot_chosen_area(data, **kwargs):
    for p, pp in zip(data, data[1:]):
        plt.plot([p[0], pp[0]], [p[1], pp[1]], ls='-', **kwargs)

def getERPAtTime(t='08:00:01', fname = 'erp2.json'):
    with open(fname, 'r') as f:
        erpDict = json.load(f)
        
    ret = []  
    for key, values in erpDict.items():
        for v0, v1 in zip(values['rates'], values['rates'][1:]):
            if t >= v0[0] and t<v1[0]:
                if float(v0[1])>0:            
                    ret += [[values['lng'], values['lat']]]
                break
    return np.array(ret)
            
def plotERP(t='08:00:01', fname = 'erp2.json', cm='b*', alpha=1., markersize=15, **kwargs):
    erps = getERPAtTime(t, fname)
    plt.plot(erps[:, 0], erps[:, 1], cm, markersize=markersize, alpha=alpha, label='ERP Gantries', **kwargs)
            
def plotERPArea(t='08:00:01', fname = 'erp3.json', cm='m*', cm_line = 'm-', alpha=1., markersize=15):
    with open(fname, 'r') as f:
        erpDict = json.load(f)
    workingERP = {}
    for key, values in erpDict.items():
        for v0, v1 in zip(values['rates'], values['rates'][1:]):
            if t >= v0[0] and t<v1[0]:
                if float(v0[1])>0:            
                    workingERP[key] = values
                break
    for key, values in workingERP.items():
        lat, lng = values['lat'], values['lng']
        plt.plot([lng], [lat], cm, markersize=markersize, alpha=alpha)
        plt.text(lng, lat, key)

        area = values['area']
        area = np.array(area + [area[0]])
        plt.plot(area[:, 1], area[:, 0], cm_line, alpha=alpha)

def readERPArea(t='08:00:01', fname = 'erp3.json'):
    with open(fname, 'r') as f:
        erpDict = json.load(f)
    ret = []
    for key, values in erpDict.items():
        for v0, v1 in zip(values['rates'], values['rates'][1:]):
            if t >= v0[0] and t<v1[0]:
                if float(v0[1])>0:            
                    erp_area = np.array(values['area'])
                    erp_area_ = np.zeros_like(erp_area)
                    erp_area_[:, 0] = erp_area[:, 1]
                    erp_area_[:, 1] = erp_area[:, 0]
                    ret += [erp_area_]
                break
    return ret

def segment_point_mindist(p, p1, p2):
    if np.sum((p1-p2)**2)==0.:
        r = p1
    else:
        u = np.sum((p-p1)*(p2-p1)) / np.sum((p1-p2)**2)
        u = np.minimum(np.maximum(u, 0.), 1.)

        r = p1 + u*(p2-p1)
    return np.sum((p-r)**2)

def dot_prod(a, b):
    return a[0]*b[1]-a[1]*b[0]

# Return true if line segments AB and CD intersect

def is_two_segment_intersect(a, b, c, d, eps=1e-9):
    return dot_prod(b-a, c-a)*dot_prod(b-a, d-a)< eps and \
            dot_prod(c-d, a-d)*dot_prod(c-d, b-d)< eps

def polygon_point_mindist(q, poly):
    poly_p = chain(poly[1:], poly[:1])
    intersect_cnt = 0
    for p, pp in zip(poly, poly_p):
        if is_two_segment_intersect(q, q+1e9, p, pp):
            intersect_cnt +=1
    if intersect_cnt%2==1:
        #not contain
        return 0.
    min_dist = 1e9
    poly_p = chain(poly[1:], poly[:1])
    for p, pp in zip(poly, poly_p):
        min_dist = np.minimum(min_dist, segment_point_mindist(q, p, pp) )
    return min_dist



def get_latest_json_filename(resname = 'vessel1', sigma=0.1, delta=2, tau=100, folder='res'):
    patstr = r'{}_\[.*_sigma=?{}_delta=?{}(.0)?_tau=?{}\]\.json'.format(resname, sigma, delta, tau)
    # print('patstr=', patstr)
    pat = re.compile(patstr)
    latestFile = None
    for file in os.listdir(folder):
        # print(file, patstr)
        if pat.match(file):
            if latestFile is None:
                latestFile = file
            elif latestFile < file:
                latestFile = file
    if latestFile is None:
        return None
    return os.path.join(folder, latestFile)


    