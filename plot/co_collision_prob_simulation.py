import numpy as np
import matplotlib.pylab as plt

nData = 10000
p = np.random.normal(size=(nData, 3))
p = p/np.linalg.norm(p, axis=-1).reshape((-1, 1))

xs = np.random.normal(size=(3, 3))
xs = xs/np.linalg.norm(xs, axis=-1).reshape((-1, 1))

projections = np.dot(p, xs.T)

def kind(projections):
    signed_proj = projections>0
    return np.sum(signed_proj, axis=-1)

pkind = kind(projections)

print(pkind)
# ptrue = p[ppred]
# pfalse = p[np.logical_not(ppred) ]

# from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(xs[:, 0], xs[:, 1], xs[:, 2], color='r', marker='*', s=200, alpha=1.)

colors = ['k', 'b', 'g', 'm']
for i in range(2, 4):
    pi = p[pkind==i]
    ax.scatter(pi[:, 0], pi[:, 1], pi[:, 2], color=colors[i], alpha=0.2*i+0.1)
plt.show()