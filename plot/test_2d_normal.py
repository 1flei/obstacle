import numpy as np
import matplotlib.pylab as plt



dim = 1000
theta = np.random.uniform(0, 2*np.pi, size=dim)
r = np.random.normal(size=dim)

x = r*np.cos(theta)

plt.hist(x, bins=30)
plt.hist(r, bins=30)
plt.show()