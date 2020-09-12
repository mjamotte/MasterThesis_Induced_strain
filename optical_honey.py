from mpl_toolkits import mplot3d

import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure()
ax = plt.axes(projection="3d")
x = np.linspace(-6, 6, 100)
y = np.linspace(-6, 6, 100)

X, Y = np.meshgrid(x, y)
Z = np.cos(X/2+np.sqrt(3)/2*Y)*np.sin(X/2-np.sqrt(3)/2*Y)*np.sin(X)

ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Z/10, rstride=1, cstride=1,cmap='bone', edgecolor='none',rasterized=True)
ax.axis('off')


plt.show()