import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

p = np.loadtxt("trayectoria_E_alpha.dat")

figure = figure()
ax = Axes3D(figure)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.plot(p[:,0], p[:,1], p[:, 2])
show()
