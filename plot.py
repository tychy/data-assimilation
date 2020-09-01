from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import sys

args = sys.argv
x = np.loadtxt('rk-F' + args[1] + '.dat')
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x[:, 0], x[:, 1], x[:, 2])
ax.set_xlabel('$x_1$')
ax.set_ylabel('$x_2$')
ax.set_zlabel('$x_3$')
plt.savefig('rk-F' + args[1] + '.png')

#plt.plot(x[:, 0])
#plt.savefig('rk-F' + args[1] + 'x0.png')

