from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

x = np.loadtxt('rk.dat')
print(type(x))
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x[:, 0], x[:, 1], x[:, 2])
ax.set_xlabel('$x_1$')
ax.set_ylabel('$x_2$')
ax.set_zlabel('$x_3$')
plt.savefig("rk.png")
