from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import sys

args = sys.argv
x = np.loadtxt('norm.dat')
fig = plt.figure()
plt.hist(x,bins=100)
plt.savefig('norm-len{}mean{:.5f}var{:.5f}.png'.format(len(x),np.mean(x),np.var(x)))

