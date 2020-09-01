import matplotlib.pyplot as plt
import numpy as np


fig = plt.figure()
x = np.loadtxt('error1.dat')
x_mean = np.zeros((x.shape[0],), dtype=float)
for i in range(5):
    j = i + 1
    x = np.loadtxt('error' + str(j) + '.dat')
    x_mean += x
plt.plot(x_mean/5, label='mean')
plt.legend()
plt.savefig('error.png')
