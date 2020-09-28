import numpy as np

x = np.loadtxt('gendata.dat')
print(x.shape)

y = np.loadtxt('randomdata.dat')
print(y.shape)

z = x + y
print(z.shape)
print(x[0])
print(y[0])
print(z[0])
np.savetxt('datawithnoise.dat', z)

#test
i = np.loadtxt('datawithnoise.dat')
print(i.shape)

