import numpy as np
import os
import random
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (20, 6)
from scipy.integrate import odeint
from sklearn.metrics import mean_squared_error
from kflib import *
from params import *
from po_kf import *

N = 40
F = 8
x0 = F * np.ones(N)  # Initial state (equilibrium)
x0[0] += 0.01  # Add small perturbation to the first variable
x = l96_step(x0, 10)
x_one = l96_step(x0, 73.0)  # spin up
x = l96_step(x_one[-1], 73.0)
x_last = x_one[-1]

if not (os.path.exists("datawithnoise.txt") and os.path.exists("gendata.txt")):
    data = x[::5]
    np.savetxt("gendata.txt", data)
    random_ls = np.random.randn(data.shape[0], x[::5].shape[1])
    datawithnoise = data + random_ls
    np.savetxt("datawithnoise.txt", datawithnoise)
else:
    data = np.loadtxt("gendata.txt")
    datawithnoise = np.loadtxt("datawithnoise.txt")

m = 40
inflation_ls = [1.0, 1.05, 1.1, 1.2, 1.3, 1.5]
fig = plt.figure()
for inflation in inflation_ls:
    xferror_before_assim, xferror_after_assim = po(
        x_last, data, datawithnoise, step=500, del_num=0, m=m, inflation=inflation
    )
    np.savetxt("po_m{}_inflation{}.txt".format(m, inflation), xferror_after_assim)

    plt.plot(xferror_after_assim, label="po_m{}_inflation{}".format(m, inflation))
plt.legend()
plt.savefig("po_run.png")