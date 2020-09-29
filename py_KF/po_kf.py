import numpy as np
import random
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (20, 6)
from scipy.integrate import odeint
from sklearn.metrics import mean_squared_error
from kflib import *
from params import *


def po_forecast(xa):
    return l96_step(xa, 0.01)[-1]


def po(x_last, data, datawithnoise, step=100, param=1.0, del_num=0, m=8):
    I = np.eye(N)
    xf_ls = []
    xa_ls = np.array([
        l96_step(x_last, 0.05 * int(np.random.randint(1, 1000)))[-1] for _ in range(m)
    ])  # 初期値を作成

    Z = np.zeros_like(I)
    H_del = np.eye(N)
    del_key = np.random.choice(N, del_num, replace=False)
    H_del = np.delete(H_del, del_key, 0)
    Z_del = np.delete(Z, del_key, 0)

    xferror_before_assim = []
    xferror_after_assim = []
    trpa = []
    yrmse = []
    isappend = False
    for i in range(step * 5):
        if i % 5 == 0:
            H = H_del
            R = np.eye(H.shape[0])
            yt = data[i // 5]
            y = np.array(datawithnoise[i // 5])
            y = np.delete(y, del_key)

            xf_mean = np.zeros(N)
            xf_ls = []
            for k in range(m):
                xa = xa_ls[k]
                xf = po_forecast(xa)
                xf_ls.append(xf)
                xf_mean = xf_mean + xf
            xf_mean = xf_mean / m
            xf_ls = np.array(xf)
            xferror_before_assim.append(np.sqrt(mean_squared_error(yt, xf)))
            dxf = xf_ls - xf_mean * np.ones(N)
            dyf = H @ xf - H @ xf_mean
            K = dxf @dyf.T np.linalg.inv(dyf @dxf.T+(m-1)*I)
            xa_ls = xf + K @ (y-H@xf)
            xferror_after_assim.append(np.sqrt(mean_squared_error(yt, xa)))
        else:
            H = Z_del
            for k in range(m):
                xa = xa_ls[k]
                xf = po_forecast(xa)
                xa_ls[k] = xf
    return xferror_before_assim, xferror_after_assim