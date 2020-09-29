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


def po(x_last, data, datawithnoise, step=100, del_num=0, m=8, inflation=1.0):
    I = np.eye(N)
    xa_ls = np.array(
        [l96_step(x_last, 0.05 * int(np.random.randint(1, 1000)))[-1] for _ in range(m)]
    )  # 初期値を作成
    xf_ls = xa_ls

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
        yt = data[i // 5]
        y = np.array(datawithnoise[i // 5])
        y = np.delete(y, del_key)
        if i % 5 == 0:
            H = H_del
            R = np.eye(H.shape[0]) * (m - 1)
            xf_mean = np.zeros(N)
            # xf_meanの計算
            for k in range(m):
                xa = xa_ls[k]
                xf = po_forecast(xa)
                xf_ls[k] = xf
                xf_mean = xf_mean + xf
            xf_mean = xf_mean / m
            # delta_xf=xf-xf_meanの計算
            xf_ls = np.array(xf_ls)
            xferror_before_assim.append(np.sqrt(mean_squared_error(yt, xf)))
            dxf = (xf_ls - xf_mean * np.ones(N)).T
            dyf = np.array([H @ xf_ls[idx] - H @ xf_mean for idx in range(m)]).T
            dxf = inflation * dxf  # inflation
            dyf = inflation * dyf  # inflation

            K = dxf @ dyf.T @ np.linalg.inv(dyf @ dxf.T + R)
            before_score = 0.0
            after_score = 0.0
            for k in range(m):
                xf = xf_ls[k]
                before_score += np.sqrt(mean_squared_error(yt, xf))
                xa = xf + K @ (y + np.random.randn() - H @ xf)
                xa_ls[k] = xa
                after_score += np.sqrt(mean_squared_error(yt, xa))
            xferror_before_assim.append(before_score / m)
            xferror_after_assim.append(after_score / m)

        else:
            H = Z_del
            for k in range(m):
                xa = xa_ls[k]
                xf = po_forecast(xa)
                xa_ls[k] = xf

    return xferror_before_assim, xferror_after_assim
