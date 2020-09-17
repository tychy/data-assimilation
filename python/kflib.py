from tqdm.notebook import tqdm_notebook as tqdm
import numpy as np
import random
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = (20 ,6)
from scipy.integrate import odeint
from sklearn.metrics import mean_squared_error
from params import *

def l96_model(x, t):
    """Lorenz 96 model with constant forcing"""
    # Setting up vector
    d = np.zeros(N)
    # Loops over indices (with operations and Python underflow indexing handling edge cases)
    for i in range(N):
        d[i] = (x[(i + 1) % N] - x[i - 2]) * x[i - 1] - x[i] + F
    return d
        
def l96_step(x, step):
    dt = 0.01
    hdt = dt /2
    t = np.arange(0.0, step, dt)
    x_ls = []
    for _ in t:
        k1 = l96_model(x, t)
        k2 = l96_model(x+k1 * hdt, t)
        k3 = l96_model(x+k2 *hdt, t)
        k4 = l96_model(x+k3 * dt, t)
        x = x + (k1+2*k2 + 2*k3 + k4)*dt/6
        x_ls.append(x)
    return np.array(x_ls)


def kf_forecast(xa, pa,R,H, param=1.0):
    #pa = pa * 1.1
    step = 0.01
    xf = l96_step(xa, 0.01)[-1]
    def getM(xc):
        m = []
        delta = 0.001
        for i in range(xc.shape[0]):
            b = np.zeros_like(xc)
            b[i] += delta
            m.append(l96_step(xc + b, step)[-1] - l96_step(xc, step)[-1])
        return np.array(m).T /delta
    M = getM(xa)
    pf = M @ pa @ M.T
    return xf, pf * param


def henbun_forecast(xa, pa):
    xf = l96_step(xa, 0.01)[-1]
    return xf, pa

def kf(x_last, data, datawithnoise,step=100, param=1.0, del_key=[]):
    # 初期値
    pa = np.eye(N) * 25
    rs = 0.05 * int(np.random.randint(1, 1000))
    xa = l96_step(x_last, rs)[-1]
    I = np.eye(N)
    Z = np.zeros_like(I)
    H_del = np.eye(N)
    H_del = np.delete(H_del, del_key, 0)
    
    xferror_before_assim =[]
    xferror_after_assim = []
    trpa = []
    yrmse = []
    for  i in tqdm(range(step * 5)):
        if i % 5 == 0:
            H = H_del
        else:
            H = Z
            xa = l96_step(xa, 0.01)[-1]
            continue
        R = np.eye(H.shape[0])
        yt = data[i//5]
        y = np.array(datawithnoise[i//5])
        y = np.delete(y, del_key)
        
        xf, pf = kf_forecast(xa, pa,R,H, param)
        xferror_before_assim.append( np.sqrt(mean_squared_error(yt,  xf)))
        K = pf @ H.T @ np.linalg.inv(H@ pf @ H.T + R)
        xa = xf + K @ (y - H @ xf)
        pa = (np.eye(N)- K @ H )@ pf
        xferror_after_assim.append( np.sqrt(mean_squared_error(yt,  xa)))
        trpa.append(np.sqrt(np.trace(pa)/N))
    return xferror_before_assim, xferror_after_assim, trpa


def henbun_forecast(xa, pa):
    xf = l96_step(xa, 0.01)[-1]
    return xf, pa


def henbun(x_last, data, datawithnoise,step=100, del_key=[]):
    # 初期値
    pa = np.eye(N)
    rs = 0.05 * int(np.random.randint(1, 1000))
    xa = l96_step(x_last, rs)[-1]
    H = np.eye(N)
    Z = np.zeros((N,N))
    
    xferror_before_assim =[]
    xferror_after_assim = []
    trpa = []
    yrmse = []
    H_del = np.eye(N)
    H_del = np.delete(H_del, del_key, 0)
    for  i in tqdm(range(step * 5)):
        if i % 5 == 0:
            H = H_del
        else:
            H = Z
            xa = l96_step(xa, 0.01)[-1]
            continue
        R = np.eye(H.shape[0])
        yt = data[i//5]
        y = np.array(datawithnoise[i//5])
        y = np.delete(y, del_key)
        
        xf, pf = henbun_forecast(xa, pa)
        xferror_before_assim.append( np.sqrt(mean_squared_error(yt,  xf)))
        K = pf @ H.T @ np.linalg.inv(H@ pf @ H.T + R)
        xa = xf + K @ (y - H @ xf)
        xferror_after_assim.append( np.sqrt(mean_squared_error(yt,  xa)))
        trpa.append(np.sqrt(np.trace(pa)/N))
    return xferror_before_assim, xferror_after_assim, trpa 