import matplotlib.pyplot as plt
import numpy as np


def parse_dat(fname):
    x = []#each contain numpy list
    with open(fname) as f:
        line = f.readline()
        data = np.array([],dtype=float)
        while True:
            line = f.readline()
            if line.find("time") > -1:
                x.append(data)
                data = np.array([],dtype=float)
                continue
            if not line:
                x.append(data)
                break
            data = np.append(data, float(line))
    return x


def main():
    fig = plt.figure()
    x = parse_dat('error.dat')
    x_mean = np.zeros((x[0].shape[0],), dtype=float)
    for i in range(len(x)):
        x_mean += x[i]
    plt.plot(x_mean/len(x), label='mean')
    plt.grid(which='major',color='black',linestyle='-')
    plt.legend()
    plt.savefig('error.png')


if __name__ == "__main__":
    main()
