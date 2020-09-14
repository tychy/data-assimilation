import matplotlib.pyplot as plt
import numpy as np


def main():
    fig = plt.figure()
    x = np.loadtxt("3derror.dat")
    plt.plot(x, label='3dhenbun')
    plt.grid(which='major',color='black',linestyle='-')
    plt.xlabel("hour/6")
    plt.ylabel("RMSE(x_f - x_t)")

    plt.legend()
    plt.savefig('error.png')


if __name__ == "__main__":
    main()
