import _ucrdtw
import sys
import numpy as np
import matplotlib.pyplot as plt
import math

if __name__ == '__main__':
    data = np.cos(np.linspace(0.0, math.pi * 6, 600)) + (np.random.uniform(-0.5, 0.5, 600) / 10)
    query = np.sin(np.linspace(0.0, math.pi * 2, 200))

    plt.figure()
    plt.plot(data)
    plt.plot(query)
    loc, dist = _ucrdtw.ucrdtw(data, query, 0.05, True)
    query = np.concatenate((np.linspace(0.0, 0.0, loc), query))
    plt.plot(query)
    plt.show()
