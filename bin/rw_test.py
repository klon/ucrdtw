import _ucrdtw
import sys
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    data = np.cumsum(np.random.uniform(-0.5, 0.5, 1000000))
    query = np.cumsum(np.random.uniform(-0.5, 0.5, 100))
    loc, dist = _ucrdtw.ucrdtw(data, query, 0.05, True)

    query = np.concatenate((np.linspace(0.0, 0.0, loc), query)) + (data[loc] - query[0])

    plt.figure()
    plt.plot(data)
    plt.plot(query)
    plt.show()
