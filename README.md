# ucrdtw

Python extension for UCR Suite highly optimized subsequence search using Dynamic Time Warping (DTW)

Based on the paper [Searching and Mining Trillions of Time Series Subsequences under Dynamic Time Warping](http://www.cs.ucr.edu/~eamonn/SIGKDD_trillion.pdf) 

More info on the UCR Suite web page http://www.cs.ucr.edu/~eamonn/UCRsuite.html

### Requirements
Python 2.7+, numpy 1.8+

### Installation

`python setup.py build && python setup.py install`

### Usage

```
import _ucrdtw
import numpy as np
import matplotlib.pyplot as plt

data = np.cumsum(np.random.uniform(-0.5, 0.5, 1000000))
query = np.cumsum(np.random.uniform(-0.5, 0.5, 100))
loc, dist = _ucrdtw.ucrdtw(data, query, 0.05, True)
query = np.concatenate((np.linspace(0.0, 0.0, loc), query)) + (data[loc] - query[0])

plt.figure()
plt.plot(data)
plt.plot(query)
plt.show()
```
