import numpy as np
import matplotlib.pyplot as plt

from benchmarks import *

plt.ion()

# _bmSphere = BMSphere( 2, BenchmarkFunction.PLOTTING_MODE_WIREFRAME )
# _bmSphere = BMAckley( 20, 0.2, 2 * np.pi, 2, BenchmarkFunction.PLOTTING_MODE_WIREFRAME )
# _bmSphere = BMSchwefel( 2, BenchmarkFunction.PLOTTING_MODE_WIREFRAME )
_bmSphere = BMShafferFcn6( 2, BenchmarkFunction.PLOTTING_MODE_SURFACE, step = 0.1 )
_bmSphere.setRange(0.0,1.0)
_bmSphere.plotBase()


plt.show()