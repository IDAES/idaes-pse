import pandas as pd
import os
from pyomo.common.fileutils import PYOMO_ROOT_DIR
from pyomo.contrib.surrogates import sampling as sp
from matplotlib import pyplot as plt

os.path.join(PYOMO_ROOT_DIR, 'contrib', 'surrogates', 'examples', 'data_files')

# Comparison of sampling methods with plots

data = pd.read_csv('six_hump_function_data.tab', sep='\t', header=0, index_col=0)


b = sp.LatinHypercubeSampling(data, 100)
td1 = b.sample_points()

c = sp.HaltonSampling(data, 100)
td2 = c.sample_points()

d = sp.HammersleySampling(data, 100)
td3 = d.sample_points()

e = sp.CVTSampling(data, 100, tolerance=1e-6)
td4 = e.sample_points()

f = sp.UniformSampling(data, [10, 10])
td5 = f.sample_points()

g = sp.UniformSampling(data, [10, 10], False)
td6 = g.sample_points()

ax1 = plt.subplot(2, 3, 1)
ax1.plot(td1.values[:, 0], td1.values[:, 1], 'o')
ax1.grid()
ax1.set_title('LHS')

ax2 = plt.subplot(2, 3, 2)
ax2.plot(td2.values[:, 0], td2.values[:, 1], 'o')
ax2.grid()
ax2.set_title('Halton')

ax3 = plt.subplot(2, 3, 3)
ax3.plot(td3.values[:, 0], td3.values[:, 1], 'o')
ax3.grid()
ax3.set_title('Hammersley')

ax4 = plt.subplot(2, 3, 4)
ax4.plot(td4.values[:, 0], td4.values[:, 1], 'o')
ax4.grid()
ax4.set_title('CVT')

ax4 = plt.subplot(2, 3, 5)
ax4.plot(td5.values[:, 0], td5.values[:, 1], 'o')
ax4.grid()
ax4.set_title('Uniform (Edges)')

ax4 = plt.subplot(2, 3, 6)
ax4.plot(td6.values[:, 0], td6.values[:, 1], 'o')
ax4.grid()
ax4.set_title('Uniform (centres)')

plt.show()
