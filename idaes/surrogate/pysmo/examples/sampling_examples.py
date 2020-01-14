##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
import pandas as pd
import os
from pyomo.common.fileutils import PYOMO_ROOT_DIR
from idaes.surrogate.pysmo import sampling as sp
from matplotlib import pyplot as plt

os.path.join(PYOMO_ROOT_DIR, 'contrib', 'surrogates', 'examples', 'data_files')

# Comparison of sampling methods with plots

data = pd.read_csv('six_hump_function_data.tab', sep='\t', header=0, index_col=0)
bounds_list = [[-1, -1], [1, 1]]

b1 = sp.LatinHypercubeSampling(data, 100, sampling_type='selection')
td11 = b1.sample_points()

b2 = sp.LatinHypercubeSampling(bounds_list, 100, sampling_type='creation')
td12 = b2.sample_points()

c1 = sp.HaltonSampling(data, 100, sampling_type='selection')
td21 = c1.sample_points()

c2 = sp.HaltonSampling(bounds_list, 100, sampling_type='creation')
td22 = c2.sample_points()

d1 = sp.HammersleySampling(data, 100, sampling_type='selection')
td31 = d1.sample_points()

d2 = sp.HammersleySampling(bounds_list, 100, sampling_type='creation')
td32 = d2.sample_points()

e1 = sp.CVTSampling(data, 100, tolerance=1e-6, sampling_type='selection')
td41 = e1.sample_points()

e2 = sp.CVTSampling(bounds_list, 100, tolerance=1e-6, sampling_type='creation')
td42 = e2.sample_points()

f1 = sp.UniformSampling(data, [10, 10], 'selection')
td51 = f1.sample_points()

f2 = sp.UniformSampling(bounds_list, [10, 10], 'creation')
td52 = f2.sample_points()

g1 = sp.UniformSampling(data, [10, 10], 'selection', False)
td61 = g1.sample_points()

g2 = sp.UniformSampling(bounds_list, [10, 10], 'creation', False)
td62 = g2.sample_points()

ax1 = plt.subplot(3, 4, 1)
ax1.plot(td11.values[:, 0], td11.values[:, 1], 'o')
ax1.grid()
ax1.set_title('LHS - selection')

ax2 = plt.subplot(3, 4, 2)
ax2.plot(td12[:, 0], td12[:, 1], 'o')
ax2.grid()
ax2.set_title('LHS - creation')

ax3 = plt.subplot(3, 4, 3)
ax3.plot(td21.values[:, 0], td21.values[:, 1], 'o')
ax3.grid()
ax3.set_title('Halton - selection')

ax4 = plt.subplot(3, 4, 4)
ax4.plot(td22[:, 0], td22[:, 1], 'o')
ax4.grid()
ax4.set_title('Halton - creation')

ax5 = plt.subplot(3, 4, 5)
ax5.plot(td31.values[:, 0], td31.values[:, 1], 'o')
ax5.grid()
ax3.set_title('Hammersley - selection')

ax6 = plt.subplot(3, 4, 6)
ax6.plot(td32[:, 0], td32[:, 1], 'o')
ax6.grid()
ax6.set_title('Hammersley - creation')

ax7 = plt.subplot(3, 4, 7)
ax7.plot(td41.values[:, 0], td41.values[:, 1], 'o')
ax7.grid()
ax7.set_title('CVT - selection')

ax8 = plt.subplot(3, 4, 8)
ax8.plot(td42[:, 0], td42[:, 1], 'o')
ax8.grid()
ax8.set_title('CVT - creation')

ax9 = plt.subplot(3, 4, 9)
ax9.plot(td51.values[:, 0], td51.values[:, 1], 'o')
ax9.grid()
ax9.set_title('Uniform (edges) - selection')

ax10 = plt.subplot(3, 4, 10)
ax10.plot(td52[:, 0], td52[:, 1], 'o')
ax10.grid()
ax10.set_title('Uniform (edges) - creation')

ax11 = plt.subplot(3, 4, 11)
ax11.plot(td61.values[:, 0], td61.values[:, 1], 'o')
ax11.grid()
ax11.set_title('Uniform (centres) - selection')

ax12 = plt.subplot(3, 4, 12)
ax12.plot(td62[:, 0], td62[:, 1], 'o')
ax12.grid()
ax12.set_title('Uniform (centres) - creation')

plt.show()
