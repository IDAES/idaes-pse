import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

from idaes.surrogate.pysmo import sampling as ns


# np_test = np.random.randn(20, 5)
# df_test = pd.DataFrame(np_test, columns=['A', 'B', 'C', 'D', 'E'])

df_test = pd.read_csv(r'rubbish_data.csv', header=0)

# Test selection
init = ns.UniformSampling(df_test, list_of_samples_per_variable=[10, 10], sampling_type='selection', xlabels=['Pressure (MPa)', 'Temperature (K)'])#, ylabels=['Enthalpy (kJ/mol)', 'Entropy (J/mol*K)'])
dt = init.sample_points()
print(dt)
plt.plot(dt['Pressure (MPa)'].values, dt['Temperature (K)'].values, 'o')
plt.show()

init2 = ns.HammersleySampling(df_test, number_of_samples=100, sampling_type='selection', xlabels=['Pressure (MPa)', 'Enthalpy (kJ/mol)'])#, ylabels=['Enthalpy (kJ/mol)', 'Entropy (J/mol*K)'])
dt2 = init2.sample_points()
print(dt2)
plt.plot(dt2['Pressure (MPa)'].values, dt2['Enthalpy (kJ/mol)'].values, 'o')
plt.show()

init3 = ns.LatinHypercubeSampling(df_test, number_of_samples=100, sampling_type='selection', xlabels=['Pressure (MPa)',  'Entropy (J/mol*K)'])#, ylabels=['Enthalpy (kJ/mol)', 'Entropy (J/mol*K)'])
dt3 = init3.sample_points()
print(dt3)
plt.plot(dt3['Pressure (MPa)'].values, dt3[ 'Entropy (J/mol*K)'].values, 'o')
plt.show()

init4 = ns.HaltonSampling(df_test, number_of_samples=100, sampling_type='selection', xlabels=['Temperature (K)',  'Entropy (J/mol*K)'])#, ylabels=['Enthalpy (kJ/mol)', 'Entropy (J/mol*K)'])
dt4 = init4.sample_points()
print(dt4)
plt.plot(dt4['Temperature (K)'].values, dt4['Entropy (J/mol*K)'].values, 'o')
plt.show()

init5 = ns.CVTSampling(df_test, number_of_samples=100, sampling_type='selection', xlabels=['Pressure (MPa)', 'Temperature (K)'], tolerance=1e-8)#, ylabels=['Enthalpy (kJ/mol)', 'Entropy (J/mol*K)'])
dt5 = init5.sample_points()
print(dt5)
plt.plot(dt5['Pressure (MPa)'].values, dt5['Temperature (K)'].values, 'o')
plt.show()




# Test selection
init = ns.UniformSampling(df_test.values, list_of_samples_per_variable=[10, 10], sampling_type='selection', xlabels=[1, 0])#, ylabels=['Enthalpy (kJ/mol)', 'Entropy (J/mol*K)'])
dt = init.sample_points()
print(dt)
plt.plot(dt[:, 0], dt[:, 1], 'o')
plt.show()


init2 = ns.HammersleySampling(df_test.values, number_of_samples=100, sampling_type='selection', xlabels=[1, 2])#, ylabels=['Enthalpy (kJ/mol)', 'Entropy (J/mol*K)'])
dt2 = init2.sample_points()
print(dt2)
plt.plot(dt2[:, 0], dt2[:, 1], 'o')
plt.show()


init3 = ns.LatinHypercubeSampling(df_test.values, number_of_samples=100, sampling_type='selection', xlabels=[1,  3])#, ylabels=['Enthalpy (kJ/mol)', 'Entropy (J/mol*K)'])
dt3 = init3.sample_points()
print(dt3)
plt.plot(dt3[:, 0], dt3[:, 1], 'o')
plt.show()


init4 = ns.HaltonSampling(df_test.values, number_of_samples=100, sampling_type='selection', xlabels=[0, 3])#, ylabels=['Enthalpy (kJ/mol)', 'Entropy (J/mol*K)'])
dt4 = init4.sample_points()
print(dt4)
plt.plot(dt4[:, 0], dt4[:, 1], 'o')
plt.show()


init5 = ns.CVTSampling(df_test.values, number_of_samples=100, sampling_type='selection', xlabels=[1, 0], tolerance=1e-8)#, ylabels=['Enthalpy (kJ/mol)', 'Entropy (J/mol*K)'])
dt5 = init5.sample_points()
print(dt5)
plt.plot(dt5[:, 0], dt5[:, 1], 'o')
plt.show()