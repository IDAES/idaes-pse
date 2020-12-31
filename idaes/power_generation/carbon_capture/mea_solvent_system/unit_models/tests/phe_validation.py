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
"""
Comparison of  plate heat exchanger model prediction with NCCC data

Detailed model equations and data can be found in the paper :
Akula, P., Eslick, J., Bhattacharyya, D. and Miller, D.C., 2019.
Modelling and Parameter Estimation of a Plate Heat Exchanger
as Part of a Solvent-Based Post-Combustion CO2 Capture System.
In Computer Aided Chemical Engineering (Vol. 47, pp. 47-52). Elsevier.

Author: Paul Akula
"""
# Import Python libraries and third-party
import sys
import os
import matplotlib.pyplot as plt

# Import Pyomo libraries
from pyomo.environ import ConcreteModel, SolverFactory, value

# Import IDAES Libraries
from idaes.core import FlowsheetBlock

# Access the mea_solvent_system dir from the current dir (tests dir)
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from unit_models.phe import PHE
from property_package.liquid_prop import LiquidParameterBlock

solver = SolverFactory('ipopt')

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
# Set up property package
m.fs.hotside_properties = LiquidParameterBlock()
m.fs.coldside_properties = LiquidParameterBlock()

# create instance of plate heat exchanger  on flowsheet
m.fs.hx = PHE(default={
    "hot_side": {
        "property_package": m.fs.hotside_properties
    },
    "cold_side":
    {
        "property_package": m.fs.coldside_properties
    }})

for t in m.fs.time:
    # hot fluid
    m.fs.hx.hot_inlet.flow_mol[t].fix(81.6322)
    m.fs.hx.hot_inlet.temperature[t].fix(392.23)
    m.fs.hx.hot_inlet.pressure[t].fix(182840)
    m.fs.hx.hot_inlet.mole_frac_comp[t, "CO2"].fix(0.01585)
    m.fs.hx.hot_inlet.mole_frac_comp[t, "H2O"].fix(0.87457)
    m.fs.hx.hot_inlet.mole_frac_comp[t, "MEA"].fix(0.10958)

    # cold fluid
    m.fs.hx.cold_inlet.flow_mol[t].fix(84.7399)
    m.fs.hx.cold_inlet.temperature[t].fix(326.36)
    m.fs.hx.cold_inlet.pressure[t].fix(182840)
    m.fs.hx.cold_inlet.mole_frac_comp[t, "CO2"].fix(0.04142)
    m.fs.hx.cold_inlet.mole_frac_comp[t, "H2O"].fix(0.85078)
    m.fs.hx.cold_inlet.mole_frac_comp[t, "MEA"].fix(0.10780)

m.fs.hx.initialize(outlvl=0)

# model results
PHE_THOUT = []
PHE_TCOUT = []

# NCCC data
hot_molar_flowrate = [60.54879, 102.07830, 60.05750, 58.86128,
                      104.30124, 105.49475, 27.53904, 60.88287,
                      60.88904, 60.04379, 77.69829, 76.66811, 59.68288]

cold_molar_flowrate = [63.01910, 104.99350, 62.53341, 60.40000,
                       106.18207, 107.38606, 28.19399, 63.60044, 63.16421,
                       61.14541, 81.36657, 79.89472, 60.39896]

hot_Temp_IN = [392.23, 389.57, 393.78, 382.42, 376.32, 392.69, 389.69,
               392.10, 392.36, 392.30, 391.19, 390.95, 392.26]
cold_Temp_IN = [326.36, 332.26, 329.12, 318.82, 319.58, 330.54, 321.42,
                327.72, 327.47, 326.72, 328.59, 325.44, 326.11]

NCCC_hot_Temp_OUT = [330.42, 336.70, 331.44, 323.41, 324.57, 334.63,
                     324.83, 331.25, 331.08, 330.43, 332.96, 329.96, 329.75]

NCCC_cold_Temp_OUT = [384.9111111, 383.2111111, 383.6944444, 376.1722222, 370.5,
                      384.9555556, 382.0388889, 384.6388889, 384.8055556, 384.5888889,
                      384.3277778, 383.5333333, 384.5277778]

# mole fraction of components in Lean Solvent/Hot fluid
xh_H2O = [0.8747, 0.8569, 0.8739, 0.8505, 0.8573, 0.8808, 0.8582,
          0.8757, 0.8758, 0.8702, 0.8646, 0.8590, 0.8676]
xh_MEA = [0.1095, 0.1148, 0.1138, 0.1110, 0.1020, 0.1033, 0.1144,
          0.1070, 0.1071, 0.1115, 0.1106, 0.1152, 0.1134]
xh_CO2 = [0.0158, 0.0284, 0.0123, 0.0385, 0.0407, 0.0160, 0.0274,
          0.0172, 0.0171, 0.0183, 0.0248, 0.0258, 0.0190]

# mole fraction of components in Rich Solvent/Cold fluid
xc_H2O = [0.8509, 0.8426, 0.8546, 0.8378, 0.8501, 0.8676, 0.8315,
          0.8554, 0.8586, 0.8490, 0.8468, 0.8408, 0.8455]
xc_MEA = [0.1077, 0.1137, 0.1123, 0.1104, 0.1019, 0.1038, 0.1143,
          0.1050, 0.1055, 0.1110, 0.1079, 0.1127, 0.1141]
xc_CO2 = [0.0414, 0.0438, .0331, 0.0518, 0.0480, 0.0286, 0.0541,
          0.0397, 0.0360, 0.0400, 0.0453, 0.0465, 0.0404]

for t in m.fs.time:
    for i in range(len(hot_molar_flowrate)):
        # hot fluid
        m.fs.hx.hot_inlet.flow_mol[t].fix(hot_molar_flowrate[i])
        m.fs.hx.hot_inlet.temperature[t].fix(hot_Temp_IN[i])
        m.fs.hx.hot_inlet.pressure[t].fix(202650)
        m.fs.hx.hot_inlet.mole_frac_comp[t, "CO2"].fix(xh_CO2[i])
        m.fs.hx.hot_inlet.mole_frac_comp[t, "H2O"].fix(xh_H2O[i])
        m.fs.hx.hot_inlet.mole_frac_comp[t, "MEA"].fix(xh_MEA[i])

        # cold fluid
        m.fs.hx.cold_inlet.flow_mol[t].fix(cold_molar_flowrate[i])
        m.fs.hx.cold_inlet.temperature[t].fix(cold_Temp_IN[i])
        m.fs.hx.cold_inlet.pressure[t].fix(202650)
        m.fs.hx.cold_inlet.mole_frac_comp[t, "CO2"].fix(xc_CO2[i])
        m.fs.hx.cold_inlet.mole_frac_comp[t, "H2O"].fix(xc_H2O[i])
        m.fs.hx.cold_inlet.mole_frac_comp[t, "MEA"].fix(xc_MEA[i])

        solver.solve(m.fs.hx, tee=False)
        PHE_THOUT.append(value(m.fs.hx.hot_outlet.temperature[0]))
        PHE_TCOUT.append(value(m.fs.hx.cold_outlet.temperature[0]))

# plot data
fontsize = 16
labelsize = 16
markersize = 12

x = [i for i in range(1, len(NCCC_cold_Temp_OUT) + 1)]

plt.figure(figsize=(8, 6))

plt.plot(x, NCCC_cold_Temp_OUT,
         color='g',
         linestyle='',
         label='Data: Rich solvent',
         mfc="None",
         marker='s',
         markersize=markersize)
plt.plot(x, PHE_TCOUT,
         color='g',
         linestyle='',
         label='Model: Rich solvent',
         mfc="None",
         marker='x',
         markersize=markersize)
plt.plot(x, NCCC_hot_Temp_OUT,
         color='r',
         linestyle='',
         label='Data: Lean solvent',
         mfc="None",
         marker='o',
         markersize=markersize)
plt.plot(x, PHE_THOUT,
         color='r',
         linestyle='',
         label='Model: Lean solvent',
         mfc="None",
         marker='+',
         markersize=markersize)

plt.ylim(278, 400)
plt.ylabel('Temperature (K)', fontsize=fontsize, fontweight='bold')
plt.xlabel('NCCC Case No.', fontsize=fontsize, fontweight='bold')
plt.legend(loc='lower right', fontsize=fontsize)
plt.tick_params(labelsize=labelsize)
plt.tight_layout()
plt.show()
