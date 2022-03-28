#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
data_dic = {
    1: "(74904.4 * b.wall_temperature_waterwall[t, 1] \
            +11301.8 * b.wall_temperature_waterwall[t, 2] \
            +2427.54 * b.wall_temperature_waterwall[t, 3] \
            +2891.35 * b.wall_temperature_waterwall[t, 4] \
            -28320.8 * b.wall_temperature_waterwall[t, 5] \
            +171.944 * b.wall_temperature_waterwall[t, 6] \
            +14462.9 * b.wall_temperature_waterwall[t, 7] \
            +677.973 * b.wall_temperature_waterwall[t, 8] \
            -122.598 * b.wall_temperature_waterwall[t, 9] \
            -300.609 * b.wall_temperature_waterwall[t, 12] \
            +1770.05 * b.wall_temperature_platen[t] \
            -169.562 * b.wall_temperature_roof[t] \
            +1.92447e+06 * b.flowrate_coal_raw[t] \
            +8.65298e+07 * b.mf_H2O_coal_raw[t] \
            +1.46803e+09 * b.SR[t] \
            -1.68637e+08 * b.SR_lf[t] \
            -6385.5 * b.secondary_air_inlet.temperature[t] \
            +952262 * b.ratio_PA2coal[t] \
            -2.10696e+07 * log(b.wall_temperature_waterwall[t, 1]) \
            -5.56721e+06 * log(b.wall_temperature_waterwall[t, 2]) \
            -1.21448e+06 * log(b.wall_temperature_waterwall[t, 3]) \
            -3.48806e+06 * log(b.wall_temperature_waterwall[t, 4]) \
            +1.02069e+07 * log(b.wall_temperature_waterwall[t, 5]) \
            -5.47775e+06 * log(b.wall_temperature_waterwall[t, 7]) \
            -279177 * log(b.wall_temperature_waterwall[t, 8]) \
            +111061 * log(b.wall_temperature_waterwall[t, 9]) \
            -646506 * log(b.wall_temperature_platen[t]) \
            +3.96004e+06 * log(b.flowrate_coal_raw[t]) \
            +45092.6 * log(b.mf_H2O_coal_raw[t]) \
            -3.9193e+08 * log(b.SR[t]) \
            +1.61406e+08 * log(b.SR_lf[t]) \
            -3.84139e+06 * log(b.secondary_air_inlet.temperature[t]) \
            +556876 * log(b.ratio_PA2coal[t]) \
            -9.88017e+07 * exp(b.mf_H2O_coal_raw[t]) \
            -2.86022e+08 * exp(b.SR[t]) \
            -35.7424 * b.wall_temperature_waterwall[t, 1]**2 \
            +10.3165 * b.wall_temperature_waterwall[t, 5]**2 \
            -4.85009 * b.wall_temperature_waterwall[t, 7]**2 \
            -21190.8 * b.flowrate_coal_raw[t]**2 \
            -5.00302e+08 * b.SR[t]**2 \
            +179.062 * b.flowrate_coal_raw[t]**3 \
            +2.33815e+08 * b.SR[t]**3 \
            -74.5573 * b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t] \
            -0.202341 * b.wall_temperature_waterwall[t, 2]*b.wall_temperature_waterwall[t, 3] \
            -52.6317 * b.wall_temperature_waterwall[t, 2]*b.flowrate_coal_raw[t] \
            +0.243508 * b.wall_temperature_waterwall[t, 2]*b.secondary_air_inlet.temperature[t] \
            -39.6689 * b.wall_temperature_waterwall[t, 2]*b.ratio_PA2coal[t] \
            -0.320972 * b.wall_temperature_waterwall[t, 3]*b.wall_temperature_waterwall[t, 4] \
            -9.03032 * b.wall_temperature_waterwall[t, 3]*b.flowrate_coal_raw[t] \
            +1.08186 * b.wall_temperature_waterwall[t, 3]*b.secondary_air_inlet.temperature[t] \
            +0.525868 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_roof[t] \
            -39.011 * b.wall_temperature_waterwall[t, 4]*b.flowrate_coal_raw[t] \
            +4489.28 * b.wall_temperature_waterwall[t, 4]*b.SR[t] \
            -629.837 * b.wall_temperature_waterwall[t, 5]*b.SR_lf[t] \
            +44.344 * b.wall_temperature_waterwall[t, 5]*b.ratio_PA2coal[t] \
            +1.72535 * b.wall_temperature_waterwall[t, 7]*b.flowrate_coal_raw[t] \
            +796.63 * b.wall_temperature_waterwall[t, 7]*b.mf_H2O_coal_raw[t] \
            +0.282511 * b.wall_temperature_waterwall[t, 7]*b.secondary_air_inlet.temperature[t] \
            +0.447853 * b.wall_temperature_waterwall[t, 8]*b.wall_temperature_waterwall[t, 12] \
            -0.343465 * b.wall_temperature_waterwall[t, 8]*b.wall_temperature_platen[t] \
            -32.6091 * b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t] \
            +1107.94 * b.wall_temperature_waterwall[t, 8]*b.mf_H2O_coal_raw[t] \
            +0.0605128 * b.wall_temperature_waterwall[t, 10]*b.wall_temperature_waterwall[t, 11] \
            -0.463551 * b.wall_temperature_platen[t]*b.wall_temperature_roof[t] \
            +6.50059 * b.wall_temperature_platen[t]*b.flowrate_coal_raw[t] \
            -1051.37 * b.wall_temperature_platen[t]*b.mf_H2O_coal_raw[t] \
            -6.1049 * b.wall_temperature_platen[t]*b.ratio_PA2coal[t] \
            -5465.56 * b.wall_temperature_roof[t]*b.mf_H2O_coal_raw[t] \
            +0.833308 * b.wall_temperature_roof[t]*b.secondary_air_inlet.temperature[t] \
            -1.75056e+06 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            -183240 * b.flowrate_coal_raw[t]*b.SR[t] \
            -310130 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            +957.014 * b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -39804.2 * b.flowrate_coal_raw[t]*b.ratio_PA2coal[t] \
            +4.49631e+06 * b.mf_H2O_coal_raw[t]*b.SR[t] \
            +2.16414e+07 * b.mf_H2O_coal_raw[t]*b.SR_lf[t] \
            -19740 * b.mf_H2O_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -4.21302e+06 * b.SR[t]*b.SR_lf[t] \
            +132221 * b.SR[t]*b.ratio_PA2coal[t] \
            +20016.4 * b.SR_lf[t]*b.secondary_air_inlet.temperature[t] \
            -2318.51 * b.secondary_air_inlet.temperature[t]*b.ratio_PA2coal[t] \
            +0.000979408 * (b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t])**2 \
            +0.000719763 * (b.wall_temperature_waterwall[t, 2]*b.flowrate_coal_raw[t])**2 \
            +0.000572793 * (b.wall_temperature_waterwall[t, 4]*b.flowrate_coal_raw[t])**2 \
            -1.27643 * (b.wall_temperature_waterwall[t, 4]*b.SR[t])**2 \
            +0.000450081 * (b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t])**2 \
            +140272 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**2 \
            +666.822 * (b.flowrate_coal_raw[t]*b.SR[t])**2 \
            -3994.96 * (b.flowrate_coal_raw[t]*b.SR_lf[t])**2 \
            -0.00660786 * (b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t])**2 \
            +85.0348 * (b.flowrate_coal_raw[t]*b.ratio_PA2coal[t])**2 \
            +0.102891 * (b.wall_temperature_roof[t]*b.mf_H2O_coal_raw[t])**3 \
            -7838.06 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**3)",
    2: "(9402.37 * b.wall_temperature_waterwall[t, 1] \
            +63618 * b.wall_temperature_waterwall[t, 2] \
            +2118.88 * b.wall_temperature_waterwall[t, 3] \
            +1821.84 * b.wall_temperature_waterwall[t, 4] \
            -24586.5 * b.wall_temperature_waterwall[t, 5] \
            +166.046 * b.wall_temperature_waterwall[t, 6] \
            +695.921 * b.wall_temperature_waterwall[t, 7] \
            +1070.12 * b.wall_temperature_waterwall[t, 8] \
            -231.135 * b.wall_temperature_waterwall[t, 9] \
            +18.6849 * b.wall_temperature_waterwall[t, 12] \
            +1794.21 * b.wall_temperature_platen[t] \
            +2842.9 * b.wall_temperature_roof[t] \
            +2.04728e+06 * b.flowrate_coal_raw[t] \
            +9.13283e+07 * b.mf_H2O_coal_raw[t] \
            +6.91793e+06 * b.SR[t] \
            +4.81359e+08 * b.SR_lf[t] \
            -40341.5 * b.secondary_air_inlet.temperature[t] \
            +1.1322e+06 * b.ratio_PA2coal[t] \
            -4.24844e+06 * log(b.wall_temperature_waterwall[t, 1]) \
            -1.67693e+07 * log(b.wall_temperature_waterwall[t, 2]) \
            -1.2861e+06 * log(b.wall_temperature_waterwall[t, 3]) \
            -4.45769e+06 * log(b.wall_temperature_waterwall[t, 4]) \
            +9.0161e+06 * log(b.wall_temperature_waterwall[t, 5]) \
            -537334 * log(b.wall_temperature_waterwall[t, 8]) \
            +181486 * log(b.wall_temperature_waterwall[t, 9]) \
            -715084 * log(b.wall_temperature_platen[t]) \
            +4.61896e+06 * log(b.flowrate_coal_raw[t]) \
            +30243.2 * log(b.mf_H2O_coal_raw[t]) \
            -8.42989e+06 * log(b.SR[t]) \
            -1.45952e+08 * log(b.SR_lf[t]) \
            +9.65255e+06 * log(b.secondary_air_inlet.temperature[t]) \
            +490346 * log(b.ratio_PA2coal[t]) \
            -1.04352e+08 * exp(b.mf_H2O_coal_raw[t]) \
            -1.24828e+08 * exp(b.SR_lf[t]) \
            -30.804 * b.wall_temperature_waterwall[t, 2]**2 \
            +9.2668 * b.wall_temperature_waterwall[t, 5]**2 \
            -21962.3 * b.flowrate_coal_raw[t]**2 \
            +185.889 * b.flowrate_coal_raw[t]**3 \
            +0.0096263 * b.secondary_air_inlet.temperature[t]**3 \
            -76.1054 * b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t] \
            -81.3516 * b.wall_temperature_waterwall[t, 2]*b.flowrate_coal_raw[t] \
            +1570.34 * b.wall_temperature_waterwall[t, 2]*b.mf_H2O_coal_raw[t] \
            -2440.82 * b.wall_temperature_waterwall[t, 2]*b.SR_lf[t] \
            +0.527827 * b.wall_temperature_waterwall[t, 2]*b.secondary_air_inlet.temperature[t] \
            -31.6806 * b.wall_temperature_waterwall[t, 2]*b.ratio_PA2coal[t] \
            -11.9902 * b.wall_temperature_waterwall[t, 3]*b.flowrate_coal_raw[t] \
            +1.17551 * b.wall_temperature_waterwall[t, 3]*b.secondary_air_inlet.temperature[t] \
            +0.712584 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_roof[t] \
            -8.98673 * b.wall_temperature_waterwall[t, 4]*b.flowrate_coal_raw[t] \
            +6307.56 * b.wall_temperature_waterwall[t, 4]*b.SR[t] \
            +1.03784 * b.wall_temperature_waterwall[t, 4]*b.secondary_air_inlet.temperature[t] \
            -3.53159 * b.wall_temperature_waterwall[t, 5]*b.flowrate_coal_raw[t] \
            -1225.87 * b.wall_temperature_waterwall[t, 5]*b.SR_lf[t] \
            +85.2083 * b.wall_temperature_waterwall[t, 5]*b.ratio_PA2coal[t] \
            +656.338 * b.wall_temperature_waterwall[t, 7]*b.mf_H2O_coal_raw[t] \
            -674.09 * b.wall_temperature_waterwall[t, 7]*b.SR_lf[t] \
            -0.265762 * b.wall_temperature_waterwall[t, 8]*b.wall_temperature_platen[t] \
            -8.30267 * b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t] \
            +1227.37 * b.wall_temperature_waterwall[t, 8]*b.mf_H2O_coal_raw[t] \
            +0.0526552 * b.wall_temperature_waterwall[t, 10]*b.wall_temperature_waterwall[t, 11] \
            -0.407706 * b.wall_temperature_platen[t]*b.wall_temperature_roof[t] \
            +7.81514 * b.wall_temperature_platen[t]*b.flowrate_coal_raw[t] \
            -1118.85 * b.wall_temperature_platen[t]*b.mf_H2O_coal_raw[t] \
            -15.2498 * b.wall_temperature_platen[t]*b.ratio_PA2coal[t] \
            -4845.51 * b.wall_temperature_roof[t]*b.mf_H2O_coal_raw[t] \
            -3234.16 * b.wall_temperature_roof[t]*b.SR_lf[t] \
            +0.69596 * b.wall_temperature_roof[t]*b.secondary_air_inlet.temperature[t] \
            -1.8204e+06 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            -159924 * b.flowrate_coal_raw[t]*b.SR[t] \
            -398436 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            +1029.33 * b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -42896.7 * b.flowrate_coal_raw[t]*b.ratio_PA2coal[t] \
            +3.55413e+06 * b.mf_H2O_coal_raw[t]*b.SR[t] \
            +2.26223e+07 * b.mf_H2O_coal_raw[t]*b.SR_lf[t] \
            -20279.8 * b.mf_H2O_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -3.77353e+06 * b.SR[t]*b.SR_lf[t] \
            +308.059 * b.SR[t]*b.secondary_air_inlet.temperature[t] \
            +107358 * b.SR[t]*b.ratio_PA2coal[t] \
            +20666.6 * b.SR_lf[t]*b.secondary_air_inlet.temperature[t] \
            -2564.62 * b.secondary_air_inlet.temperature[t]*b.ratio_PA2coal[t] \
            +0.00107543 * (b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t])**2 \
            +0.000875275 * (b.wall_temperature_waterwall[t, 2]*b.flowrate_coal_raw[t])**2 \
            -1.75451 * (b.wall_temperature_waterwall[t, 4]*b.SR[t])**2 \
            +149896 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**2 \
            +654.854 * (b.flowrate_coal_raw[t]*b.SR[t])**2 \
            -3875.11 * (b.flowrate_coal_raw[t]*b.SR_lf[t])**2 \
            -0.00723196 * (b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t])**2 \
            +91.9401 * (b.flowrate_coal_raw[t]*b.ratio_PA2coal[t])**2 \
            +0.0961783 * (b.wall_temperature_roof[t]*b.mf_H2O_coal_raw[t])**3 \
            -8608.57 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**3)",
    3: "(2825.25 * b.wall_temperature_waterwall[t, 1] \
            +3496.16 * b.wall_temperature_waterwall[t, 2] \
            +39611.4 * b.wall_temperature_waterwall[t, 3] \
            +2064.54 * b.wall_temperature_waterwall[t, 4] \
            -19353.7 * b.wall_temperature_waterwall[t, 5] \
            +140.98 * b.wall_temperature_waterwall[t, 6] \
            -95.4723 * b.wall_temperature_waterwall[t, 7] \
            +180.791 * b.wall_temperature_waterwall[t, 8] \
            +20.6861 * b.wall_temperature_waterwall[t, 12] \
            -77.5609 * b.wall_temperature_platen[t] \
            +1.76659e+06 * b.flowrate_coal_raw[t] \
            +2.73809e+07 * b.mf_H2O_coal_raw[t] \
            +1.30314e+07 * b.SR[t] \
            +2.99323e+08 * b.SR_lf[t] \
            -8992.63 * b.secondary_air_inlet.temperature[t] \
            +743538 * b.ratio_PA2coal[t] \
            -1.21925e+06 * log(b.wall_temperature_waterwall[t, 1]) \
            -1.70788e+06 * log(b.wall_temperature_waterwall[t, 2]) \
            -1.11283e+07 * log(b.wall_temperature_waterwall[t, 3]) \
            -1.15967e+06 * log(b.wall_temperature_waterwall[t, 4]) \
            +6.95686e+06 * log(b.wall_temperature_waterwall[t, 5]) \
            +3.12918e+06 * log(b.flowrate_coal_raw[t]) \
            +86928.8 * log(b.mf_H2O_coal_raw[t]) \
            -1.11133e+07 * log(b.SR[t]) \
            -8.27815e+07 * log(b.SR_lf[t]) \
            +317858 * log(b.ratio_PA2coal[t]) \
            -4.84182e+07 * exp(b.mf_H2O_coal_raw[t]) \
            -969324 * exp(b.SR[t]) \
            -8.1132e+07 * exp(b.SR_lf[t]) \
            -20.0424 * b.wall_temperature_waterwall[t, 3]**2 \
            +7.21081 * b.wall_temperature_waterwall[t, 5]**2 \
            +0.177218 * b.wall_temperature_waterwall[t, 8]**2 \
            +0.469474 * b.wall_temperature_platen[t]**2 \
            -19318.5 * b.flowrate_coal_raw[t]**2 \
            +3.33513 * b.secondary_air_inlet.temperature[t]**2 \
            +149.518 * b.flowrate_coal_raw[t]**3 \
            -29.4453 * b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t] \
            -32.348 * b.wall_temperature_waterwall[t, 2]*b.flowrate_coal_raw[t] \
            +0.37361 * b.wall_temperature_waterwall[t, 2]*b.secondary_air_inlet.temperature[t] \
            -49.1179 * b.wall_temperature_waterwall[t, 2]*b.ratio_PA2coal[t] \
            -0.0559492 * b.wall_temperature_waterwall[t, 3]*b.wall_temperature_waterwall[t, 8] \
            -24.6599 * b.wall_temperature_waterwall[t, 3]*b.flowrate_coal_raw[t] \
            +1.0472 * b.wall_temperature_waterwall[t, 3]*b.secondary_air_inlet.temperature[t] \
            +0.051371 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_waterwall[t, 8] \
            +0.291925 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_roof[t] \
            -7.09556 * b.wall_temperature_waterwall[t, 4]*b.flowrate_coal_raw[t] \
            -579.937 * b.wall_temperature_waterwall[t, 5]*b.SR_lf[t] \
            +30.4256 * b.wall_temperature_waterwall[t, 5]*b.ratio_PA2coal[t] \
            +0.187627 * b.wall_temperature_waterwall[t, 7]*b.wall_temperature_roof[t] \
            +506.243 * b.wall_temperature_waterwall[t, 7]*b.mf_H2O_coal_raw[t] \
            -0.209174 * b.wall_temperature_waterwall[t, 8]*b.wall_temperature_platen[t] \
            -27.8211 * b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t] \
            +808.58 * b.wall_temperature_waterwall[t, 8]*b.mf_H2O_coal_raw[t] \
            +0.0378872 * b.wall_temperature_waterwall[t, 10]*b.wall_temperature_waterwall[t, 11] \
            -0.427984 * b.wall_temperature_platen[t]*b.wall_temperature_roof[t] \
            +3.97146 * b.wall_temperature_platen[t]*b.flowrate_coal_raw[t] \
            -2.56524 * b.wall_temperature_platen[t]*b.ratio_PA2coal[t] \
            -1.1316e+06 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            -149208 * b.flowrate_coal_raw[t]*b.SR[t] \
            -620763 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            +692.429 * b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -27383.3 * b.flowrate_coal_raw[t]*b.ratio_PA2coal[t] \
            +2.74367e+06 * b.mf_H2O_coal_raw[t]*b.SR[t] \
            +2.19307e+07 * b.mf_H2O_coal_raw[t]*b.SR_lf[t] \
            -12923.6 * b.mf_H2O_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -1.89277e+06 * b.SR[t]*b.SR_lf[t] \
            +538.573 * b.SR[t]*b.secondary_air_inlet.temperature[t] \
            +91705.4 * b.SR[t]*b.ratio_PA2coal[t] \
            +9695.45 * b.SR_lf[t]*b.secondary_air_inlet.temperature[t] \
            -1698.89 * b.secondary_air_inlet.temperature[t]*b.ratio_PA2coal[t] \
            +0.000466813 * (b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t])**2 \
            +0.000489048 * (b.wall_temperature_waterwall[t, 2]*b.flowrate_coal_raw[t])**2 \
            -0.0629385 * (b.wall_temperature_waterwall[t, 4]*b.SR[t])**2 \
            +0.000414 * (b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t])**2 \
            +94334.8 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**2 \
            +647.348 * (b.flowrate_coal_raw[t]*b.SR[t])**2 \
            -0.00513905 * (b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t])**2 \
            +57.4901 * (b.flowrate_coal_raw[t]*b.ratio_PA2coal[t])**2 \
            -5271.02 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**3)",
    4: "(30998 * b.wall_temperature_waterwall[t, 1] \
            +2457.99 * b.wall_temperature_waterwall[t, 2] \
            +1137.13 * b.wall_temperature_waterwall[t, 3] \
            +19608.1 * b.wall_temperature_waterwall[t, 4] \
            -19918.2 * b.wall_temperature_waterwall[t, 5] \
            +919.01 * b.wall_temperature_waterwall[t, 6] \
            +369.947 * b.wall_temperature_waterwall[t, 7] \
            +646.454 * b.wall_temperature_waterwall[t, 8] \
            +102.897 * b.wall_temperature_waterwall[t, 10] \
            +24.8228 * b.wall_temperature_waterwall[t, 12] \
            +1739.59 * b.wall_temperature_platen[t] \
            +1.89204e+06 * b.flowrate_coal_raw[t] \
            +6.68693e+07 * b.mf_H2O_coal_raw[t] \
            +4.63505e+06 * b.SR[t] \
            +4.53751e+08 * b.SR_lf[t] \
            +1612.6 * b.secondary_air_inlet.temperature[t] \
            +924842 * b.ratio_PA2coal[t] \
            -1.07435e+07 * log(b.wall_temperature_waterwall[t, 1]) \
            -1.48504e+06 * log(b.wall_temperature_waterwall[t, 2]) \
            -750102 * log(b.wall_temperature_waterwall[t, 3]) \
            -7.28634e+06 * log(b.wall_temperature_waterwall[t, 4]) \
            +7.19016e+06 * log(b.wall_temperature_waterwall[t, 5]) \
            -342011 * log(b.wall_temperature_waterwall[t, 6]) \
            -465932 * log(b.wall_temperature_waterwall[t, 7]) \
            -275633 * log(b.wall_temperature_waterwall[t, 8]) \
            -500395 * log(b.wall_temperature_platen[t]) \
            +3.33092e+06 * log(b.flowrate_coal_raw[t]) \
            -1.06661e+07 * log(b.SR[t]) \
            -1.71633e+08 * log(b.SR_lf[t]) \
            -3.18827e+06 * log(b.secondary_air_inlet.temperature[t]) \
            -8.32439e+07 * exp(b.mf_H2O_coal_raw[t]) \
            -1.05477e+08 * exp(b.SR_lf[t]) \
            -10.3499 * b.wall_temperature_waterwall[t, 1]**2 \
            -12.7144 * b.wall_temperature_waterwall[t, 4]**2 \
            +7.8685 * b.wall_temperature_waterwall[t, 5]**2 \
            -19144.9 * b.flowrate_coal_raw[t]**2 \
            +143.33 * b.flowrate_coal_raw[t]**3 \
            -24.1862 * b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t] \
            -28.4264 * b.wall_temperature_waterwall[t, 1]*b.mf_H2O_coal_raw[t] \
            -312.37 * b.wall_temperature_waterwall[t, 1]*b.SR_lf[t] \
            -3.62923 * b.wall_temperature_waterwall[t, 2]*b.flowrate_coal_raw[t] \
            +0.58587 * b.wall_temperature_waterwall[t, 2]*b.secondary_air_inlet.temperature[t] \
            -20.9999 * b.wall_temperature_waterwall[t, 2]*b.ratio_PA2coal[t] \
            -4.26059 * b.wall_temperature_waterwall[t, 3]*b.flowrate_coal_raw[t] \
            +0.849726 * b.wall_temperature_waterwall[t, 3]*b.secondary_air_inlet.temperature[t] \
            +0.383348 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_roof[t] \
            -31.1716 * b.wall_temperature_waterwall[t, 4]*b.flowrate_coal_raw[t] \
            +6795.15 * b.wall_temperature_waterwall[t, 4]*b.SR[t] \
            -0.334707 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 6] \
            +0.0296533 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_roof[t] \
            -6.75789 * b.wall_temperature_waterwall[t, 5]*b.flowrate_coal_raw[t] \
            -887.887 * b.wall_temperature_waterwall[t, 5]*b.SR_lf[t] \
            +54.121 * b.wall_temperature_waterwall[t, 5]*b.ratio_PA2coal[t] \
            +0.390219 * b.wall_temperature_waterwall[t, 7]*b.wall_temperature_roof[t] \
            +617.963 * b.wall_temperature_waterwall[t, 7]*b.mf_H2O_coal_raw[t] \
            -354.236 * b.wall_temperature_waterwall[t, 7]*b.SR_lf[t] \
            +0.668885 * b.wall_temperature_waterwall[t, 7]*b.secondary_air_inlet.temperature[t] \
            +0.368313 * b.wall_temperature_waterwall[t, 8]*b.wall_temperature_waterwall[t, 10] \
            -0.224167 * b.wall_temperature_waterwall[t, 8]*b.wall_temperature_platen[t] \
            -32.5371 * b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t] \
            +944.267 * b.wall_temperature_waterwall[t, 8]*b.mf_H2O_coal_raw[t] \
            -0.444565 * b.wall_temperature_waterwall[t, 10]*b.wall_temperature_platen[t] \
            +482.957 * b.wall_temperature_waterwall[t, 11]*b.mf_H2O_coal_raw[t] \
            -0.448545 * b.wall_temperature_platen[t]*b.wall_temperature_roof[t] \
            -92.9318 * b.wall_temperature_roof[t]*b.ratio_PA2coal[t] \
            -965339 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            -238814 * b.flowrate_coal_raw[t]*b.SR[t] \
            -656822 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            +742.242 * b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -26766.3 * b.flowrate_coal_raw[t]*b.ratio_PA2coal[t] \
            -1.27395e+06 * b.mf_H2O_coal_raw[t]*b.SR[t] \
            +2.30581e+07 * b.mf_H2O_coal_raw[t]*b.SR_lf[t] \
            -13952.1 * b.mf_H2O_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -269639 * b.mf_H2O_coal_raw[t]*b.ratio_PA2coal[t] \
            +841.556 * b.SR[t]*b.secondary_air_inlet.temperature[t] \
            +84459.4 * b.SR[t]*b.ratio_PA2coal[t] \
            +7494.08 * b.SR_lf[t]*b.secondary_air_inlet.temperature[t] \
            -1672.75 * b.secondary_air_inlet.temperature[t]*b.ratio_PA2coal[t] \
            +0.000370175 * (b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t])**2 \
            -1.84302 * (b.wall_temperature_waterwall[t, 4]*b.SR[t])**2 \
            +0.000463608 * (b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t])**2 \
            +45660.1 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**2 \
            +978.377 * (b.flowrate_coal_raw[t]*b.SR[t])**2 \
            -0.00577998 * (b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t])**2 \
            +51.3598 * (b.flowrate_coal_raw[t]*b.ratio_PA2coal[t])**2 \
            +1.03552e+07 * (b.mf_H2O_coal_raw[t]*b.SR[t])**2)",
    5: "(1381.36 * b.wall_temperature_waterwall[t, 1] \
            +2323.67 * b.wall_temperature_waterwall[t, 2] \
            +16.5773 * b.wall_temperature_waterwall[t, 3] \
            +1352.22 * b.wall_temperature_waterwall[t, 4] \
            -16466.5 * b.wall_temperature_waterwall[t, 5] \
            +899.714 * b.wall_temperature_waterwall[t, 6] \
            +3.93998 * b.wall_temperature_waterwall[t, 7] \
            +1266.5 * b.wall_temperature_waterwall[t, 8] \
            +27.0207 * b.wall_temperature_waterwall[t, 12] \
            +1538.23 * b.wall_temperature_platen[t] \
            +2873.09 * b.wall_temperature_roof[t] \
            +1.97605e+06 * b.flowrate_coal_raw[t] \
            -1.23146e+08 * b.mf_H2O_coal_raw[t] \
            +7.90727e+06 * b.SR[t] \
            -8.28156e+07 * b.SR_lf[t] \
            +8409.65 * b.secondary_air_inlet.temperature[t] \
            +771488 * b.ratio_PA2coal[t] \
            -512140 * log(b.wall_temperature_waterwall[t, 1]) \
            -960400 * log(b.wall_temperature_waterwall[t, 2]) \
            -1.69396e+06 * log(b.wall_temperature_waterwall[t, 5]) \
            -578390 * log(b.wall_temperature_waterwall[t, 8]) \
            -500164 * log(b.wall_temperature_platen[t]) \
            +3.27151e+06 * log(b.flowrate_coal_raw[t]) \
            -207751 * log(b.mf_H2O_coal_raw[t]) \
            -1.3498e+07 * log(b.SR[t]) \
            +7.59142e+07 * log(b.SR_lf[t]) \
            -3.4224e+06 * log(b.secondary_air_inlet.temperature[t]) \
            +1.08294e+08 * exp(b.mf_H2O_coal_raw[t]) \
            +0.728811 * b.wall_temperature_waterwall[t, 4]**2 \
            -18833.1 * b.flowrate_coal_raw[t]**2 \
            -9.85614e+07 * b.mf_H2O_coal_raw[t]**2 \
            +120.533 * b.flowrate_coal_raw[t]**3 \
            -0.0492326 * b.wall_temperature_waterwall[t, 1]*b.wall_temperature_waterwall[t, 11] \
            -22.6501 * b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t] \
            -0.453682 * b.wall_temperature_waterwall[t, 2]*b.wall_temperature_waterwall[t, 5] \
            -0.372273 * b.wall_temperature_waterwall[t, 2]*b.wall_temperature_platen[t] \
            +0.420252 * b.wall_temperature_waterwall[t, 2]*b.secondary_air_inlet.temperature[t] \
            -78.644 * b.wall_temperature_waterwall[t, 2]*b.ratio_PA2coal[t] \
            -0.548922 * b.wall_temperature_waterwall[t, 3]*b.wall_temperature_waterwall[t, 4] \
            +0.782376 * b.wall_temperature_waterwall[t, 3]*b.wall_temperature_waterwall[t, 11] \
            +176.927 * b.wall_temperature_waterwall[t, 3]*b.SR[t] \
            -0.0216351 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_waterwall[t, 8] \
            +0.521495 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_roof[t] \
            -23.7278 * b.wall_temperature_waterwall[t, 4]*b.flowrate_coal_raw[t] \
            -354.669 * b.wall_temperature_waterwall[t, 4]*b.SR[t] \
            -1800.24 * b.wall_temperature_waterwall[t, 4]*b.SR_lf[t] \
            +1.00785 * b.wall_temperature_waterwall[t, 4]*b.secondary_air_inlet.temperature[t] \
            -0.320146 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 6] \
            -36.3443 * b.wall_temperature_waterwall[t, 5]*b.flowrate_coal_raw[t] \
            -276.22 * b.wall_temperature_waterwall[t, 5]*b.SR[t] \
            +31986.7 * b.wall_temperature_waterwall[t, 5]*b.SR_lf[t] \
            -33.1281 * b.wall_temperature_waterwall[t, 6]*b.flowrate_coal_raw[t] \
            +0.337834 * b.wall_temperature_waterwall[t, 7]*b.secondary_air_inlet.temperature[t] \
            -22.6854 * b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t] \
            +0.0949332 * b.wall_temperature_waterwall[t, 9]*b.wall_temperature_waterwall[t, 10] \
            -0.836476 * b.wall_temperature_waterwall[t, 11]*b.secondary_air_inlet.temperature[t] \
            -0.524954 * b.wall_temperature_platen[t]*b.wall_temperature_roof[t] \
            +3.06373 * b.wall_temperature_platen[t]*b.flowrate_coal_raw[t] \
            +424.529 * b.wall_temperature_roof[t]*b.mf_H2O_coal_raw[t] \
            -2856.66 * b.wall_temperature_roof[t]*b.SR_lf[t] \
            -983765 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            -407939 * b.flowrate_coal_raw[t]*b.SR[t] \
            -532340 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            +744.637 * b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -24675.1 * b.flowrate_coal_raw[t]*b.ratio_PA2coal[t] \
            +6.59033e+06 * b.mf_H2O_coal_raw[t]*b.SR[t] \
            +1.76578e+07 * b.mf_H2O_coal_raw[t]*b.SR_lf[t] \
            -14012.2 * b.mf_H2O_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            +2047.81 * b.SR[t]*b.secondary_air_inlet.temperature[t] \
            +136654 * b.SR[t]*b.ratio_PA2coal[t] \
            -1520.71 * b.secondary_air_inlet.temperature[t]*b.ratio_PA2coal[t] \
            +0.000356587 * (b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t])**2 \
            +0.000348091 * (b.wall_temperature_waterwall[t, 4]*b.flowrate_coal_raw[t])**2 \
            -11.5928 * (b.wall_temperature_waterwall[t, 5]*b.SR_lf[t])**2 \
            +0.000500883 * (b.wall_temperature_waterwall[t, 6]*b.flowrate_coal_raw[t])**2 \
            +0.000224347 * (b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t])**2 \
            -0.0108347 * (b.wall_temperature_roof[t]*b.ratio_PA2coal[t])**2 \
            +46556.5 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**2 \
            +1577.66 * (b.flowrate_coal_raw[t]*b.SR[t])**2 \
            -0.00554362 * (b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t])**2 \
            +43.5639 * (b.flowrate_coal_raw[t]*b.ratio_PA2coal[t])**2)",
    6: "(23766.1 * b.wall_temperature_waterwall[t, 1] \
            +1109.67 * b.wall_temperature_waterwall[t, 2] \
            -411.801 * b.wall_temperature_waterwall[t, 3] \
            +1277.07 * b.wall_temperature_waterwall[t, 4] \
            +1664.61 * b.wall_temperature_waterwall[t, 5] \
            -144143 * b.wall_temperature_waterwall[t, 6] \
            +321.352 * b.wall_temperature_waterwall[t, 7] \
            +2667.15 * b.wall_temperature_waterwall[t, 8] \
            +812.223 * b.wall_temperature_waterwall[t, 9] \
            +278.07 * b.wall_temperature_waterwall[t, 10] \
            -357.846 * b.wall_temperature_waterwall[t, 11] \
            +955.339 * b.wall_temperature_waterwall[t, 12] \
            +1101.44 * b.wall_temperature_platen[t] \
            -1083.91 * b.wall_temperature_roof[t] \
            +1.87747e+06 * b.flowrate_coal_raw[t] \
            +6.52201e+07 * b.mf_H2O_coal_raw[t] \
            +1.15786e+07 * b.SR[t] \
            +78146.6 * b.SR_lf[t] \
            +8886.98 * b.secondary_air_inlet.temperature[t] \
            +349860 * b.ratio_PA2coal[t] \
            -1.03956e+07 * log(b.wall_temperature_waterwall[t, 1]) \
            +410736 * log(b.wall_temperature_waterwall[t, 3]) \
            +2.84904e+07 * log(b.wall_temperature_waterwall[t, 6]) \
            -982762 * log(b.wall_temperature_waterwall[t, 8]) \
            -638942 * log(b.wall_temperature_platen[t]) \
            +2.51419e+06 * log(b.flowrate_coal_raw[t]) \
            -355831 * log(b.mf_H2O_coal_raw[t]) \
            -3.22971e+07 * log(b.SR[t]) \
            -3.53833e+06 * log(b.secondary_air_inlet.temperature[t]) \
            -7.27122e+07 * exp(b.mf_H2O_coal_raw[t]) \
            +121.792 * b.wall_temperature_waterwall[t, 6]**2 \
            -17167.1 * b.flowrate_coal_raw[t]**2 \
            -0.00511324 * b.wall_temperature_waterwall[t, 1]**3 \
            -0.0481974 * b.wall_temperature_waterwall[t, 6]**3 \
            +74.4938 * b.flowrate_coal_raw[t]**3 \
            -0.919796 * b.wall_temperature_waterwall[t, 1]*b.wall_temperature_waterwall[t, 8] \
            +0.0606529 * b.wall_temperature_waterwall[t, 1]*b.wall_temperature_waterwall[t, 10] \
            -46.0797 * b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t] \
            -0.499405 * b.wall_temperature_waterwall[t, 2]*b.wall_temperature_waterwall[t, 5] \
            -181.556 * b.wall_temperature_waterwall[t, 2]*b.ratio_PA2coal[t] \
            +116.97 * b.wall_temperature_waterwall[t, 3]*b.SR[t] \
            +0.393401 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_waterwall[t, 6] \
            +0.675614 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_waterwall[t, 10] \
            -753.901 * b.wall_temperature_waterwall[t, 4]*b.SR[t] \
            -2042.45 * b.wall_temperature_waterwall[t, 4]*b.SR_lf[t] \
            +2.15852 * b.wall_temperature_waterwall[t, 4]*b.secondary_air_inlet.temperature[t] \
            -0.224207 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 8] \
            +0.275077 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 10] \
            -620.705 * b.wall_temperature_waterwall[t, 5]*b.SR[t] \
            -99.9039 * b.wall_temperature_waterwall[t, 6]*b.flowrate_coal_raw[t] \
            -14.5881 * b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t] \
            -534.118 * b.wall_temperature_waterwall[t, 9]*b.SR[t] \
            -1.31088 * b.wall_temperature_waterwall[t, 10]*b.wall_temperature_waterwall[t, 12] \
            -1.41959 * b.wall_temperature_waterwall[t, 11]*b.flowrate_coal_raw[t] \
            +177.07 * b.wall_temperature_waterwall[t, 11]*b.ratio_PA2coal[t] \
            +0.753534 * b.wall_temperature_platen[t]*b.flowrate_coal_raw[t] \
            +7.17536 * b.wall_temperature_roof[t]*b.flowrate_coal_raw[t] \
            +1.88795 * b.wall_temperature_roof[t]*b.secondary_air_inlet.temperature[t] \
            -83.1177 * b.wall_temperature_roof[t]*b.ratio_PA2coal[t] \
            -796480 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            -1.00571e+06 * b.flowrate_coal_raw[t]*b.SR[t] \
            +101563 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            +662.665 * b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -12441.5 * b.flowrate_coal_raw[t]*b.ratio_PA2coal[t] \
            +1.61179e+07 * b.mf_H2O_coal_raw[t]*b.SR[t] \
            -11495.5 * b.mf_H2O_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            +8.29731e+06 * b.SR[t]*b.SR_lf[t] \
            +6537.87 * b.SR[t]*b.secondary_air_inlet.temperature[t] \
            +160251 * b.SR[t]*b.ratio_PA2coal[t] \
            -9390.64 * b.SR_lf[t]*b.secondary_air_inlet.temperature[t] \
            -957.772 * b.secondary_air_inlet.temperature[t]*b.ratio_PA2coal[t] \
            +0.000819586 * (b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t])**2 \
            +0.00102671 * (b.wall_temperature_waterwall[t, 6]*b.flowrate_coal_raw[t])**2 \
            -2.36557e-06 * (b.wall_temperature_waterwall[t, 7]*b.flowrate_coal_raw[t])**2 \
            +37045.4 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**2 \
            +3479.38 * (b.flowrate_coal_raw[t]*b.SR[t])**2 \
            -0.00390226 * (b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t])**2)",
    7: "(353.147 * b.wall_temperature_waterwall[t, 1] \
            +944.788 * b.wall_temperature_waterwall[t, 2] \
            -715.67 * b.wall_temperature_waterwall[t, 3] \
            -113.773 * b.wall_temperature_waterwall[t, 4] \
            +1714.33 * b.wall_temperature_waterwall[t, 5] \
            +1186.75 * b.wall_temperature_waterwall[t, 6] \
            -41011.8 * b.wall_temperature_waterwall[t, 7] \
            +2313.62 * b.wall_temperature_waterwall[t, 8] \
            +1044.74 * b.wall_temperature_waterwall[t, 9] \
            -193.394 * b.wall_temperature_waterwall[t, 10] \
            -1258.23 * b.wall_temperature_waterwall[t, 11] \
            +971.027 * b.wall_temperature_waterwall[t, 12] \
            +1359.51 * b.wall_temperature_platen[t] \
            -694.308 * b.wall_temperature_roof[t] \
            +1.69143e+06 * b.flowrate_coal_raw[t] \
            -1.14413e+07 * b.mf_H2O_coal_raw[t] \
            +8.24903e+06 * b.SR[t] \
            -157814 * b.SR_lf[t] \
            +2768.6 * b.secondary_air_inlet.temperature[t] \
            +477829 * b.ratio_PA2coal[t] \
            +146588 * log(b.wall_temperature_waterwall[t, 1]) \
            +581939 * log(b.wall_temperature_waterwall[t, 3]) \
            -224896 * log(b.wall_temperature_waterwall[t, 4]) \
            -1.10027e+06 * log(b.wall_temperature_waterwall[t, 8]) \
            -785548 * log(b.wall_temperature_platen[t]) \
            +2.17579e+06 * log(b.flowrate_coal_raw[t]) \
            -254748 * log(b.mf_H2O_coal_raw[t]) \
            -2.27807e+07 * log(b.SR[t]) \
            +59.8201 * b.wall_temperature_waterwall[t, 7]**2 \
            -13890.7 * b.flowrate_coal_raw[t]**2 \
            -3.2383e+07 * b.mf_H2O_coal_raw[t]**2 \
            -0.0319362 * b.wall_temperature_waterwall[t, 7]**3 \
            +52.1526 * b.flowrate_coal_raw[t]**3 \
            -0.830502 * b.wall_temperature_waterwall[t, 1]*b.wall_temperature_waterwall[t, 8] \
            +0.226652 * b.wall_temperature_waterwall[t, 1]*b.wall_temperature_waterwall[t, 10] \
            -76.2281 * b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t] \
            +779.831 * b.wall_temperature_waterwall[t, 1]*b.SR_lf[t] \
            -0.348643 * b.wall_temperature_waterwall[t, 2]*b.wall_temperature_waterwall[t, 5] \
            -196.262 * b.wall_temperature_waterwall[t, 2]*b.ratio_PA2coal[t] \
            +127.532 * b.wall_temperature_waterwall[t, 3]*b.SR[t] \
            +0.125889 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_waterwall[t, 6] \
            -0.214732 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_waterwall[t, 8] \
            +0.810439 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_waterwall[t, 10] \
            -698.238 * b.wall_temperature_waterwall[t, 4]*b.SR[t] \
            +1.85785 * b.wall_temperature_waterwall[t, 4]*b.secondary_air_inlet.temperature[t] \
            -0.436732 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 8] \
            +0.407706 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 10] \
            -4.17827 * b.wall_temperature_waterwall[t, 5]*b.flowrate_coal_raw[t] \
            -703.956 * b.wall_temperature_waterwall[t, 5]*b.SR[t] \
            -0.85661 * b.wall_temperature_waterwall[t, 6]*b.wall_temperature_waterwall[t, 12] \
            -9.60464 * b.wall_temperature_waterwall[t, 6]*b.flowrate_coal_raw[t] \
            -34.2924 * b.wall_temperature_waterwall[t, 7]*b.flowrate_coal_raw[t] \
            +1107.63 * b.wall_temperature_waterwall[t, 7]*b.mf_H2O_coal_raw[t] \
            +1.26848 * b.wall_temperature_waterwall[t, 8]*b.wall_temperature_waterwall[t, 11] \
            -17.4004 * b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t] \
            -658.546 * b.wall_temperature_waterwall[t, 9]*b.SR[t] \
            -1.01075 * b.wall_temperature_waterwall[t, 10]*b.wall_temperature_waterwall[t, 12] \
            -2.06046 * b.wall_temperature_waterwall[t, 11]*b.flowrate_coal_raw[t] \
            +177.497 * b.wall_temperature_waterwall[t, 11]*b.ratio_PA2coal[t] \
            +52.1763 * b.wall_temperature_waterwall[t, 12]*b.flowrate_coal_raw[t] \
            +0.450986 * b.wall_temperature_platen[t]*b.flowrate_coal_raw[t] \
            +7.72431 * b.wall_temperature_roof[t]*b.flowrate_coal_raw[t] \
            +1.57739 * b.wall_temperature_roof[t]*b.secondary_air_inlet.temperature[t] \
            -173.703 * b.wall_temperature_roof[t]*b.ratio_PA2coal[t] \
            -718200 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            -942875 * b.flowrate_coal_raw[t]*b.SR[t] \
            +175527 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            +441.476 * b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -10706.6 * b.flowrate_coal_raw[t]*b.ratio_PA2coal[t] \
            +1.62442e+07 * b.mf_H2O_coal_raw[t]*b.SR[t] \
            -8678.04 * b.mf_H2O_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            +5.56242e+06 * b.SR[t]*b.SR_lf[t] \
            +7113.57 * b.SR[t]*b.secondary_air_inlet.temperature[t] \
            -9979.24 * b.SR_lf[t]*b.secondary_air_inlet.temperature[t] \
            -650.626 * b.secondary_air_inlet.temperature[t]*b.ratio_PA2coal[t] \
            +0.00138792 * (b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t])**2 \
            -0.00105797 * (b.wall_temperature_waterwall[t, 12]*b.flowrate_coal_raw[t])**2 \
            +33941.2 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**2 \
            +2634.71 * (b.flowrate_coal_raw[t]*b.SR[t])**2)",
    8: "(30738.3 * b.wall_temperature_waterwall[t, 1] \
            +635.315 * b.wall_temperature_waterwall[t, 2] \
            -737.628 * b.wall_temperature_waterwall[t, 3] \
            +613.096 * b.wall_temperature_waterwall[t, 4] \
            +1198.21 * b.wall_temperature_waterwall[t, 5] \
            +28980.1 * b.wall_temperature_waterwall[t, 6] \
            -86338.9 * b.wall_temperature_waterwall[t, 7] \
            -16399.5 * b.wall_temperature_waterwall[t, 8] \
            +1028.24 * b.wall_temperature_waterwall[t, 9] \
            +288.999 * b.wall_temperature_waterwall[t, 10] \
            -1919.18 * b.wall_temperature_waterwall[t, 11] \
            +1183.65 * b.wall_temperature_waterwall[t, 12] \
            +1613.28 * b.wall_temperature_platen[t] \
            -809.056 * b.wall_temperature_roof[t] \
            +1.29475e+06 * b.flowrate_coal_raw[t] \
            +2.75866e+07 * b.mf_H2O_coal_raw[t] \
            +6.67449e+06 * b.SR[t] \
            +5.28363e+07 * b.SR_lf[t] \
            -3717.86 * b.secondary_air_inlet.temperature[t] \
            +385624 * b.ratio_PA2coal[t] \
            -1.00557e+07 * log(b.wall_temperature_waterwall[t, 1]) \
            +559241 * log(b.wall_temperature_waterwall[t, 3]) \
            -1.01411e+07 * log(b.wall_temperature_waterwall[t, 6]) \
            +1.14643e+07 * log(b.wall_temperature_waterwall[t, 7]) \
            +9.46655e+06 * log(b.wall_temperature_waterwall[t, 8]) \
            -894841 * log(b.wall_temperature_platen[t]) \
            +1.54449e+06 * log(b.flowrate_coal_raw[t]) \
            -122538 * log(b.mf_H2O_coal_raw[t]) \
            -1.6811e+07 * log(b.SR[t]) \
            -5.55471e+07 * log(b.SR_lf[t]) \
            -2.83514e+06 * log(b.secondary_air_inlet.temperature[t]) \
            -4.07619e+07 * exp(b.mf_H2O_coal_raw[t]) \
            -10.6528 * b.wall_temperature_waterwall[t, 1]**2 \
            +0.272176 * b.wall_temperature_waterwall[t, 4]**2 \
            -9.34661 * b.wall_temperature_waterwall[t, 6]**2 \
            +86.226 * b.wall_temperature_waterwall[t, 7]**2 \
            -6869.16 * b.flowrate_coal_raw[t]**2 \
            -0.0341736 * b.wall_temperature_waterwall[t, 7]**3 \
            -0.908054 * b.wall_temperature_waterwall[t, 1]*b.wall_temperature_waterwall[t, 8] \
            -64.9931 * b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t] \
            -0.1575 * b.wall_temperature_waterwall[t, 2]*b.wall_temperature_waterwall[t, 5] \
            -145.333 * b.wall_temperature_waterwall[t, 2]*b.ratio_PA2coal[t] \
            +124.842 * b.wall_temperature_waterwall[t, 3]*b.SR[t] \
            +0.652536 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_waterwall[t, 10] \
            -697.56 * b.wall_temperature_waterwall[t, 4]*b.SR[t] \
            -1290.22 * b.wall_temperature_waterwall[t, 4]*b.SR_lf[t] \
            +1.54101 * b.wall_temperature_waterwall[t, 4]*b.secondary_air_inlet.temperature[t] \
            -0.340243 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 8] \
            +0.186733 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 10] \
            +0.264209 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 12] \
            -600.521 * b.wall_temperature_waterwall[t, 5]*b.SR[t] \
            -1.06736 * b.wall_temperature_waterwall[t, 6]*b.wall_temperature_waterwall[t, 12] \
            -8.99471 * b.wall_temperature_waterwall[t, 6]*b.flowrate_coal_raw[t] \
            +2.71176 * b.wall_temperature_waterwall[t, 7]*b.flowrate_coal_raw[t] \
            +0.727491 * b.wall_temperature_waterwall[t, 8]*b.wall_temperature_waterwall[t, 11] \
            -50.9128 * b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t] \
            -666.436 * b.wall_temperature_waterwall[t, 8]*b.SR[t] \
            -0.100745 * b.wall_temperature_waterwall[t, 9]*b.wall_temperature_waterwall[t, 10] \
            -530.861 * b.wall_temperature_waterwall[t, 9]*b.SR[t] \
            -0.824416 * b.wall_temperature_waterwall[t, 10]*b.wall_temperature_waterwall[t, 12] \
            -2.53957 * b.wall_temperature_waterwall[t, 11]*b.flowrate_coal_raw[t] \
            +166.069 * b.wall_temperature_waterwall[t, 11]*b.ratio_PA2coal[t] \
            +0.638522 * b.wall_temperature_platen[t]*b.flowrate_coal_raw[t] \
            +6.65548 * b.wall_temperature_roof[t]*b.flowrate_coal_raw[t] \
            +1.65814 * b.wall_temperature_roof[t]*b.secondary_air_inlet.temperature[t] \
            -132.374 * b.wall_temperature_roof[t]*b.ratio_PA2coal[t] \
            -636864 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            -703357 * b.flowrate_coal_raw[t]*b.SR[t] \
            +191874 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            +531.415 * b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -9232.96 * b.flowrate_coal_raw[t]*b.ratio_PA2coal[t] \
            +1.59172e+07 * b.mf_H2O_coal_raw[t]*b.SR[t] \
            -8594.17 * b.mf_H2O_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            +3.30615e+06 * b.SR[t]*b.SR_lf[t] \
            +6573.55 * b.SR[t]*b.secondary_air_inlet.temperature[t] \
            -555.335 * b.secondary_air_inlet.temperature[t]*b.ratio_PA2coal[t] \
            +0.00118832 * (b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t])**2 \
            +0.846319 * (b.wall_temperature_waterwall[t, 11]*b.SR_lf[t])**2 \
            +28243.5 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**2 \
            -834.918 * (b.flowrate_coal_raw[t]*b.SR[t])**2 \
            -0.00281944 * (b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t])**2 \
            +23.2033 * (b.flowrate_coal_raw[t]*b.SR[t])**3)",
    9: "(1554.78 * b.wall_temperature_waterwall[t, 1] \
            +292.513 * b.wall_temperature_waterwall[t, 2] \
            -686.842 * b.wall_temperature_waterwall[t, 3] \
            +504.056 * b.wall_temperature_waterwall[t, 4] \
            +1072 * b.wall_temperature_waterwall[t, 5] \
            +33410.5 * b.wall_temperature_waterwall[t, 6] \
            +52811.6 * b.wall_temperature_waterwall[t, 7] \
            +2682.22 * b.wall_temperature_waterwall[t, 8] \
            +59476.3 * b.wall_temperature_waterwall[t, 9] \
            +1980.99 * b.wall_temperature_waterwall[t, 10] \
            -418.882 * b.wall_temperature_waterwall[t, 11] \
            +979.077 * b.wall_temperature_waterwall[t, 12] \
            +2613.03 * b.wall_temperature_platen[t] \
            -325.642 * b.wall_temperature_roof[t] \
            +1.11231e+06 * b.flowrate_coal_raw[t] \
            +2.63584e+07 * b.mf_H2O_coal_raw[t] \
            +5.72152e+06 * b.SR[t] \
            +3.15685e+08 * b.SR_lf[t] \
            -6126.97 * b.secondary_air_inlet.temperature[t] \
            -349285 * b.ratio_PA2coal[t] \
            +417256 * log(b.wall_temperature_waterwall[t, 1]) \
            +3.53244e+06 * log(b.wall_temperature_waterwall[t, 3]) \
            -360670 * log(b.wall_temperature_waterwall[t, 4]) \
            +4.18922e+06 * log(b.wall_temperature_waterwall[t, 5]) \
            -1.17359e+07 * log(b.wall_temperature_waterwall[t, 6]) \
            -1.90854e+07 * log(b.wall_temperature_waterwall[t, 7]) \
            -1.19489e+06 * log(b.wall_temperature_waterwall[t, 8]) \
            -1.73885e+07 * log(b.wall_temperature_waterwall[t, 9]) \
            -848457 * log(b.wall_temperature_waterwall[t, 10]) \
            -1.39216e+06 * log(b.wall_temperature_platen[t]) \
            +20679.5 * log(b.wall_temperature_roof[t]) \
            +1.18208e+06 * log(b.flowrate_coal_raw[t]) \
            -138694 * log(b.mf_H2O_coal_raw[t]) \
            -3.86313e+06 * log(b.SR[t]) \
            -1.78945e+08 * log(b.SR_lf[t]) \
            -3.73556e+07 * exp(b.mf_H2O_coal_raw[t]) \
            -11.0194 * b.wall_temperature_waterwall[t, 6]**2 \
            -17.719 * b.wall_temperature_waterwall[t, 7]**2 \
            -26.8571 * b.wall_temperature_waterwall[t, 9]**2 \
            -3871.37 * b.flowrate_coal_raw[t]**2 \
            -6.88157e+07 * b.SR_lf[t]**2 \
            -24.6314 * b.flowrate_coal_raw[t]**3 \
            -0.898065 * b.wall_temperature_waterwall[t, 1]*b.wall_temperature_waterwall[t, 8] \
            -50.8152 * b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t] \
            -1.50253 * b.wall_temperature_waterwall[t, 1]*b.secondary_air_inlet.temperature[t] \
            -0.176067 * b.wall_temperature_waterwall[t, 2]*b.wall_temperature_waterwall[t, 5] \
            -0.045357 * b.wall_temperature_waterwall[t, 3]*b.wall_temperature_waterwall[t, 12] \
            -6290.31 * b.wall_temperature_waterwall[t, 3]*b.SR[t] \
            +0.0280148 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_waterwall[t, 6] \
            -0.0760966 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_platen[t] \
            -543.634 * b.wall_temperature_waterwall[t, 4]*b.SR[t] \
            +1.55669 * b.wall_temperature_waterwall[t, 4]*b.secondary_air_inlet.temperature[t] \
            -0.462604 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 8] \
            +0.552709 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 10] \
            -0.187582 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 11] \
            +0.311448 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 12] \
            -9598.22 * b.wall_temperature_waterwall[t, 5]*b.SR[t] \
            -0.893275 * b.wall_temperature_waterwall[t, 6]*b.wall_temperature_waterwall[t, 12] \
            -0.0718352 * b.wall_temperature_waterwall[t, 6]*b.wall_temperature_roof[t] \
            -7.00461 * b.wall_temperature_waterwall[t, 6]*b.flowrate_coal_raw[t] \
            -0.347099 * b.wall_temperature_waterwall[t, 7]*b.wall_temperature_waterwall[t, 10] \
            +0.986772 * b.wall_temperature_waterwall[t, 7]*b.flowrate_coal_raw[t] \
            +0.874394 * b.wall_temperature_waterwall[t, 8]*b.wall_temperature_waterwall[t, 11] \
            -13.7583 * b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t] \
            -0.111941 * b.wall_temperature_waterwall[t, 9]*b.wall_temperature_waterwall[t, 10] \
            -38.4122 * b.wall_temperature_waterwall[t, 9]*b.flowrate_coal_raw[t] \
            -926.033 * b.wall_temperature_waterwall[t, 9]*b.SR[t] \
            -0.719251 * b.wall_temperature_waterwall[t, 10]*b.wall_temperature_waterwall[t, 12] \
            +0.311257 * b.wall_temperature_platen[t]*b.flowrate_coal_raw[t] \
            +0.973601 * b.wall_temperature_roof[t]*b.secondary_air_inlet.temperature[t] \
            -86.0127 * b.wall_temperature_roof[t]*b.ratio_PA2coal[t] \
            -573877 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            -518729 * b.flowrate_coal_raw[t]*b.SR[t] \
            +177185 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            +356.318 * b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -8350.35 * b.flowrate_coal_raw[t]*b.ratio_PA2coal[t] \
            +1.40375e+07 * b.mf_H2O_coal_raw[t]*b.SR[t] \
            -7543.82 * b.mf_H2O_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            +6147.42 * b.SR[t]*b.secondary_air_inlet.temperature[t] \
            +695872 * b.SR_lf[t]*b.ratio_PA2coal[t] \
            -464.964 * b.secondary_air_inlet.temperature[t]*b.ratio_PA2coal[t] \
            +0.000944236 * (b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t])**2 \
            +1.68544 * (b.wall_temperature_waterwall[t, 3]*b.SR[t])**2 \
            +2.37726 * (b.wall_temperature_waterwall[t, 5]*b.SR[t])**2 \
            -2.11316e-05 * (b.wall_temperature_waterwall[t, 11]*b.flowrate_coal_raw[t])**2 \
            +23969.4 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**2 \
            -2523.54 * (b.flowrate_coal_raw[t]*b.SR[t])**2 \
            +31.1715 * (b.flowrate_coal_raw[t]*b.SR[t])**3)",
    10: "(1837.36 * b.wall_temperature_waterwall[t, 1] \
            +729.997 * b.wall_temperature_waterwall[t, 2] \
            -230.755 * b.wall_temperature_waterwall[t, 3] \
            +862.309 * b.wall_temperature_waterwall[t, 4] \
            -167.942 * b.wall_temperature_waterwall[t, 5] \
            +2060.72 * b.wall_temperature_waterwall[t, 6] \
            +47434.8 * b.wall_temperature_waterwall[t, 7] \
            +2552 * b.wall_temperature_waterwall[t, 8] \
            +2565.36 * b.wall_temperature_waterwall[t, 9] \
            +12188 * b.wall_temperature_waterwall[t, 10] \
            -383.468 * b.wall_temperature_waterwall[t, 11] \
            -125.387 * b.wall_temperature_waterwall[t, 12] \
            +5141.33 * b.wall_temperature_platen[t] \
            -312.616 * b.wall_temperature_roof[t] \
            +1.17132e+06 * b.flowrate_coal_raw[t] \
            +4.60141e+07 * b.mf_H2O_coal_raw[t] \
            +2.69988e+06 * b.SR[t] \
            +2.93429e+08 * b.SR_lf[t] \
            -7395.61 * b.secondary_air_inlet.temperature[t] \
            +367885 * b.ratio_PA2coal[t] \
            -112755 * log(b.wall_temperature_waterwall[t, 1]) \
            -816926 * log(b.wall_temperature_waterwall[t, 4]) \
            -627025 * log(b.wall_temperature_waterwall[t, 6]) \
            -1.72497e+07 * log(b.wall_temperature_waterwall[t, 7]) \
            -1.30207e+06 * log(b.wall_temperature_waterwall[t, 8]) \
            -833056 * log(b.wall_temperature_waterwall[t, 9]) \
            -2.73476e+06 * log(b.wall_temperature_platen[t]) \
            +790545 * log(b.flowrate_coal_raw[t]) \
            -289014 * log(b.mf_H2O_coal_raw[t]) \
            -8.88643e+06 * log(b.SR[t]) \
            -1.68062e+08 * log(b.SR_lf[t]) \
            -5.5004e+07 * exp(b.mf_H2O_coal_raw[t]) \
            -4.59464e+07 * exp(b.SR_lf[t]) \
            -15.9244 * b.wall_temperature_waterwall[t, 7]**2 \
            -12.4671 * b.wall_temperature_waterwall[t, 10]**2 \
            -0.0548467 * b.wall_temperature_platen[t]**2 \
            -3084.74 * b.flowrate_coal_raw[t]**2 \
            -20.7712 * b.flowrate_coal_raw[t]**3 \
            -3.52178 * b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t] \
            -2.37475 * b.wall_temperature_waterwall[t, 1]*b.secondary_air_inlet.temperature[t] \
            -0.208333 * b.wall_temperature_waterwall[t, 2]*b.wall_temperature_waterwall[t, 5] \
            -171.455 * b.wall_temperature_waterwall[t, 2]*b.ratio_PA2coal[t] \
            +306.026 * b.wall_temperature_waterwall[t, 3]*b.SR[t] \
            -0.15677 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_waterwall[t, 6] \
            +1.33655 * b.wall_temperature_waterwall[t, 4]*b.secondary_air_inlet.temperature[t] \
            +0.485426 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 10] \
            +0.374308 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 12] \
            -0.828451 * b.wall_temperature_waterwall[t, 6]*b.wall_temperature_waterwall[t, 12] \
            -0.180578 * b.wall_temperature_waterwall[t, 6]*b.wall_temperature_roof[t] \
            -13.9066 * b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t] \
            -671.924 * b.wall_temperature_waterwall[t, 9]*b.SR[t] \
            +0.538709 * b.wall_temperature_waterwall[t, 10]*b.wall_temperature_waterwall[t, 11] \
            -56.0447 * b.wall_temperature_waterwall[t, 10]*b.flowrate_coal_raw[t] \
            -712.621 * b.wall_temperature_waterwall[t, 10]*b.SR[t] \
            +44.9712 * b.wall_temperature_waterwall[t, 12]*b.flowrate_coal_raw[t] \
            +75.2145 * b.wall_temperature_waterwall[t, 12]*b.ratio_PA2coal[t] \
            -5.25712 * b.wall_temperature_platen[t]*b.flowrate_coal_raw[t] \
            +0.836281 * b.wall_temperature_roof[t]*b.secondary_air_inlet.temperature[t] \
            -693841 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            -489211 * b.flowrate_coal_raw[t]*b.SR[t] \
            +199410 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            +409.786 * b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -9506.34 * b.flowrate_coal_raw[t]*b.ratio_PA2coal[t] \
            +1.57634e+07 * b.mf_H2O_coal_raw[t]*b.SR[t] \
            -9269.71 * b.mf_H2O_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            +7874.35 * b.SR[t]*b.secondary_air_inlet.temperature[t] \
            -575.815 * b.secondary_air_inlet.temperature[t]*b.ratio_PA2coal[t] \
            +0.019545 * (b.wall_temperature_waterwall[t, 11]*b.ratio_PA2coal[t])**2 \
            -0.000956923 * (b.wall_temperature_waterwall[t, 12]*b.flowrate_coal_raw[t])**2 \
            +26897.5 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**2 \
            -3926.81 * (b.flowrate_coal_raw[t]*b.SR[t])**2 \
            -6.04217e-05 * (b.wall_temperature_waterwall[t, 4]*b.SR[t])**3 \
            +39.5236 * (b.flowrate_coal_raw[t]*b.SR[t])**3)",
    11: "(265.482 * b.wall_temperature_waterwall[t, 1] \
            +24.2531 * b.wall_temperature_waterwall[t, 2] \
            -632.925 * b.wall_temperature_waterwall[t, 3] \
            +114.247 * b.wall_temperature_waterwall[t, 4] \
            +28.3431 * b.wall_temperature_waterwall[t, 5] \
            +445.719 * b.wall_temperature_waterwall[t, 6] \
            +1363.91 * b.wall_temperature_waterwall[t, 7] \
            +595.684 * b.wall_temperature_waterwall[t, 8] \
            +1378.1 * b.wall_temperature_waterwall[t, 9] \
            +1082.57 * b.wall_temperature_waterwall[t, 10] \
            +29176.4 * b.wall_temperature_waterwall[t, 11] \
            +385.741 * b.wall_temperature_waterwall[t, 12] \
            -9989.8 * b.wall_temperature_platen[t] \
            -387.904 * b.wall_temperature_roof[t] \
            +459422 * b.flowrate_coal_raw[t] \
            -8.16755e+06 * b.mf_H2O_coal_raw[t] \
            -2.7349e+06 * b.SR[t] \
            +2.68906e+07 * b.SR_lf[t] \
            -3539.12 * b.secondary_air_inlet.temperature[t] \
            +341694 * b.ratio_PA2coal[t] \
            +143448 * log(b.wall_temperature_waterwall[t, 1]) \
            +137163 * log(b.wall_temperature_waterwall[t, 3]) \
            -174890 * log(b.wall_temperature_waterwall[t, 4]) \
            -283768 * log(b.wall_temperature_waterwall[t, 6]) \
            -313430 * log(b.wall_temperature_waterwall[t, 7]) \
            -421769 * log(b.wall_temperature_waterwall[t, 8]) \
            -498950 * log(b.wall_temperature_waterwall[t, 9]) \
            -628469 * log(b.wall_temperature_waterwall[t, 10]) \
            -4.79246e+06 * log(b.wall_temperature_waterwall[t, 11]) \
            +2.7541e+06 * log(b.wall_temperature_platen[t]) \
            -30902.5 * log(b.flowrate_coal_raw[t]) \
            +423000 * log(b.SR[t]) \
            -2.62773e+07 * log(b.SR_lf[t]) \
            -856849 * log(b.secondary_air_inlet.temperature[t]) \
            -976839 * exp(b.SR_lf[t]) \
            -18.8883 * b.wall_temperature_waterwall[t, 11]**2 \
            -397.68 * b.flowrate_coal_raw[t]**2 \
            +0.00509064 * b.wall_temperature_platen[t]**3 \
            +0.000109843 * b.wall_temperature_roof[t]**3 \
            -13.0398 * b.flowrate_coal_raw[t]**3 \
            -7.10244e+06 * b.mf_H2O_coal_raw[t]**3 \
            -0.475103 * b.wall_temperature_waterwall[t, 1]*b.wall_temperature_waterwall[t, 12] \
            -1.28357 * b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t] \
            +86.2882 * b.wall_temperature_waterwall[t, 1]*b.SR[t] \
            +446.141 * b.wall_temperature_waterwall[t, 1]*b.SR_lf[t] \
            -0.996705 * b.wall_temperature_waterwall[t, 1]*b.secondary_air_inlet.temperature[t] \
            +0.0280906 * b.wall_temperature_waterwall[t, 2]*b.wall_temperature_waterwall[t, 5] \
            +0.232787 * b.wall_temperature_waterwall[t, 2]*b.secondary_air_inlet.temperature[t] \
            -48.0401 * b.wall_temperature_waterwall[t, 2]*b.ratio_PA2coal[t] \
            +162.375 * b.wall_temperature_waterwall[t, 3]*b.SR[t] \
            +297.717 * b.wall_temperature_waterwall[t, 3]*b.SR_lf[t] \
            +0.232195 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_roof[t] \
            +2.02509 * b.wall_temperature_waterwall[t, 4]*b.flowrate_coal_raw[t] \
            +0.0771506 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 12] \
            -2.17805 * b.wall_temperature_waterwall[t, 6]*b.flowrate_coal_raw[t] \
            +2186.86 * b.wall_temperature_waterwall[t, 6]*b.mf_H2O_coal_raw[t] \
            -2.97757 * b.wall_temperature_waterwall[t, 6]*b.ratio_PA2coal[t] \
            -0.0238492 * b.wall_temperature_waterwall[t, 7]*b.wall_temperature_waterwall[t, 11] \
            +226.047 * b.wall_temperature_waterwall[t, 7]*b.SR[t] \
            -1085.84 * b.wall_temperature_waterwall[t, 7]*b.SR_lf[t] \
            +0.236826 * b.wall_temperature_waterwall[t, 8]*b.wall_temperature_waterwall[t, 10] \
            -2.73927 * b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t] \
            +539.145 * b.wall_temperature_waterwall[t, 8]*b.mf_H2O_coal_raw[t] \
            -0.104024 * b.wall_temperature_waterwall[t, 9]*b.wall_temperature_roof[t] \
            -291.787 * b.wall_temperature_waterwall[t, 9]*b.SR[t] \
            +13.4721 * b.wall_temperature_waterwall[t, 10]*b.flowrate_coal_raw[t] \
            -124.929 * b.wall_temperature_waterwall[t, 10]*b.SR[t] \
            -15.757 * b.wall_temperature_waterwall[t, 11]*b.flowrate_coal_raw[t] \
            -4817.97 * b.wall_temperature_waterwall[t, 11]*b.SR[t] \
            +1892.85 * b.wall_temperature_waterwall[t, 11]*b.SR_lf[t] \
            +1.63258 * b.wall_temperature_platen[t]*b.flowrate_coal_raw[t] \
            +2676.98 * b.wall_temperature_platen[t]*b.SR[t] \
            +0.49742 * b.wall_temperature_platen[t]*b.secondary_air_inlet.temperature[t] \
            -72.8604 * b.wall_temperature_platen[t]*b.ratio_PA2coal[t] \
            +0.430275 * b.wall_temperature_roof[t]*b.secondary_air_inlet.temperature[t] \
            -310586 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            -142967 * b.flowrate_coal_raw[t]*b.SR[t] \
            +88449.6 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            +217.961 * b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -7117.74 * b.flowrate_coal_raw[t]*b.ratio_PA2coal[t] \
            +6.1753e+06 * b.mf_H2O_coal_raw[t]*b.SR[t] \
            -3411.28 * b.mf_H2O_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -106414 * b.mf_H2O_coal_raw[t]*b.ratio_PA2coal[t] \
            +845856 * b.SR[t]*b.SR_lf[t] \
            +4202.28 * b.SR[t]*b.secondary_air_inlet.temperature[t] \
            -340.792 * b.secondary_air_inlet.temperature[t]*b.ratio_PA2coal[t] \
            -5.1335e-05 * (b.wall_temperature_waterwall[t, 9]*b.flowrate_coal_raw[t])**2 \
            -0.000310581 * (b.wall_temperature_waterwall[t, 10]*b.flowrate_coal_raw[t])**2 \
            +1.18396 * (b.wall_temperature_waterwall[t, 11]*b.SR[t])**2 \
            -0.716778 * (b.wall_temperature_platen[t]*b.SR[t])**2 \
            +11739 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**2 \
            -2422.02 * (b.flowrate_coal_raw[t]*b.SR[t])**2 \
            -0.000797508 * (b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t])**2 \
            +16.1006 * (b.flowrate_coal_raw[t]*b.ratio_PA2coal[t])**2 \
            -0.0417793 * (b.wall_temperature_waterwall[t, 6]*b.mf_H2O_coal_raw[t])**3 \
            +19.1183 * (b.flowrate_coal_raw[t]*b.SR[t])**3)",
    12: "(220.386 * b.wall_temperature_waterwall[t, 1] \
            +67.4508 * b.wall_temperature_waterwall[t, 2] \
            -53.5201 * b.wall_temperature_waterwall[t, 3] \
            +168.186 * b.wall_temperature_waterwall[t, 4] \
            +164.298 * b.wall_temperature_waterwall[t, 5] \
            +274.74 * b.wall_temperature_waterwall[t, 6] \
            +159.494 * b.wall_temperature_waterwall[t, 7] \
            +274.631 * b.wall_temperature_waterwall[t, 8] \
            +589.546 * b.wall_temperature_waterwall[t, 9] \
            +692.692 * b.wall_temperature_waterwall[t, 10] \
            -125.502 * b.wall_temperature_waterwall[t, 11] \
            +43494.6 * b.wall_temperature_waterwall[t, 12] \
            -33880.3 * b.wall_temperature_platen[t] \
            -9863.55 * b.wall_temperature_roof[t] \
            +241200 * b.flowrate_coal_raw[t] \
            +2.77276e+07 * b.mf_H2O_coal_raw[t] \
            -2.88054e+06 * b.SR[t] \
            +1.06869e+07 * b.SR_lf[t] \
            -14610.7 * b.secondary_air_inlet.temperature[t] \
            +172036 * b.ratio_PA2coal[t] \
            +47797.2 * log(b.wall_temperature_waterwall[t, 1]) \
            -11594.2 * log(b.wall_temperature_waterwall[t, 3]) \
            -125501 * log(b.wall_temperature_waterwall[t, 4]) \
            -144576 * log(b.wall_temperature_waterwall[t, 6]) \
            -198449 * log(b.wall_temperature_waterwall[t, 7]) \
            -221594 * log(b.wall_temperature_waterwall[t, 8]) \
            -199032 * log(b.wall_temperature_waterwall[t, 9]) \
            -258237 * log(b.wall_temperature_waterwall[t, 10]) \
            -387498 * log(b.wall_temperature_waterwall[t, 11]) \
            -1.13543e+07 * log(b.wall_temperature_waterwall[t, 12]) \
            +8.627e+06 * log(b.wall_temperature_platen[t]) \
            +3.05064e+06 * log(b.wall_temperature_roof[t]) \
            -246256 * log(b.flowrate_coal_raw[t]) \
            -79998 * log(b.mf_H2O_coal_raw[t]) \
            +2.16333e+06 * log(b.SR[t]) \
            -1.1065e+07 * log(b.SR_lf[t]) \
            +4.91018e+06 * log(b.secondary_air_inlet.temperature[t]) \
            -2.68698e+07 * exp(b.mf_H2O_coal_raw[t]) \
            -23.3046 * b.wall_temperature_waterwall[t, 12]**2 \
            +16.9006 * b.wall_temperature_platen[t]**2 \
            +4.26692 * b.wall_temperature_roof[t]**2 \
            +1231.98 * b.flowrate_coal_raw[t]**2 \
            +4.60765 * b.secondary_air_inlet.temperature[t]**2 \
            -20.1062 * b.flowrate_coal_raw[t]**3 \
            -0.0208005 * b.wall_temperature_waterwall[t, 1]*b.wall_temperature_platen[t] \
            +0.211335 * b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t] \
            +142.559 * b.wall_temperature_waterwall[t, 1]*b.SR[t] \
            -0.708488 * b.wall_temperature_waterwall[t, 1]*b.secondary_air_inlet.temperature[t] \
            +0.104243 * b.wall_temperature_waterwall[t, 2]*b.wall_temperature_waterwall[t, 5] \
            -43.2466 * b.wall_temperature_waterwall[t, 2]*b.ratio_PA2coal[t] \
            +84.0522 * b.wall_temperature_waterwall[t, 3]*b.SR[t] \
            -0.097681 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_waterwall[t, 6] \
            +0.168211 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_roof[t] \
            +0.0311415 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_platen[t] \
            -201.961 * b.wall_temperature_waterwall[t, 5]*b.SR_lf[t] \
            +0.0979395 * b.wall_temperature_waterwall[t, 6]*b.wall_temperature_platen[t] \
            +152.533 * b.wall_temperature_waterwall[t, 7]*b.SR[t] \
            +0.138216 * b.wall_temperature_waterwall[t, 8]*b.wall_temperature_waterwall[t, 12] \
            +358.824 * b.wall_temperature_waterwall[t, 8]*b.mf_H2O_coal_raw[t] \
            -1.73455 * b.wall_temperature_waterwall[t, 9]*b.flowrate_coal_raw[t] \
            -132.667 * b.wall_temperature_waterwall[t, 9]*b.SR[t] \
            -0.128919 * b.wall_temperature_waterwall[t, 10]*b.wall_temperature_roof[t] \
            -1.61439 * b.wall_temperature_waterwall[t, 10]*b.flowrate_coal_raw[t] \
            -51.0787 * b.wall_temperature_waterwall[t, 10]*b.SR[t] \
            +0.415955 * b.wall_temperature_waterwall[t, 11]*b.flowrate_coal_raw[t] \
            -109.65 * b.wall_temperature_waterwall[t, 11]*b.SR[t] \
            +976.739 * b.wall_temperature_waterwall[t, 11]*b.SR_lf[t] \
            -3.17061 * b.wall_temperature_waterwall[t, 12]*b.flowrate_coal_raw[t] \
            -220.696 * b.wall_temperature_waterwall[t, 12]*b.SR[t] \
            +31.6233 * b.wall_temperature_waterwall[t, 12]*b.ratio_PA2coal[t] \
            +3.78535 * b.wall_temperature_platen[t]*b.flowrate_coal_raw[t] \
            +2053.18 * b.wall_temperature_platen[t]*b.SR[t] \
            +0.387234 * b.wall_temperature_platen[t]*b.secondary_air_inlet.temperature[t] \
            -46.3757 * b.wall_temperature_platen[t]*b.ratio_PA2coal[t] \
            +0.097768 * b.wall_temperature_roof[t]*b.secondary_air_inlet.temperature[t] \
            +2.52659 * b.wall_temperature_roof[t]*b.ratio_PA2coal[t] \
            -249015 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            -15466 * b.flowrate_coal_raw[t]*b.SR[t] \
            +58532.1 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            +121.404 * b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -3138.31 * b.flowrate_coal_raw[t]*b.ratio_PA2coal[t] \
            +2.05815e+06 * b.mf_H2O_coal_raw[t]*b.SR[t] \
            -164857 * b.mf_H2O_coal_raw[t]*b.SR_lf[t] \
            -2532.86 * b.mf_H2O_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -861.46 * b.SR[t]*b.secondary_air_inlet.temperature[t] \
            -209.784 * b.secondary_air_inlet.temperature[t]*b.ratio_PA2coal[t] \
            +0.00576341 * (b.wall_temperature_waterwall[t, 11]*b.ratio_PA2coal[t])**2 \
            -0.000253398 * (b.wall_temperature_waterwall[t, 12]*b.flowrate_coal_raw[t])**2 \
            -0.592335 * (b.wall_temperature_platen[t]*b.SR[t])**2 \
            +9328.86 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**2 \
            -2333.23 * (b.flowrate_coal_raw[t]*b.SR[t])**2 \
            +3.15989e+06 * (b.mf_H2O_coal_raw[t]*b.SR[t])**2 \
            +1.25506 * (b.SR[t]*b.secondary_air_inlet.temperature[t])**2 \
            -0.000219682 * (b.wall_temperature_waterwall[t, 12]*b.SR_lf[t])**3 \
            +14.8923 * (b.flowrate_coal_raw[t]*b.SR[t])**3)",
    "pl": "(411941 * b.wall_temperature_waterwall[t, 1] \
            +2599.73 * b.wall_temperature_waterwall[t, 2] \
            -946.452 * b.wall_temperature_waterwall[t, 3] \
            +3366.98 * b.wall_temperature_waterwall[t, 4] \
            +208.286 * b.wall_temperature_waterwall[t, 5] \
            +9022.67 * b.wall_temperature_waterwall[t, 6] \
            +1152.84 * b.wall_temperature_waterwall[t, 7] \
            +5576.24 * b.wall_temperature_waterwall[t, 8] \
            +9303.39 * b.wall_temperature_waterwall[t, 9] \
            +8901.32 * b.wall_temperature_waterwall[t, 10] \
            +2840.75 * b.wall_temperature_waterwall[t, 11] \
            -110147 * b.wall_temperature_waterwall[t, 12] \
            +259512 * b.wall_temperature_platen[t] \
            -54879.5 * b.wall_temperature_roof[t] \
            +4.15488e+06 * b.flowrate_coal_raw[t] \
            -3.15854e+07 * b.mf_H2O_coal_raw[t] \
            -2.10085e+07 * b.SR[t] \
            +1.70259e+08 * b.SR_lf[t] \
            -34982.2 * b.secondary_air_inlet.temperature[t] \
            +1.74171e+06 * b.ratio_PA2coal[t] \
            -1.05527e+08 * log(b.wall_temperature_waterwall[t, 1]) \
            -9.8091e+06 * log(b.wall_temperature_waterwall[t, 4]) \
            -2.12136e+06 * log(b.wall_temperature_waterwall[t, 6]) \
            -1.71391e+06 * log(b.wall_temperature_waterwall[t, 7]) \
            -2.30845e+06 * log(b.wall_temperature_waterwall[t, 8]) \
            -2.87745e+06 * log(b.wall_temperature_waterwall[t, 9]) \
            -4.14628e+06 * log(b.wall_temperature_waterwall[t, 10]) \
            +2.58819e+06 * log(b.wall_temperature_waterwall[t, 11]) \
            +3.70654e+07 * log(b.wall_temperature_waterwall[t, 12]) \
            +1.33648e+07 * log(b.wall_temperature_platen[t]) \
            +1.645e+07 * log(b.wall_temperature_roof[t]) \
            -2.78402e+06 * log(b.flowrate_coal_raw[t]) \
            +1.22061e+08 * log(b.SR[t]) \
            -1.80973e+08 * log(b.SR_lf[t]) \
            +2.24892e+07 * log(b.secondary_air_inlet.temperature[t]) \
            -256.702 * b.wall_temperature_waterwall[t, 1]**2 \
            +44.0874 * b.wall_temperature_waterwall[t, 12]**2 \
            -127.481 * b.wall_temperature_platen[t]**2 \
            +22.8917 * b.wall_temperature_roof[t]**2 \
            -1180.45 * b.flowrate_coal_raw[t]**2 \
            +0.0704504 * b.wall_temperature_waterwall[t, 1]**3 \
            -95.2406 * b.flowrate_coal_raw[t]**3 \
            -3.57179e+08 * b.mf_H2O_coal_raw[t]**3 \
            -1.32576 * b.wall_temperature_waterwall[t, 1]*b.wall_temperature_waterwall[t, 8] \
            +0.503933 * b.wall_temperature_waterwall[t, 1]*b.wall_temperature_waterwall[t, 10] \
            -2.70951 * b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t] \
            -6.71086 * b.wall_temperature_waterwall[t, 1]*b.secondary_air_inlet.temperature[t] \
            -2.22919 * b.wall_temperature_waterwall[t, 2]*b.wall_temperature_waterwall[t, 4] \
            +0.766666 * b.wall_temperature_waterwall[t, 2]*b.wall_temperature_waterwall[t, 5] \
            -441.236 * b.wall_temperature_waterwall[t, 2]*b.ratio_PA2coal[t] \
            +1018.32 * b.wall_temperature_waterwall[t, 3]*b.SR[t] \
            +0.631224 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_platen[t] \
            +1.77501 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_roof[t] \
            +12571.3 * b.wall_temperature_waterwall[t, 4]*b.SR[t] \
            +0.0318959 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 12] \
            -2.82968 * b.wall_temperature_waterwall[t, 6]*b.wall_temperature_waterwall[t, 12] \
            -5.14737 * b.wall_temperature_waterwall[t, 6]*b.secondary_air_inlet.temperature[t] \
            +10.3483 * b.wall_temperature_waterwall[t, 7]*b.flowrate_coal_raw[t] \
            +1588.54 * b.wall_temperature_waterwall[t, 7]*b.SR[t] \
            -24.8692 * b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t] \
            +4266.37 * b.wall_temperature_waterwall[t, 8]*b.mf_H2O_coal_raw[t] \
            -1.72834 * b.wall_temperature_waterwall[t, 9]*b.wall_temperature_roof[t] \
            -20.1069 * b.wall_temperature_waterwall[t, 9]*b.flowrate_coal_raw[t] \
            -1743.83 * b.wall_temperature_waterwall[t, 9]*b.SR[t] \
            +72.025 * b.wall_temperature_waterwall[t, 9]*b.ratio_PA2coal[t] \
            -21.9304 * b.wall_temperature_waterwall[t, 10]*b.flowrate_coal_raw[t] \
            -604.996 * b.wall_temperature_waterwall[t, 10]*b.SR[t] \
            -0.803033 * b.wall_temperature_waterwall[t, 11]*b.flowrate_coal_raw[t] \
            -25055.4 * b.wall_temperature_waterwall[t, 11]*b.SR[t] \
            +15039.3 * b.wall_temperature_waterwall[t, 11]*b.SR_lf[t] \
            -6.70118 * b.wall_temperature_waterwall[t, 12]*b.flowrate_coal_raw[t] \
            +259.828 * b.wall_temperature_waterwall[t, 12]*b.ratio_PA2coal[t] \
            -882.85 * b.wall_temperature_platen[t]*b.flowrate_coal_raw[t] \
            +13754.4 * b.wall_temperature_platen[t]*b.mf_H2O_coal_raw[t] \
            -264793 * b.wall_temperature_platen[t]*b.SR[t] \
            +4.22028 * b.wall_temperature_platen[t]*b.secondary_air_inlet.temperature[t] \
            -466.574 * b.wall_temperature_platen[t]*b.ratio_PA2coal[t] \
            +1.01125 * b.wall_temperature_roof[t]*b.flowrate_coal_raw[t] \
            +4891.17 * b.wall_temperature_roof[t]*b.mf_H2O_coal_raw[t] \
            +2.58075 * b.wall_temperature_roof[t]*b.secondary_air_inlet.temperature[t] \
            -2.74627e+06 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            -600924 * b.flowrate_coal_raw[t]*b.SR[t] \
            +688577 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            +1399.54 * b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -32498.3 * b.flowrate_coal_raw[t]*b.ratio_PA2coal[t] \
            +1.81724e+07 * b.mf_H2O_coal_raw[t]*b.SR[t] \
            -25011.8 * b.mf_H2O_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -26275 * b.SR[t]*b.secondary_air_inlet.temperature[t] \
            -2285.29 * b.secondary_air_inlet.temperature[t]*b.ratio_PA2coal[t] \
            +6.03477 * (b.wall_temperature_waterwall[t, 11]*b.SR[t])**2 \
            +0.0812934 * (b.wall_temperature_waterwall[t, 11]*b.ratio_PA2coal[t])**2 \
            +0.00133992 * (b.wall_temperature_platen[t]*b.flowrate_coal_raw[t])**2 \
            +131.081 * (b.wall_temperature_platen[t]*b.SR[t])**2 \
            +74778.6 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**2 \
            -18994.8 * (b.flowrate_coal_raw[t]*b.SR[t])**2 \
            +3.50946e+07 * (b.mf_H2O_coal_raw[t]*b.SR[t])**2 \
            +19.4424 * (b.SR[t]*b.secondary_air_inlet.temperature[t])**2 \
            -0.0016092 * (b.wall_temperature_waterwall[t, 4]*b.SR[t])**3 \
            -0.0295869 * (b.wall_temperature_platen[t]*b.SR[t])**3 \
            +127.257 * (b.flowrate_coal_raw[t]*b.SR[t])**3)",
    "roof": "(279.354 * b.wall_temperature_waterwall[t, 1] \
            +142.292 * b.wall_temperature_waterwall[t, 2] \
            -82.9421 * b.wall_temperature_waterwall[t, 3] \
            +273.335 * b.wall_temperature_waterwall[t, 4] \
            +892.852 * b.wall_temperature_waterwall[t, 5] \
            +348.631 * b.wall_temperature_waterwall[t, 6] \
            +1436.58 * b.wall_temperature_waterwall[t, 7] \
            +905.29 * b.wall_temperature_waterwall[t, 8] \
            +979.254 * b.wall_temperature_waterwall[t, 9] \
            +1221.17 * b.wall_temperature_waterwall[t, 10] \
            -688.783 * b.wall_temperature_waterwall[t, 11] \
            +2284.85 * b.wall_temperature_waterwall[t, 12] \
            -40634 * b.wall_temperature_platen[t] \
            +34289.7 * b.wall_temperature_roof[t] \
            +358502 * b.flowrate_coal_raw[t] \
            +4.78726e+07 * b.mf_H2O_coal_raw[t] \
            +4.59178e+07 * b.SR[t] \
            +2.01021e+07 * b.SR_lf[t] \
            -14947.3 * b.secondary_air_inlet.temperature[t] \
            +317109 * b.ratio_PA2coal[t] \
            +100628 * log(b.wall_temperature_waterwall[t, 1]) \
            -259749 * log(b.wall_temperature_waterwall[t, 4]) \
            -246590 * log(b.wall_temperature_waterwall[t, 6]) \
            -442193 * log(b.wall_temperature_waterwall[t, 7]) \
            -431568 * log(b.wall_temperature_waterwall[t, 8]) \
            -441819 * log(b.wall_temperature_waterwall[t, 9]) \
            -668967 * log(b.wall_temperature_waterwall[t, 10]) \
            +1.31422e+06 * log(b.wall_temperature_waterwall[t, 11]) \
            -1.2249e+06 * log(b.wall_temperature_waterwall[t, 12]) \
            +1.01753e+07 * log(b.wall_temperature_platen[t]) \
            -7.76406e+06 * log(b.wall_temperature_roof[t]) \
            -168868 * log(b.flowrate_coal_raw[t]) \
            -155788 * log(b.mf_H2O_coal_raw[t]) \
            +3.56858e+06 * log(b.SR[t]) \
            +1.09389e+07 * log(b.SR_lf[t]) \
            +5.70227e+06 * log(b.secondary_air_inlet.temperature[t]) \
            -4.72384e+07 * exp(b.mf_H2O_coal_raw[t]) \
            +20.1487 * b.wall_temperature_platen[t]**2 \
            -20.9093 * b.wall_temperature_roof[t]**2 \
            +2180.83 * b.flowrate_coal_raw[t]**2 \
            -9.65879e+06 * b.SR[t]**2 \
            +4.36136 * b.secondary_air_inlet.temperature[t]**2 \
            -33.3553 * b.flowrate_coal_raw[t]**3 \
            -0.290778 * b.wall_temperature_waterwall[t, 1]*b.flowrate_coal_raw[t] \
            +221.494 * b.wall_temperature_waterwall[t, 1]*b.SR[t] \
            -1.07102 * b.wall_temperature_waterwall[t, 1]*b.secondary_air_inlet.temperature[t] \
            +0.0884997 * b.wall_temperature_waterwall[t, 2]*b.wall_temperature_waterwall[t, 5] \
            -56.5688 * b.wall_temperature_waterwall[t, 2]*b.ratio_PA2coal[t] \
            +112.139 * b.wall_temperature_waterwall[t, 3]*b.SR[t] \
            +0.228812 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_roof[t] \
            -0.233352 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 8] \
            -689.204 * b.wall_temperature_waterwall[t, 5]*b.SR_lf[t] \
            +0.183533 * b.wall_temperature_waterwall[t, 6]*b.wall_temperature_platen[t] \
            +2.40284 * b.wall_temperature_waterwall[t, 7]*b.flowrate_coal_raw[t] \
            +324.697 * b.wall_temperature_waterwall[t, 7]*b.SR[t] \
            -1168.41 * b.wall_temperature_waterwall[t, 7]*b.SR_lf[t] \
            +492.766 * b.wall_temperature_waterwall[t, 8]*b.mf_H2O_coal_raw[t] \
            -2.8134 * b.wall_temperature_waterwall[t, 9]*b.flowrate_coal_raw[t] \
            +0.256616 * b.wall_temperature_waterwall[t, 10]*b.wall_temperature_waterwall[t, 11] \
            -0.25843 * b.wall_temperature_waterwall[t, 10]*b.wall_temperature_roof[t] \
            +19.2408 * b.wall_temperature_waterwall[t, 10]*b.flowrate_coal_raw[t] \
            -138.355 * b.wall_temperature_waterwall[t, 10]*b.SR[t] \
            -4483.2 * b.wall_temperature_waterwall[t, 11]*b.SR[t] \
            +1985.73 * b.wall_temperature_waterwall[t, 11]*b.SR_lf[t] \
            -3.09845 * b.wall_temperature_waterwall[t, 12]*b.flowrate_coal_raw[t] \
            +60.349 * b.wall_temperature_waterwall[t, 12]*b.ratio_PA2coal[t] \
            +3.7545 * b.wall_temperature_platen[t]*b.flowrate_coal_raw[t] \
            +3376.5 * b.wall_temperature_platen[t]*b.SR[t] \
            -81.1028 * b.wall_temperature_platen[t]*b.ratio_PA2coal[t] \
            -10.9151 * b.wall_temperature_roof[t]*b.flowrate_coal_raw[t] \
            +638.803 * b.wall_temperature_roof[t]*b.mf_H2O_coal_raw[t] \
            +0.212785 * b.wall_temperature_roof[t]*b.secondary_air_inlet.temperature[t] \
            -388386 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            -55185.6 * b.flowrate_coal_raw[t]*b.SR[t] \
            +84891 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            +226.151 * b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -8510.32 * b.flowrate_coal_raw[t]*b.ratio_PA2coal[t] \
            +3.64716e+06 * b.mf_H2O_coal_raw[t]*b.SR[t] \
            -4313.24 * b.mf_H2O_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            -4.91759e+07 * b.SR[t]*b.SR_lf[t] \
            -2591.82 * b.SR[t]*b.secondary_air_inlet.temperature[t] \
            -343.32 * b.secondary_air_inlet.temperature[t]*b.ratio_PA2coal[t] \
            -0.0526747 * (b.wall_temperature_waterwall[t, 9]*b.SR[t])**2 \
            -0.000385934 * (b.wall_temperature_waterwall[t, 10]*b.flowrate_coal_raw[t])**2 \
            +1.16051 * (b.wall_temperature_waterwall[t, 11]*b.SR[t])**2 \
            -0.947626 * (b.wall_temperature_platen[t]*b.SR[t])**2 \
            +13089.9 * (b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t])**2 \
            -3406.54 * (b.flowrate_coal_raw[t]*b.SR[t])**2 \
            -0.000648143 * (b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t])**2 \
            +20.2802 * (b.flowrate_coal_raw[t]*b.ratio_PA2coal[t])**2 \
            +5.86754e+06 * (b.mf_H2O_coal_raw[t]*b.SR[t])**2 \
            +9.71079e+06 * (b.SR[t]*b.SR_lf[t])**2 \
            +2.31621 * (b.SR[t]*b.secondary_air_inlet.temperature[t])**2 \
            +22.1927 * (b.flowrate_coal_raw[t]*b.SR[t])**3)",
    "flyash": "(exp(7.78102e-05 * b.wall_temperature_waterwall[t, 1] \
            -7.54006e-05 * b.wall_temperature_waterwall[t, 2] \
            -3.89992e-05 * b.wall_temperature_waterwall[t, 3] \
            +0.000219719 * b.wall_temperature_waterwall[t, 4] \
            -3.75494e-05 * b.wall_temperature_waterwall[t, 5] \
            -0.000963424 * b.wall_temperature_waterwall[t, 6] \
            -4.89079e-05 * b.wall_temperature_waterwall[t, 7] \
            +0.000204467 * b.wall_temperature_waterwall[t, 8] \
            -0.000143756 * b.wall_temperature_waterwall[t, 10] \
            -0.000389332 * b.wall_temperature_waterwall[t, 11] \
            -0.000338076 * b.wall_temperature_platen[t] \
            +0.241386 * b.flowrate_coal_raw[t] \
            +2.67141 * b.mf_H2O_coal_raw[t] \
            +910.531 * b.SR[t] \
            +115.082 * b.SR_lf[t] \
            +0.00275081 * b.secondary_air_inlet.temperature[t] \
            -0.10997 * b.ratio_PA2coal[t] \
            -2.10237 * log(b.flowrate_coal_raw[t]) \
            -535.077 * log(b.SR[t]) \
            -36.5477 * log(b.SR_lf[t]) \
            +90.074 * exp(b.SR[t]) \
            -667.684 * exp(b.SR_lf[t]) \
            -4.24234e-08 * b.wall_temperature_roof[t]**2 \
            -0.00369243 * b.flowrate_coal_raw[t]**2 \
            -317.41 * b.SR[t]**2 \
            +865.673 * b.SR_lf[t]**2 \
            +9.64582e-06 * b.flowrate_coal_raw[t]**3 \
            -1.69795e-07 * b.wall_temperature_waterwall[t, 1]*b.wall_temperature_waterwall[t, 9] \
            -2.65157e-07 * b.wall_temperature_waterwall[t, 2]*b.wall_temperature_waterwall[t, 5] \
            +0.000163992 * b.wall_temperature_waterwall[t, 2]*b.SR[t] \
            -2.36504e-07 * b.wall_temperature_waterwall[t, 4]*b.wall_temperature_waterwall[t, 8] \
            +5.23042e-06 * b.wall_temperature_waterwall[t, 4]*b.flowrate_coal_raw[t] \
            -0.000165057 * b.wall_temperature_waterwall[t, 4]*b.SR[t] \
            +5.85872e-08 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_waterwall[t, 11] \
            +4.41736e-06 * b.wall_temperature_waterwall[t, 5]*b.flowrate_coal_raw[t] \
            +4.67193e-07 * b.wall_temperature_waterwall[t, 6]*b.wall_temperature_waterwall[t, 11] \
            +7.09971e-06 * b.wall_temperature_waterwall[t, 6]*b.flowrate_coal_raw[t] \
            +6.36398e-07 * b.wall_temperature_waterwall[t, 6]*b.secondary_air_inlet.temperature[t] \
            -2.41267e-07 * b.wall_temperature_waterwall[t, 7]*b.wall_temperature_waterwall[t, 8] \
            +6.65347e-06 * b.wall_temperature_waterwall[t, 7]*b.flowrate_coal_raw[t] \
            +3.4501e-06 * b.wall_temperature_waterwall[t, 8]*b.flowrate_coal_raw[t] \
            +2.60247e-06 * b.wall_temperature_waterwall[t, 9]*b.flowrate_coal_raw[t] \
            +4.71651e-06 * b.wall_temperature_waterwall[t, 10]*b.flowrate_coal_raw[t] \
            +1.7254e-05 * b.wall_temperature_platen[t]*b.flowrate_coal_raw[t] \
            +2.48343e-06 * b.wall_temperature_roof[t]*b.flowrate_coal_raw[t] \
            -0.0571714 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            +0.0921271 * b.flowrate_coal_raw[t]*b.SR[t] \
            -0.233398 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            -6.86695e-05 * b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t] \
            +0.000995532 * b.flowrate_coal_raw[t]*b.ratio_PA2coal[t] \
            -0.58226 * b.mf_H2O_coal_raw[t]*b.SR[t] \
            -0.00447073 * b.SR[t]*b.secondary_air_inlet.temperature[t] \
            +0.000246125 * b.secondary_air_inlet.temperature[t]*b.ratio_PA2coal[t] \
            -2.22952e-10 * (b.wall_temperature_platen[t]*b.flowrate_coal_raw[t])**2 \
            +1.80036e-06 * (b.wall_temperature_platen[t]*b.mf_H2O_coal_raw[t])**2 \
            -0.000659596 * (b.flowrate_coal_raw[t]*b.SR[t])**2 \
            +0.00333862 * (b.flowrate_coal_raw[t]*b.SR_lf[t])**2 \
            +1.22954e-09 * (b.flowrate_coal_raw[t]*b.secondary_air_inlet.temperature[t])**2))",
    "NOx": "(-0.00436267 * b.wall_temperature_waterwall[t, 1] \
            -0.0254073 * b.wall_temperature_waterwall[t, 2] \
            +0.0510658 * b.wall_temperature_waterwall[t, 3] \
            -0.0639424 * b.wall_temperature_waterwall[t, 4] \
            -0.172523 * b.wall_temperature_waterwall[t, 5] \
            -0.000482641 * b.wall_temperature_waterwall[t, 6] \
            +0.355125 * b.wall_temperature_waterwall[t, 7] \
            +0.00348034 * b.wall_temperature_waterwall[t, 8] \
            -0.97655 * b.wall_temperature_waterwall[t, 9] \
            +1.31514 * b.wall_temperature_waterwall[t, 10] \
            +0.0262232 * b.wall_temperature_waterwall[t, 11] \
            -0.169549 * b.wall_temperature_waterwall[t, 12] \
            -0.000935994 * b.wall_temperature_platen[t] \
            -0.187802 * b.wall_temperature_roof[t] \
            -285.589 * b.flowrate_coal_raw[t] \
            +3715.18 * b.mf_H2O_coal_raw[t] \
            +4412.39 * b.SR[t] \
            +3191.45 * b.SR_lf[t] \
            -2.80808 * b.secondary_air_inlet.temperature[t] \
            -29.6263 * b.ratio_PA2coal[t] \
            -8.19895 * log(b.wall_temperature_waterwall[t, 7]) \
            -632.752 * log(b.wall_temperature_waterwall[t, 10]) \
            -154.422 * log(b.flowrate_coal_raw[t]) \
            -286.045 * log(b.SR[t]) \
            +2.8111 * b.flowrate_coal_raw[t]**2 \
            -2.51074e-07 * b.wall_temperature_waterwall[t, 10]**3 \
            +1290.32 * b.SR_lf[t]**3 \
            +7.41824e-05 * b.wall_temperature_waterwall[t, 3]*b.wall_temperature_waterwall[t, 4] \
            +0.00014509 * b.wall_temperature_waterwall[t, 3]*b.wall_temperature_waterwall[t, 12] \
            +0.000150378 * b.wall_temperature_waterwall[t, 5]*b.wall_temperature_roof[t] \
            +0.0248537 * b.wall_temperature_waterwall[t, 5]*b.ratio_PA2coal[t] \
            -0.353374 * b.wall_temperature_waterwall[t, 7]*b.SR_lf[t] \
            +1.00281 * b.wall_temperature_waterwall[t, 9]*b.SR_lf[t] \
            -1.81866e-05 * b.wall_temperature_waterwall[t, 10]*b.wall_temperature_waterwall[t, 12] \
            -0.279407 * b.wall_temperature_waterwall[t, 10]*b.mf_H2O_coal_raw[t] \
            +0.000775333 * b.wall_temperature_waterwall[t, 11]*b.flowrate_coal_raw[t] \
            -0.415888 * b.wall_temperature_waterwall[t, 11]*b.mf_H2O_coal_raw[t] \
            +0.0836883 * b.wall_temperature_waterwall[t, 12]*b.SR_lf[t] \
            +0.0307674 * b.wall_temperature_roof[t]*b.ratio_PA2coal[t] \
            -16.1433 * b.flowrate_coal_raw[t]*b.mf_H2O_coal_raw[t] \
            +13.5804 * b.flowrate_coal_raw[t]*b.SR[t] \
            +287.31 * b.flowrate_coal_raw[t]*b.SR_lf[t] \
            -3021.62 * b.mf_H2O_coal_raw[t]*b.SR_lf[t] \
            -4546.15 * b.SR[t]*b.SR_lf[t] \
            +2.61054 * b.SR_lf[t]*b.secondary_air_inlet.temperature[t] \
            -0.000160638 * (b.wall_temperature_waterwall[t, 3]*b.SR_lf[t])**2 \
            -0.161015 * (b.flowrate_coal_raw[t]*b.SR[t])**2 \
            -2.63321 * (b.flowrate_coal_raw[t]*b.SR_lf[t])**2)",
}
