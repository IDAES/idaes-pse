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
"""
Dictionary of surrogate models

@author: Boiler Subsystem Team (J. Ma, M. Zamarripa)
"""
# A dictionary of surrogate models needs to be provided for the boiler fire
# side model to work. Note that the user needs to provide a surrogate model or
# a fixed value for the heat flux for each water wall section, platen super
# heater (if installed), roof (if installed), flyash, and NOx prediction.
# in this test case the model builds 12 water wall zones, a platen sh, and roof
# usually a surrogate model is a function of the input variables:
# b.wall_temperature_waterwall[t, 1]
# b.wall_temperature_waterwall[t, 2]
# b.wall_temperature_waterwall[t, 3]
# b.wall_temperature_waterwall[t, 4]
# b.wall_temperature_waterwall[t, 5]
# b.wall_temperature_waterwall[t, 6]
# b.wall_temperature_waterwall[t, 7]
# b.wall_temperature_waterwall[t, 8]
# b.wall_temperature_waterwall[t, 9]
# b.wall_temperature_waterwall[t, 10] - or as many water wall zones
# b.wall_temperature_platen[t] - if has_platen_superheater=True
# b.wall_temperature_roof[t] - if has_roof_superheater=True
# b.flowrate_coal_raw[t]
# b.mf_H2O_coal_raw[t] - moisture content of raw coal
# b.SR[t] - stoichiometric ratio
# b.SR_lf[t] - lower furnace stoichiometric ratio
# b.secondary_air_inlet.temperature[t]
# b.ratio_PA2coal[t] - ratio primary air to coal flowrate calculated by model

# Since the development of surrogate models requires a specific case and
# the fire side must be connected with other unit models. This example aims
# testing the model. Therefore, the following dictionary is a simple example
# using fixed values for heat flux to the water wall, platen SH, and roof SH.
data_dic = {
    1: "2.0e7",  # heat flux to water wall zone 1 from fire side in W
    2: "1.0e7",  # heat flux to water wall zone 2 from fire side in W
    3: "1.0e7",  # heat flux to water wall zone 3 from fire side in W
    4: "1.0e7",  # heat flux to water wall zone 4 from fire side in W
    5: "1.5e7",  # heat flux to water wall zone 5 from fire side in W
    6: "1.0e7",  # heat flux to water wall zone 6 from fire side in W
    7: "1.2e7",  # heat flux to water wall zone 7 from fire side in W
    8: "3e7",  # heat flux to water wall zone 8 from fire side in W
    9: "2.5e7",  # heat flux to water wall zone 9 from fire side in W
    10: "2.0e7",  # heat flux to water wall zone 10 from fire side in W
    11: "1.8e7",  # heat flux to water wall zone 11 from fire side in W
    12: "1.0e7",  # heat flux to water wall zone 12 from fire side in W
    "pl": "5.0e7",  # heat flux to platen sh from fire side in W
    "roof": "6.5e7",  # heat flux to roof from fire side in W
    "flyash": "0.0001",  # flyash or unburned carbon mass fraction
    "NOx": "140",
}  # NOx PPM
