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
Test for plate heat exchanger model

Author: Paul Akula
"""
# Import Python libraries
import pytest
import sys
import os

# Import Pyomo libraries
from pyomo.environ import ConcreteModel, value,  SolverFactory

# Import IDAES Libraries
from idaes.core import FlowsheetBlock

# Access mea_column_files dir from the current dir (tests dir)
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from unit_models.phe import PHE
from properties.liquid_prop import LiquidParameterBlock

solver = SolverFactory('ipopt')


m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})
# Set up property package
m.fs.hotside_properties = LiquidParameterBlock()
m.fs.coldside_properties = LiquidParameterBlock()

# create instance of plate heat exchanger  on flowsheet
m.fs.unit = PHE(default={'passes': 4,
                       'channel_list': [12, 12, 12, 12],
                       'divider_plate_number':2,
                       "hot_side": {
                           "property_package": m.fs.hotside_properties
                       },
                       "cold_side":
                       {
                           "property_package": m.fs.coldside_properties
                       }})
for t in m.fs.time:
    # hot fluid
    m.fs.unit.hot_inlet.flow_mol[t].fix(60.54879)
    m.fs.unit.hot_inlet.temperature[t].fix(392.23)
    m.fs.unit.hot_inlet.pressure[t].fix(202650)
    m.fs.unit.hot_inlet.mole_frac_comp[t, "CO2"].fix(0.0158)
    m.fs.unit.hot_inlet.mole_frac_comp[t, "H2O"].fix(0.8747)
    m.fs.unit.hot_inlet.mole_frac_comp[t, "MEA"].fix(0.1095)

    # cold fluid
    m.fs.unit.cold_inlet.flow_mol[t].fix(63.01910)
    m.fs.unit.cold_inlet.temperature[t].fix(326.36)
    m.fs.unit.cold_inlet.pressure[t].fix(202650)
    m.fs.unit.cold_inlet.mole_frac_comp[t, "CO2"].fix(0.0414)
    m.fs.unit.cold_inlet.mole_frac_comp[t, "H2O"].fix(0.1077)
    m.fs.unit.cold_inlet.mole_frac_comp[t, "MEA"].fix(0.8509)

m.fs.unit.initialize()

print('Passes = {}'.format(value(m.fs.unit.P)))

