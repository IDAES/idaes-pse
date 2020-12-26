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
import sys
import os
from pyomo.environ import ConcreteModel, SolverFactory, value
from idaes.core import FlowsheetBlock

# Access mea_column_files dir from the current dir (tests dir)
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from unit_models.phe import PHE
from property_package.liquid_prop import  LiquidParameterBlock


solver = SolverFactory('ipopt')

def test_build(run=False):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    # Set up property package
    m.fs.hotside_properties  = LiquidParameterBlock()
    m.fs.coldside_properties = LiquidParameterBlock()

    #create instance of plate heat exchanger  on flowsheet
    m.fs.hx = PHE(default={'passes':4,
                           'channel_list': [12,12,12,12],
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
        m.fs.hx.hot_inlet.mole_frac[t,"CO2"].fix(0.01585)
        m.fs.hx.hot_inlet.mole_frac[t,"H2O"].fix(0.87457)
        m.fs.hx.hot_inlet.mole_frac[t,"MEA"].fix(0.10958)

        #cold fluid
        m.fs.hx.cold_inlet.flow_mol[t].fix(84.7399)
        m.fs.hx.cold_inlet.temperature[t].fix(326.36)
        m.fs.hx.cold_inlet.pressure[t].fix(182840)
        m.fs.hx.cold_inlet.mole_frac[t,"CO2"].fix(0.04142)
        m.fs.hx.cold_inlet.mole_frac[t,"H2O"].fix(0.85078)
        m.fs.hx.cold_inlet.mole_frac[t,"MEA"].fix(0.10780)

    if run:
       m.fs.hx.initialize(outlvl=0)
       # Testing PHE
       assert value(m.fs.hx.QH[0]) == value(m.fs.hx.QC[0])

if __name__ == "__main__":
          test_build(run=True)


