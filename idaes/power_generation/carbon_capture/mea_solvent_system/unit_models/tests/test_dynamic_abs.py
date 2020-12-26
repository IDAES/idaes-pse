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
Tests for dynamic MEA Column

Author: Paul Akula
"""
import sys
import os
from pyomo.environ import ConcreteModel, SolverFactory,TransformationFactory,\
                          units as pyunits
from idaes.core import FlowsheetBlock

from idaes.core.util.model_statistics  import  report_statistics

# Access mea_column_files dir from the current dir (tests dir)
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from unit_models.column import PackedColumn
from property_package.vapor_prop  import  VaporParameterBlock
from property_package.liquid_prop import  LiquidParameterBlock

# -----------------------------------------------------------------------------
solver = SolverFactory('ipopt')

# spacial domain finite elemets and finite element list
x_nfe = 10
x_nfe_list = [i/x_nfe for i in range(x_nfe+1)]

#time horizon
t_nfe = 2
time_set = [0, 4]
#------------------------------------------------------------------------------

def test_build(run=False):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": True,
                                   'time_units':pyunits.s,
                                   "time_set": time_set})
    # Set up property package
    m.fs.vapor_properties  = VaporParameterBlock()
    m.fs.liquid_properties = LiquidParameterBlock()

    #create instance of column for absorption process on flowsheet
    m.fs.abs = PackedColumn(default={
                      "process_type": "Absorber",
                      "finite_elements": x_nfe,
                      "length_domain_set":x_nfe_list,
                      "transformation_method": "dae.finite_difference",
                      "flow_type":"counter_current",
                      "vapor_side": {
                                   "transformation_scheme": "BACKWARD",
                                   "property_package": m.fs.vapor_properties,
                                   "has_pressure_change": False,
                                   "pressure_drop_type": None},
                      "liquid_side":
                                   {
                                   "transformation_scheme": "FORWARD",
                                   "property_package": m.fs.liquid_properties
                                     }})
    # Time discretization
    discretizer = TransformationFactory('dae.finite_difference')
    discretizer.apply_to(m.fs, wrt=m.fs.time, nfe=t_nfe,scheme='BACKWARD')

    # Fix  input variables
    m.fs.abs.dia_col.fix(0.64135)
    m.fs.abs.length_col.fix(18.15)
    for t in m.fs.time:
        #vapor
        m.fs.abs.vap_in_flow[t].fix(21.48)
        m.fs.abs.vap_in_temperature[t].fix(317.88)
        m.fs.abs.bot_pressure[t].fix(107650)
        m.fs.abs.vap_in_mole_frac[t,"CO2"].fix(0.11453)
        m.fs.abs.vap_in_mole_frac[t,"H2O"].fix(0.08526)
        m.fs.abs.vap_in_mole_frac[t,"N2"].fix(0.73821)
        m.fs.abs.vap_in_mole_frac[t,"O2"].fix(0.06200)
        #liquid
        m.fs.abs.liq_in_flow[t].fix(37.55)
        m.fs.abs.liq_in_temperature[t].fix(319.87)
        m.fs.abs.liq_in_mole_frac[t,"CO2"].fix(0.00963)
        m.fs.abs.liq_in_mole_frac[t,"H2O"].fix(0.87435)
        m.fs.abs.liq_in_mole_frac[t,"MEA"].fix(0.11602)

    assert m.fs.abs.config.liquid_side.transformation_scheme == 'FORWARD'
    assert m.fs.abs.config.vapor_side.transformation_scheme ==  'BACKWARD'
    assert m.fs.abs.config.process_type == 'Absorber'

    if run:
      m.fs.abs.initialize(outlvl=0)
      report_statistics(m.fs.abs)

if __name__ == "__main__":
          test_build(run=True)


