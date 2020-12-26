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
Tests for steady-state absorber MEA Column
Data Reference No.:  K4
Data Source:
Pilot Solvent Test Unit(PSTU) MEA data at the National Carbon Capture Center(NCCC)
(2014 Test campaign).

Author: Paul Akula, Anuja Deshpande
"""
import sys
import os
from pyomo.environ import ConcreteModel, SolverFactory
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics  import  report_statistics
import pytest


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

#------------------------------------------------------------------------------

def test_build(run=False):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
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
       # Outlet Stream Condition Testing
       assert m.fs.abs.vapor_phase.properties[0, 1].temperature.value == pytest.approx(346.17684550190495,abs=1e-4)
       assert m.fs.abs.vapor_phase.properties[0, 1].pressure.value == pytest.approx(107650,abs=1e-4)
       assert m.fs.abs.vapor_phase.properties[0, 1].flow_mol_comp['CO2'].value == pytest.approx(0.6751702735775524,abs=1e-4)
       assert m.fs.abs.vapor_phase.properties[0, 1].flow_mol_comp['H2O'].value == pytest.approx(6.674303577233841,abs=1e-4)
       assert m.fs.abs.vapor_phase.properties[0, 1].flow_mol_comp['N2'].value == pytest.approx(15.8567508,abs=1e-4)
       assert m.fs.abs.vapor_phase.properties[0, 1].flow_mol_comp['O2'].value == pytest.approx(1.33176,abs=1e-4)

       assert m.fs.abs.liquid_phase.properties[0, 0].temperature.value == pytest.approx(322.95279074281416,abs=1e-4)
       assert m.fs.abs.liquid_phase.properties[0, 0].pressure.value == pytest.approx(107650,abs=1e-4)
       assert m.fs.abs.liquid_phase.properties[0, 0].flow_mol_comp['CO2'].value == pytest.approx(2.146540626422447,abs=1e-4)
       assert m.fs.abs.liquid_phase.properties[0, 0].flow_mol_comp['H2O'].value == pytest.approx(27.988923722766152,abs=1e-4)
       assert m.fs.abs.liquid_phase.properties[0, 0].flow_mol_comp['MEA'].value == pytest.approx(4.356551,abs=1e-4)

       # Performance Testing - CO2 Capture %
       co2_cap_percentage = (1 - (m.fs.abs.vapor_phase.properties[0, 1].flow_mol_comp['CO2'].value/(m.fs.abs.vap_in_mole_frac[0,"CO2"].value*m.fs.abs.vap_in_flow[0].value)))*100
       assert co2_cap_percentage == pytest.approx(72.55521864935682,abs=1e-4)
       report_statistics(m.fs.abs)

if __name__ == "__main__":
          test_build(run=True)


