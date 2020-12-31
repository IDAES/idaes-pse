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
Tests for dynamic stripper  MEA Column

Author: Paul Akula, Anuja Deshpande
"""
# Import Python libraries
import sys
import os

# Import Pyomo libraries
from pyomo.environ import ConcreteModel, SolverFactory, TransformationFactory,\
    units as pyunits

# Import IDAES Libraries
from idaes.core import FlowsheetBlock

# Access mea_column_files dir from the current dir (tests dir)
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from unit_models.column import PackedColumn
from property_package.vapor_prop import VaporParameterBlock
from property_package.liquid_prop import LiquidParameterBlock

# -----------------------------------------------------------------------------
solver = SolverFactory('ipopt')

# spacial domain finite elemets and finite element list
x_nfe = 10
x_nfe_list = [i / x_nfe for i in range(x_nfe + 1)]

# time horizon
t_nfe = 2
time_set = [0, 4]
# ------------------------------------------------------------------------------


def test_build(run=False):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": True,
                                   'time_units': pyunits.s,
                                   "time_set": time_set})
    # Set up property package
    m.fs.vapor_properties = VaporParameterBlock(
        default={'process_type': 'Stripper'})
    m.fs.liquid_properties = LiquidParameterBlock()

    # create instance of column for stripping process on flowsheet
    m.fs.reg = PackedColumn(default={
        "process_type": "Stripper",
        "finite_elements": x_nfe,
        "length_domain_set": x_nfe_list,
        "transformation_method": "dae.finite_difference",
        "flow_type": "counter_current",
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
    discretizer.apply_to(m.fs, wrt=m.fs.time, nfe=t_nfe, scheme='BACKWARD')

    # Fix  input variables
    m.fs.reg.dia_col.fix(0.64135)
    m.fs.reg.length_col.fix(12.1)
    for t in m.fs.time:
        # vapor
        m.fs.reg.vap_in_flow[t].fix(17.496)
        m.fs.reg.vap_in_temperature[t].fix(396.6)
        m.fs.reg.bot_pressure[t].fix(183430)
        m.fs.reg.vap_in_mole_frac[t, "CO2"].fix(0.0145)
        m.fs.reg.vap_in_mole_frac[t, "H2O"].fix(0.9855)

        # liquid
        m.fs.reg.liq_in_flow[t].fix(84.48)
        m.fs.reg.liq_in_temperature[t].fix(382.15)
        m.fs.reg.liq_in_mole_frac[t, "CO2"].fix(0.0331)
        m.fs.reg.liq_in_mole_frac[t, "H2O"].fix(0.8547)
        m.fs.reg.liq_in_mole_frac[t, "MEA"].fix(0.1122)

    assert m.fs.reg.config.liquid_side.transformation_scheme == 'FORWARD'
    assert m.fs.reg.config.vapor_side.transformation_scheme == 'BACKWARD'
    assert m.fs.reg.config.process_type == 'Stripper'

    if run:
        m.fs.reg.initialize(outlvl=0,
                            homotopy_steps_h=[0.1, 0.2, 0.4, 0.6, 0.8, 1])


if __name__ == "__main__":
    test_build(run=True)
