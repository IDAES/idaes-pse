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
Test for Kettle reboiler model

Author: Paul Akula
"""
# Import Python libraries
import sys
import os
import pytest

# Import Pyomo libraries
from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           units,
                           value,
                           Var)

# TODO
# Remove this after the ModuleNotFoundError is resolved
# Access mea_solvent_system dir from the current dir (tests dir)
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
from unit_models.kettlereboiler import KettleReboiler


# Import IDAES Libraries
from idaes.core import FlowsheetBlock
# from idaes.power_generation.carbon_capture.mea_solvent_system.unit_models.kettlereboiler\
#     import KettleReboiler
from idaes.power_generation.carbon_capture.mea_solvent_system.unit_models.column\
    import ProcessType
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.vapor_prop \
    import VaporParameterBlock
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.liquid_prop \
    import LiquidParameterBlock

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_solver
from pyomo.util.check_units import assert_units_equivalent

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------
class TestKettleReboiler(object):
    @pytest.fixture(scope="class")
    def reb(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up property package
        m.fs.vapor_properties = VaporParameterBlock(
            default={'process_type': ProcessType.stripper})
        m.fs.liquid_properties = LiquidParameterBlock(
            default={'process_type': ProcessType.stripper})

        # create instance of kettle reboiler  on flowsheet
        m.fs.unit = KettleReboiler(default={
                                 "vapor_side": {
                                     "property_package": m.fs.vapor_properties
                                 },
                                 "liquid_side": {
                                     "property_package": m.fs.liquid_properties
                                 }})

        # Fix  input variables
        for t in m.fs.time:
            m.fs.unit.liquid_inlet.flow_mol[t].fix(83.89)
            m.fs.unit.liquid_inlet.temperature[t].fix(392.5)
            m.fs.unit.liquid_inlet.pressure[t].fix(183700)
            m.fs.unit.liquid_inlet.mole_frac_comp[t, "CO2"].fix(0.0326)
            m.fs.unit.liquid_inlet.mole_frac_comp[t, "H2O"].fix(0.8589)
            m.fs.unit.liquid_inlet.mole_frac_comp[t, "MEA"].fix(0.1085)
            m.fs.unit.heat_duty[t].fix(430.61)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, reb):

        assert hasattr(reb.fs.unit, "liquid_inlet")
        assert len(reb.fs.unit.liquid_inlet.vars) == 4
        assert hasattr(reb.fs.unit.liquid_inlet, "flow_mol")
        assert hasattr(reb.fs.unit.liquid_inlet, "mole_frac_comp")
        assert hasattr(reb.fs.unit.liquid_inlet, "temperature")
        assert hasattr(reb.fs.unit.liquid_inlet, "pressure")

        assert hasattr(reb.fs.unit, "liquid_outlet")
        assert len(reb.fs.unit.liquid_outlet.vars) == 4
        assert hasattr(reb.fs.unit.liquid_outlet, "flow_mol")
        assert hasattr(reb.fs.unit.liquid_outlet, "mole_frac_comp")
        assert hasattr(reb.fs.unit.liquid_outlet, "temperature")
        assert hasattr(reb.fs.unit.liquid_outlet, "pressure")

        assert not hasattr(reb.fs.unit, "vapor_inlet")

        assert hasattr(reb.fs.unit, "vapor_outlet")
        assert len(reb.fs.unit.vapor_outlet.vars) == 4
        assert hasattr(reb.fs.unit.vapor_outlet, "flow_mol")
        assert hasattr(reb.fs.unit.vapor_outlet, "mole_frac_comp")
        assert hasattr(reb.fs.unit.vapor_outlet, "temperature")
        assert hasattr(reb.fs.unit.vapor_outlet, "pressure")

        assert hasattr(reb.fs.unit, "heat_duty")


    @pytest.mark.component
    def test_units(self, reb):
        assert_units_equivalent(reb.fs.unit.liquid_inlet.flow_mol[0], units.mol/units.s)
        assert_units_equivalent(reb.fs.unit.liquid_outlet.flow_mol[0], units.mol/units.s)
        assert_units_equivalent(reb.fs.unit.vapor_outlet.flow_mol[0], units.mol/units.s)
        assert_units_equivalent(reb.fs.unit.heat_duty[0], units.kW)

    @pytest.mark.unit
    def test_dof(self, reb):
        assert degrees_of_freedom(reb) == 0


    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, reb):
        initialization_tester(reb)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, reb):
        results = solver.solve(reb)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, reb):
        assert (pytest.approx(183700, abs=1e-2) ==
                value(reb.fs.unit.liquid_outlet.pressure[0]))
        assert (pytest.approx(183700, abs=1e-2) ==
                value(reb.fs.unit.vapor_outlet.pressure[0]))
        assert (pytest.approx(393.8499, abs=1e-1) ==
                value(reb.fs.unit.liquid_outlet.temperature[0]))
        assert (pytest.approx(393.8499, abs=1e-1) ==
                value(reb.fs.unit.vapor_outlet.temperature[0]))
        assert (pytest.approx(9.5571, abs=1e-2) ==
                value(reb.fs.unit.vapor_outlet.flow_mol[0]))
        assert (pytest.approx(74.33, abs=1e-2) ==
                value(reb.fs.unit.liquid_outlet.flow_mol[0]))
        assert (pytest.approx(0.1139, abs=1e-2) ==
                value(reb.fs.unit.vf[0]))
        assert (pytest.approx(0.0285, abs=1e-3) ==
                value(reb.fs.unit.liquid_outlet.mole_frac_comp[0,'CO2']))
        assert (pytest.approx(0.1224, abs=1e-3) ==
                value(reb.fs.unit.liquid_outlet.mole_frac_comp[0,'MEA']))
        assert (pytest.approx(0.8491, abs=1e-3) ==
                value(reb.fs.unit.liquid_outlet.mole_frac_comp[0,'H2O']))
        assert (pytest.approx(0.0645, abs=1e-3) ==
                value(reb.fs.unit.vapor_outlet.mole_frac_comp[0,'CO2']))
        assert (pytest.approx(0.9355, abs=1e-3) ==
                value(reb.fs.unit.vapor_outlet.mole_frac_comp[0,'H2O']))

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, reb):
        assert abs(value(reb.fs.unit.liquid_inlet.flow_mol[0] -
                         reb.fs.unit.liquid_outlet.flow_mol[0] -
                         reb.fs.unit.vapor_outlet.flow_mol[0])) <= 1e-5

        assert abs(value(reb.fs.unit.heat_duty[0]  -
                (reb.fs.unit.heat_latent[0] + reb.fs.unit.heat_sens[0])*
                 reb.fs.unit.liquid_inlet.flow_mol[0]) <= 1e-5)


