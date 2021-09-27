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
from unit_models.partialcondenser import PartialCondenser


# Import IDAES Libraries
from idaes.core import FlowsheetBlock
# from idaes.power_generation.carbon_capture.mea_solvent_system.unit_models.partialcondenser\
#     import PartialCondenser
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
class TestPartialCondenser(object):
    @pytest.fixture(scope="class")
    def cond(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up property package
        m.fs.vapor_properties = VaporParameterBlock(
            default={'process_type': ProcessType.stripper})
        m.fs.liquid_properties = LiquidParameterBlock(
            default={'process_type': ProcessType.stripper})

        # create instance of kettle reboiler  on flowsheet
        m.fs.unit = PartialCondenser(default={
                                 "vapor_side": {
                                     "property_package": m.fs.vapor_properties
                                 },
                                 "liquid_side": {
                                     "property_package": m.fs.liquid_properties
                                 }})

        # Fix  input variables
        for t in m.fs.time:
            m.fs.unit.vapor_inlet.flow_mol[t].fix(1.1117)
            m.fs.unit.vapor_inlet.temperature[t].fix(339.33)
            m.fs.unit.vapor_inlet.pressure[t].fix(184360)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "CO2"].fix(0.8817)
            m.fs.unit.vapor_inlet.mole_frac_comp[t, "H2O"].fix(0.1183)


        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, cond):

        assert hasattr(cond.fs.unit, "vapor_inlet")
        assert len(cond.fs.unit.vapor_inlet.vars) == 4
        assert hasattr(cond.fs.unit.vapor_inlet, "flow_mol")
        assert hasattr(cond.fs.unit.vapor_inlet, "mole_frac_comp")
        assert hasattr(cond.fs.unit.vapor_inlet, "temperature")
        assert hasattr(cond.fs.unit.vapor_inlet, "pressure")

        assert hasattr(cond.fs.unit, "vapor_outlet")
        assert len(cond.fs.unit.vapor_outlet.vars) == 4
        assert hasattr(cond.fs.unit.vapor_outlet, "flow_mol")
        assert hasattr(cond.fs.unit.vapor_outlet, "mole_frac_comp")
        assert hasattr(cond.fs.unit.vapor_outlet, "temperature")
        assert hasattr(cond.fs.unit.vapor_outlet, "pressure")

        assert not hasattr(cond.fs.unit, "liquid_inlet")

        assert hasattr(cond.fs.unit, "liquid_outlet")
        assert len(cond.fs.unit.liquid_outlet.vars) == 4
        assert hasattr(cond.fs.unit.liquid_outlet, "flow_mol")
        assert hasattr(cond.fs.unit.liquid_outlet, "mole_frac_comp")
        assert hasattr(cond.fs.unit.liquid_outlet, "temperature")
        assert hasattr(cond.fs.unit.liquid_outlet, "pressure")

        assert hasattr(cond.fs.unit, "condenser_duty")


    @pytest.mark.component
    def test_units(self, cond):
        assert_units_equivalent(cond.fs.unit.vapor_inlet.flow_mol[0],
                                units.mol/units.s)
        assert_units_equivalent(cond.fs.unit.vapor_outlet.flow_mol[0],
                                units.mol/units.s)
        assert_units_equivalent(cond.fs.unit.liquid_outlet.flow_mol[0],
                                units.mol/units.s)
        assert_units_equivalent(cond.fs.unit.condenser_duty[0], units.kW)

    @pytest.mark.unit
    def test_dof(self, cond):
        assert degrees_of_freedom(cond) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, cond):
        initialization_tester(cond)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, cond):
        results = solver.solve(cond)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, cond):
        assert (pytest.approx(184360, abs=1e-2) ==
                value(cond.fs.unit.liquid_outlet.pressure[0]))
        assert (pytest.approx(184360, abs=1e-2) ==
                value(cond.fs.unit.vapor_outlet.pressure[0]))
        assert (pytest.approx(303.15, abs=1e-2) ==
                value(cond.fs.unit.liquid_outlet.temperature[0]))
        assert (pytest.approx(-6.28, abs=1e-2) ==
                value(cond.fs.unit.condenser_duty[0]))
        assert (pytest.approx(303.15, abs=1e-2) ==
                value(cond.fs.unit.vapor_outlet.temperature[0]))
        assert (pytest.approx(1.0034, abs=1e-4) ==
                value(cond.fs.unit.vapor_outlet.flow_mol[0]))
        assert (pytest.approx(0.1083, abs=1e-4) ==
                value(cond.fs.unit.liquid_outlet.flow_mol[0]))
        assert (pytest.approx(0.9026, abs=1e-4) ==
                value(cond.fs.unit.vf[0]))
        assert (pytest.approx(0.0, abs=1e-4) ==
                value(cond.fs.unit.liquid_outlet.mole_frac_comp[0,'CO2']))
        assert (pytest.approx(0.0, abs=1e-4) ==
                value(cond.fs.unit.liquid_outlet.mole_frac_comp[0,'MEA']))
        assert (pytest.approx(1.0, abs=1e-4) ==
                value(cond.fs.unit.liquid_outlet.mole_frac_comp[0,'H2O']))
        assert (pytest.approx(0.9769, abs=1e-4) ==
                value(cond.fs.unit.vapor_outlet.mole_frac_comp[0,'CO2']))
        assert (pytest.approx(0.0231, abs=1e-4) ==
                value(cond.fs.unit.vapor_outlet.mole_frac_comp[0,'H2O']))

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, cond):
        assert abs(value(cond.fs.unit.vapor_inlet.flow_mol[0] -
                         cond.fs.unit.vapor_outlet.flow_mol[0] -
                         cond.fs.unit.liquid_outlet.flow_mol[0])) <= 1e-5

        assert abs(value(cond.fs.unit.condenser_duty[0]  -
                (cond.fs.unit.heat_sens[0]-cond.fs.unit.heat_latent[0])*
                 cond.fs.unit.vapor_inlet.flow_mol[0]*1e-3) <= 1e-5)


