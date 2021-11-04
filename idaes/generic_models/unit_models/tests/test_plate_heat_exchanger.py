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
Tests for Plate Heat Exchnager  unit model.
Author: Akula Paul
"""

import pytest
from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           units,
                           value,
                           Var)
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)
from idaes.generic_models.unit_models.plate_heat_exchanger import (
    PlateHeatExchanger as PHE)

from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock)
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.MEA_solvent \
    import configuration as aqueous_mea
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import (PhysicalParameterTestBlock,
                                     ReactionParameterTestBlock,
                                     initialization_tester)
from idaes.core.util import get_solver
from pyomo.util.check_units import (assert_units_consistent,
                                    assert_units_equivalent)


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.hotside_properties = GenericParameterBlock(default=aqueous_mea)
    m.fs.coldside_properties = GenericParameterBlock(default=aqueous_mea)

    m.fs.unit = PHE(default={'passes': 4,
                             'channel_list': [12, 12, 12, 12],
                             'divider_plate_number': 2,
                             "hot_side": {
                                 "property_package": m.fs.hotside_properties
                             },
                             "cold_side": {
                                 "property_package": m.fs.coldside_properties
                             }})

    # Check unit config arguments
    assert len(m.fs.unit.config) == 16
    assert len(m.fs.unit.config.hot_side) == 2
    assert len(m.fs.unit.config.cold_side) == 2


# -----------------------------------------------------------------------------
class TestPHE(object):
    @pytest.fixture(scope="class")
    def phe(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.hotside_properties = GenericParameterBlock(default=aqueous_mea)
        m.fs.coldside_properties = GenericParameterBlock(default=aqueous_mea)

        m.fs.unit = PHE(default={'passes': 4,
                                 'channel_list': [12, 12, 12, 12],
                                 'divider_plate_number': 2,
                                 "hot_side": {
                                     "property_package": m.fs.hotside_properties
                                 },
                                 "cold_side": {
                                     "property_package": m.fs.coldside_properties
                                 }})
        # hot fluid
        m.fs.unit.hot_inlet.flow_mol[0].fix(60.54879)
        m.fs.unit.hot_inlet.temperature[0].fix(392.23)
        m.fs.unit.hot_inlet.pressure[0].fix(202650)
        m.fs.unit.hot_inlet.mole_frac_comp[0, "CO2"].fix(0.0158)
        m.fs.unit.hot_inlet.mole_frac_comp[0, "H2O"].fix(0.8747)
        m.fs.unit.hot_inlet.mole_frac_comp[0, "MEA"].fix(0.1095)

        # cold fluid
        m.fs.unit.cold_inlet.flow_mol[0].fix(63.01910)
        m.fs.unit.cold_inlet.temperature[0].fix(326.36)
        m.fs.unit.cold_inlet.pressure[0].fix(202650)
        m.fs.unit.cold_inlet.mole_frac_comp[0, "CO2"].fix(0.0414)
        m.fs.unit.cold_inlet.mole_frac_comp[0, "H2O"].fix(0.8509)
        m.fs.unit.cold_inlet.mole_frac_comp[0, "MEA"].fix(0.1077)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, phe):

        assert hasattr(phe.fs.unit, "hot_inlet")
        assert len(phe.fs.unit.hot_inlet.vars) == 4
        assert hasattr(phe.fs.unit.hot_inlet, "flow_mol")
        assert hasattr(phe.fs.unit.hot_inlet, "mole_frac_comp")
        assert hasattr(phe.fs.unit.hot_inlet, "temperature")
        assert hasattr(phe.fs.unit.hot_inlet, "pressure")

        assert hasattr(phe.fs.unit, "hot_outlet")
        assert len(phe.fs.unit.hot_outlet.vars) == 4
        assert hasattr(phe.fs.unit.hot_outlet, "flow_mol")
        assert hasattr(phe.fs.unit.hot_outlet, "mole_frac_comp")
        assert hasattr(phe.fs.unit.hot_outlet, "temperature")
        assert hasattr(phe.fs.unit.hot_outlet, "pressure")

        assert hasattr(phe.fs.unit, "cold_inlet")
        assert len(phe.fs.unit.cold_inlet.vars) == 4
        assert hasattr(phe.fs.unit.cold_inlet, "flow_mol")
        assert hasattr(phe.fs.unit.cold_inlet, "mole_frac_comp")
        assert hasattr(phe.fs.unit.cold_inlet, "temperature")
        assert hasattr(phe.fs.unit.cold_inlet, "pressure")

        assert hasattr(phe.fs.unit, "cold_outlet")
        assert len(phe.fs.unit.cold_outlet.vars) == 4
        assert hasattr(phe.fs.unit.cold_outlet, "flow_mol")
        assert hasattr(phe.fs.unit.cold_outlet, "mole_frac_comp")
        assert hasattr(phe.fs.unit.cold_outlet, "temperature")
        assert hasattr(phe.fs.unit.cold_outlet, "pressure")

        assert hasattr(phe.fs.unit.cold_fluid, "deltaP")
        assert hasattr(phe.fs.unit.hot_fluid, "deltaP")

    @pytest.mark.component
    def test_units(self, phe):
        #assert_units_consistent(phe)
        assert_units_equivalent(phe.fs.unit.plate_length, units.m)
        assert_units_equivalent(phe.fs.unit.cold_fluid.deltaP[0], units.Pa)
        assert_units_equivalent(phe.fs.unit.hot_fluid.deltaP[0], units.Pa)


    @pytest.mark.unit
    def test_dof(self, phe):
        assert degrees_of_freedom(phe) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, phe):
        initialization_tester(phe)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, phe):
        results = solver.solve(phe)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, phe):
        assert (pytest.approx(182282.48, abs=1e-2) ==
                value(phe.fs.unit.hot_outlet.pressure[0]))
        assert (pytest.approx(177774.85, abs=1e-2) ==
                value(phe.fs.unit.cold_outlet.pressure[0]))

        assert (pytest.approx(329.54, abs=1e-2) ==
                value(phe.fs.unit.hot_outlet.temperature[0]))
        assert (pytest.approx(385.32, abs=1e-2) ==
                value(phe.fs.unit.cold_outlet.temperature[0]))

        assert (pytest.approx(0.015800, abs=1e-4) ==
                value(phe.fs.unit.hot_outlet.mole_frac_comp[0, "CO2"]))
        assert (pytest.approx(0.10950, abs=1e-4) ==
                value(phe.fs.unit.hot_outlet.mole_frac_comp[0, "MEA"]))
        assert (pytest.approx(0.87470, abs=1e-4) ==
                value(phe.fs.unit.hot_outlet.mole_frac_comp[0, "H2O"]))

        assert (pytest.approx(0.041400, abs=1e-4) ==
                value(phe.fs.unit.cold_outlet.mole_frac_comp[0, "CO2"]))
        assert (pytest.approx(0.10770, abs=1e-4) ==
                value(phe.fs.unit.cold_outlet.mole_frac_comp[0, "MEA"]))
        assert (pytest.approx(0.85090, abs=1e-4) ==
                value(phe.fs.unit.cold_outlet.mole_frac_comp[0, "H2O"]))

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, phe):
        # Mass conservation test
        assert abs(value(phe.fs.unit.hot_inlet.flow_mol[0] -
                         phe.fs.unit.hot_outlet.flow_mol[0])) <= 1e-6

        assert abs(value(phe.fs.unit.cold_inlet.flow_mol[0] -
                         phe.fs.unit.cold_outlet.flow_mol[0])) <= 1e-6

        # Energy conservation test
        assert abs(value(phe.fs.unit.heat_lost[0] -
                         phe.fs.unit.heat_gain[0])) <=1e-6

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, phe):
        phe.fs.unit.report()
