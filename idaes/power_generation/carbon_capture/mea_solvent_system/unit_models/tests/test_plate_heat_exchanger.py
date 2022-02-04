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
from pyomo.environ import (check_optimal_termination,
                           ConcreteModel,
                           Param,
                           RangeSet,
                           units as pyunits,
                           value)
from idaes.core import FlowsheetBlock
from idaes.power_generation.carbon_capture.mea_solvent_system.unit_models.plate_heat_exchanger import (
    PlateHeatExchanger as PHE)

from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock)
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.MEA_solvent \
    import configuration as aqueous_mea
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
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
                             'channels_per_pass': 12,
                             'number_of_divider_plates': 2,
                             "hot_side": {
                                 "property_package": m.fs.hotside_properties
                             },
                             "cold_side": {
                                 "property_package": m.fs.coldside_properties
                             }})

    # Check unit config arguments
    assert len(m.fs.unit.config) == 7


# -----------------------------------------------------------------------------
class TestPHE(object):
    @pytest.fixture(scope="class")
    def phe(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.hotside_properties = GenericParameterBlock(default=aqueous_mea)
        m.fs.coldside_properties = GenericParameterBlock(default=aqueous_mea)

        m.fs.unit = PHE(default={
            'passes': 4,
            'channels_per_pass': 12,
            'number_of_divider_plates': 2,
            "hot_side": {"property_package": m.fs.hotside_properties},
            "cold_side": {"property_package": m.fs.coldside_properties}})

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

        # Fix unit geometry - default values should be correct
        m.fs.unit.plate_length.fix()
        m.fs.unit.plate_width.fix()
        m.fs.unit.plate_thickness.fix()
        m.fs.unit.plate_pact_length.fix()
        m.fs.unit.port_diameter.fix()
        m.fs.unit.plate_therm_cond.fix()
        m.fs.unit.area.fix()

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

        assert hasattr(phe.fs.unit.cold_side, "deltaP")
        assert hasattr(phe.fs.unit.hot_side, "deltaP")

        assert isinstance(phe.fs.unit.number_of_passes, Param)
        assert isinstance(phe.fs.unit.channels_per_pass, Param)

    @pytest.mark.component
    def test_units(self, phe):
        assert_units_equivalent(phe.fs.unit.Re_hot[0], pyunits.dimensionless)
        assert_units_equivalent(phe.fs.unit.Re_cold[0], pyunits.dimensionless)
        assert_units_consistent(phe)

    @pytest.mark.unit
    def test_dof(self, phe):
        assert degrees_of_freedom(phe) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, phe):
        initialization_tester(
            phe,
            duty=(245000, pyunits.W),
            optarg={"bound_push": 1e-8, "mu_init": 1e-8})

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, phe):
        results = solver.solve(phe)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, phe):
        # phe.fs.unit.display()
        assert (pytest.approx(182282.48, rel=1e-5) ==
                value(phe.fs.unit.hot_outlet.pressure[0]))
        assert (pytest.approx(177774.85, rel=1e-5) ==
                value(phe.fs.unit.cold_outlet.pressure[0]))

        assert pytest.approx(701.548, rel=1e-5) == value(
            phe.fs.unit.Re_hot[0])
        assert pytest.approx(217.745, rel=1e-5) == value(
            phe.fs.unit.Re_cold[0])
        assert pytest.approx(4.45171, rel=1e-5) == value(
            phe.fs.unit.Pr_hot[0])
        assert pytest.approx(13.5542, rel=1e-5) == value(
            phe.fs.unit.Pr_cold[0])
        assert pytest.approx(3719.23, rel=1e-5) == value(
            phe.fs.unit.heat_transfer_coefficient_hot_side[0])
        assert pytest.approx(2702.97, rel=1e-5) == value(
            phe.fs.unit.heat_transfer_coefficient_cold_side[0])
        assert pytest.approx(1191.38, rel=1e-5) == value(
            phe.fs.unit.heat_transfer_coefficient[0])

        assert pytest.approx(24.8507, rel=1e-5) == value(
            phe.fs.unit.NTU[0])
        assert pytest.approx(0.971227, rel=1e-5) == value(
            phe.fs.unit.Cratio[0])
        assert pytest.approx(0.680902, rel=1e-5) == value(
            phe.fs.unit.effectiveness[0])

        assert pytest.approx(245771, rel=1e-5) == value(
            phe.fs.unit.heat_duty[0])

        assert (pytest.approx(365.120, rel=1e-5) ==
                value(phe.fs.unit.hot_outlet.temperature[0]))
        assert (pytest.approx(357.123, rel=1e-5) ==
                value(phe.fs.unit.cold_outlet.temperature[0]))

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, phe):
        # Mass conservation test
        assert abs(value(phe.fs.unit.hot_inlet.flow_mol[0] -
                         phe.fs.unit.hot_outlet.flow_mol[0])) <= 1e-6

        for j in phe.fs.hotside_properties.apparent_species_set:
            assert abs(value(
                phe.fs.unit.hot_inlet.flow_mol[0] *
                phe.fs.unit.hot_inlet.mole_frac_comp[0, j] -
                phe.fs.unit.hot_outlet.flow_mol[0] *
                phe.fs.unit.hot_outlet.mole_frac_comp[0, j])) <= 1e-6

        assert abs(value(phe.fs.unit.cold_inlet.flow_mol[0] -
                         phe.fs.unit.cold_outlet.flow_mol[0])) <= 1e-6

        for j in phe.fs.coldside_properties.apparent_species_set:
            assert abs(value(
                phe.fs.unit.cold_inlet.flow_mol[0] *
                phe.fs.unit.cold_inlet.mole_frac_comp[0, j] -
                phe.fs.unit.cold_outlet.flow_mol[0] *
                phe.fs.unit.cold_outlet.mole_frac_comp[0, j])) <= 1e-6

        # Energy conservation test
        assert abs(value(
            phe.fs.unit.hot_side.properties_in[
                0]._enthalpy_flow_term["Liq"] +
            phe.fs.unit.cold_side.properties_in[
                0]._enthalpy_flow_term["Liq"] -
            phe.fs.unit.hot_side.properties_out[
                0]._enthalpy_flow_term["Liq"] -
            phe.fs.unit.cold_side.properties_out[
                0]._enthalpy_flow_term["Liq"])) <= 1e-6

    # @pytest.mark.ui
    # @pytest.mark.unit
    # def test_report(self, phe):
    #     phe.fs.unit.report()
