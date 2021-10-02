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
Author: Paul Akula, Anuja Deshpande
"""
# Import Python libraries
import pytest

# Import Pyomo libraries
from pyomo.environ import ConcreteModel, value, Param, TransformationFactory,\
    TerminationCondition, units as pyunits, Var, Constraint, SolverStatus
from pyomo.util.check_units import assert_units_consistent

# Import IDAES Libraries
from idaes.core import FlowsheetBlock
from idaes.power_generation.carbon_capture.mea_solvent_system.unit_models.column \
    import PackedColumn, ProcessType
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.vapor_prop \
    import VaporParameterBlock
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.liquid_prop \
    import LiquidParameterBlock

from idaes.core import ControlVolume1DBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_solver


# -----------------------------------------------------------------------------
solver = get_solver()


class TestAbsorber:
    @pytest.fixture(scope="class")
    def model(self):
        """Setup for steady-state absorption column"""

        # Spatial domain finite elemets and finite element list
        x_nfe = 10
        x_nfe_list = [i / x_nfe for i in range(x_nfe + 1)]

        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up property package
        m.fs.vapor_properties = VaporParameterBlock(
            default={'process_type': ProcessType.absorber})
        m.fs.liquid_properties = LiquidParameterBlock(
            default={'process_type': ProcessType.absorber})

        # Create instance of column on flowsheet
        m.fs.unit = PackedColumn(default={
            "process_type": ProcessType.absorber,
            "finite_elements": x_nfe,
            "length_domain_set": x_nfe_list,
            "transformation_method": "dae.finite_difference",
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

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert degrees_of_freedom(model) == 14

        assert isinstance(model.fs.unit.liquid_phase, ControlVolume1DBlock)
        assert isinstance(model.fs.unit.vapor_phase, ControlVolume1DBlock)

        var_list = [
            "diameter_column", "area_column", "length_column", "velocity_vap",
            "velocity_liq", "pressure_equil", "interphase_mass_transfer",
            "enhancement_factor", "yi_solvent", "yeq_solute", "heat_vap", "heat_liq"]
        for v in var_list:
            assert isinstance(getattr(model.fs.unit, v), Var)

        const_list = [
            "column_cross_section_area", "vapor_side_length",
            "liquid_side_length", "eq_velocity_vap", "eq_velocity_liq",
            "pressure_at_interface", "mass_transfer",
            "liquid_phase_mass_transfer_handle",
            "vapor_phase_mass_transfer_handle",
            "vapor_phase_heat_transfer", "liquid_phase_heat_transfer"]
        for c in const_list:
            assert isinstance(getattr(model.fs.unit, c), Constraint)

    @pytest.mark.component
    def test_units(self, model):
        assert_units_consistent(model)

    @pytest.mark.component
    def test_initialize(self, model):
        # Set DoF
        model.fs.unit.diameter_column.fix(0.64135)
        model.fs.unit.length_column.fix(18.15)

        model.fs.unit.vapor_inlet.flow_mol[0].fix(21.48)
        model.fs.unit.vapor_inlet.temperature[0].fix(317.88)
        model.fs.unit.vapor_inlet.pressure[0].fix(107650)
        model.fs.unit.vapor_inlet.mole_frac_comp[0, "CO2"].fix(0.11453)
        model.fs.unit.vapor_inlet.mole_frac_comp[0, "H2O"].fix(0.08526)
        model.fs.unit.vapor_inlet.mole_frac_comp[0, "N2"].fix(0.73821)
        model.fs.unit.vapor_inlet.mole_frac_comp[0, "O2"].fix(0.06200)

        model.fs.unit.liquid_inlet.flow_mol[0].fix(37.55)
        model.fs.unit.liquid_inlet.temperature[0].fix(319.87)
        model.fs.unit.liquid_inlet.mole_frac_comp[0, "CO2"].fix(0.00963)
        model.fs.unit.liquid_inlet.mole_frac_comp[0, "H2O"].fix(0.87435)
        model.fs.unit.liquid_inlet.mole_frac_comp[0, "MEA"].fix(0.11602)

        assert degrees_of_freedom(model) == 0
        initialization_tester(model)

        assert pytest.approx(346.177, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.temperature[0])
        assert pytest.approx(322.953, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.temperature[0])

    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

        assert pytest.approx(0.0275159, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "CO2"])
        assert pytest.approx(0.0622325, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "CO2"])

        assert pytest.approx(346.176, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.temperature[0])
        assert pytest.approx(322.953, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.temperature[0])

    @pytest.mark.component
    def test_conservation(self, model):
        vi = model.fs.unit.vapor_inlet
        vo = model.fs.unit.vapor_outlet
        li = model.fs.unit.liquid_inlet
        lo = model.fs.unit.liquid_outlet

        assert 0 == pytest.approx(value(
            vi.flow_mol[0]*vi.mole_frac_comp[0, "N2"] -
            vo.flow_mol[0]*vo.mole_frac_comp[0, "N2"]), abs=1e-3)
        assert 0 == pytest.approx(value(
            vi.flow_mol[0]*vi.mole_frac_comp[0, "O2"] -
            vo.flow_mol[0]*vo.mole_frac_comp[0, "O2"]), abs=1e-3)
        assert 0 == pytest.approx(value(
            vi.flow_mol[0]*vi.mole_frac_comp[0, "CO2"] +
            li.flow_mol[0]*li.mole_frac_comp[0, "CO2"] -
            vo.flow_mol[0]*vo.mole_frac_comp[0, "CO2"] -
            lo.flow_mol[0]*lo.mole_frac_comp[0, "CO2"]), abs=1e-3)
        assert 0 == pytest.approx(value(
            vi.flow_mol[0]*vi.mole_frac_comp[0, "H2O"] +
            li.flow_mol[0]*li.mole_frac_comp[0, "H2O"] -
            vo.flow_mol[0]*vo.mole_frac_comp[0, "H2O"] -
            lo.flow_mol[0]*lo.mole_frac_comp[0, "H2O"]), abs=1e-3)
        assert 0 == pytest.approx(value(
            li.flow_mol[0]*li.mole_frac_comp[0, "MEA"] -
            lo.flow_mol[0]*lo.mole_frac_comp[0, "MEA"]), abs=1e-3)

        # TODO: Add energy conservation equations


class TestStripper:
    @pytest.fixture(scope="class")
    def model(self):
        """Setup for steady-state stripper column"""

        # Spatial domain finite elemets and finite element list
        x_nfe = 10
        x_nfe_list = [i / x_nfe for i in range(x_nfe + 1)]

        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up property package
        m.fs.vapor_properties = VaporParameterBlock(
            default={'process_type': ProcessType.stripper})
        m.fs.liquid_properties = LiquidParameterBlock(
            default={'process_type': ProcessType.stripper})

        # Create instance of column on flowsheet
        m.fs.unit = PackedColumn(default={
            "process_type": ProcessType.stripper,
            "finite_elements": x_nfe,
            "length_domain_set": x_nfe_list,
            "column_pressure": 183430,
            "transformation_method": "dae.finite_difference",
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

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert degrees_of_freedom(model) == 12

        assert isinstance(model.fs.unit.liquid_phase, ControlVolume1DBlock)
        assert isinstance(model.fs.unit.vapor_phase, ControlVolume1DBlock)

        var_list = [
            "diameter_column", "area_column", "length_column", "velocity_vap",
            "velocity_liq", "pressure_equil", "interphase_mass_transfer",
            "enhancement_factor", "yi_solvent", "yeq_solute", "heat_vap", "heat_liq"]
        for v in var_list:
            assert isinstance(getattr(model.fs.unit, v), Var)

        const_list = [
            "column_cross_section_area", "vapor_side_length",
            "liquid_side_length", "eq_velocity_vap", "eq_velocity_liq",
            "pressure_at_interface", "mass_transfer",
            "liquid_phase_mass_transfer_handle",
            "vapor_phase_mass_transfer_handle",
            "vapor_phase_heat_transfer", "liquid_phase_heat_transfer"]
        for c in const_list:
            assert isinstance(getattr(model.fs.unit, c), Constraint)

    @pytest.mark.component
    def test_units(self, model):
        assert_units_consistent(model)

    @pytest.mark.component
    def test_initialize(self, model):
        # Set DoF
        model.fs.unit.diameter_column.fix(0.64135)
        model.fs.unit.length_column.fix(12.1)

        model.fs.unit.vapor_inlet.flow_mol[0].fix(17.496)
        model.fs.unit.vapor_inlet.temperature[0].fix(396.6)
        model.fs.unit.vapor_inlet.pressure[0].fix(183430)
        model.fs.unit.vapor_inlet.mole_frac_comp[0, "CO2"].fix(0.0145)
        model.fs.unit.vapor_inlet.mole_frac_comp[0, "H2O"].fix(0.9855)

        model.fs.unit.liquid_inlet.flow_mol[0].fix(84.48)
        model.fs.unit.liquid_inlet.temperature[0].fix(382.15)
        model.fs.unit.liquid_inlet.mole_frac_comp[0, "CO2"].fix(0.0331)
        model.fs.unit.liquid_inlet.mole_frac_comp[0, "H2O"].fix(0.8547)
        model.fs.unit.liquid_inlet.mole_frac_comp[0, "MEA"].fix(0.1122)

        assert degrees_of_freedom(model) == 0
        initialization_tester(model)

        assert pytest.approx(396.516, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.temperature[0])
        assert pytest.approx(393.866, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.temperature[0])

    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

        assert pytest.approx(0.143431, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "CO2"])
        assert pytest.approx(0.0184088, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "CO2"])

        assert pytest.approx(396.516, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.temperature[0])
        assert pytest.approx(393.866, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.temperature[0])

    @pytest.mark.component
    def test_conservation(self, model):
        model.fs.unit.vapor_outlet.display()
        model.fs.unit.liquid_outlet.display()

        vi = model.fs.unit.vapor_inlet
        vo = model.fs.unit.vapor_outlet
        li = model.fs.unit.liquid_inlet
        lo = model.fs.unit.liquid_outlet

        assert 0 == pytest.approx(value(
            vi.flow_mol[0]*vi.mole_frac_comp[0, "CO2"] +
            li.flow_mol[0]*li.mole_frac_comp[0, "CO2"] -
            vo.flow_mol[0]*vo.mole_frac_comp[0, "CO2"] -
            lo.flow_mol[0]*lo.mole_frac_comp[0, "CO2"]), abs=1e-3)
        assert 0 == pytest.approx(value(
            vi.flow_mol[0]*vi.mole_frac_comp[0, "H2O"] +
            li.flow_mol[0]*li.mole_frac_comp[0, "H2O"] -
            vo.flow_mol[0]*vo.mole_frac_comp[0, "H2O"] -
            lo.flow_mol[0]*lo.mole_frac_comp[0, "H2O"]), abs=1e-3)
        assert 0 == pytest.approx(value(
            li.flow_mol[0]*li.mole_frac_comp[0, "MEA"] -
            lo.flow_mol[0]*lo.mole_frac_comp[0, "MEA"]), abs=1e-3)

        # TODO: Add energy conservation equations


@pytest.mark.integration
class TestDynamic:
    @pytest.fixture(scope="module",
                    params=[ProcessType.absorber, ProcessType.stripper])
    def column_model_dyn(self, request):
        """ Setup for dynamic column"""

        Fpar = request.param
        # Spacial domain finite elemets and finite element list
        x_nfe = 10
        x_nfe_list = [i / x_nfe for i in range(x_nfe + 1)]

        # Time horizon
        t_nfe = 2
        time_set = [0, 4]

        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": True,
                                       'time_units': pyunits.s,
                                       "time_set": time_set})
        # Set up property package

        m.fs.vapor_properties = VaporParameterBlock(
            default={'process_type': Fpar})
        m.fs.liquid_properties = LiquidParameterBlock(
            default={'process_type': Fpar})

        if Fpar == ProcessType.absorber:
            col_pressure = 107650
        elif Fpar == ProcessType.stripper:
            col_pressure = 183430

        # Create instance of column  on flowsheet
        m.fs.unit = PackedColumn(default={
            "process_type": Fpar,
            "finite_elements": x_nfe,
            "length_domain_set": x_nfe_list,
            "column_pressure": col_pressure,
            "transformation_method": "dae.finite_difference",
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

        # Fix inputs variables
        m.fs.unit.diameter_column.fix(0.64135)
        if Fpar == ProcessType.absorber:
            m.fs.unit.length_column.fix(18.15)
            # Rich loading at the bottom of column @ final time
            m.fs.unit.loading = Param(initialize=0.4927155969073804)
            for t in m.fs.time:
                # Vapor
                m.fs.unit.vapor_inlet.flow_mol[t].fix(21.48)
                m.fs.unit.vapor_inlet.temperature[t].fix(317.88)
                m.fs.unit.vapor_inlet.pressure[t].fix(col_pressure)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "CO2"].fix(0.11453)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "H2O"].fix(0.08526)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "N2"].fix(0.73821)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "O2"].fix(0.06200)
                # Liquid
                m.fs.unit.liquid_inlet.flow_mol[t].fix(37.55)
                m.fs.unit.liquid_inlet.temperature[t].fix(319.87)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "CO2"].fix(0.00963)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "H2O"].fix(0.87435)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "MEA"].fix(0.11602)

            # Initialize column
            m.fs.unit.initialize()
        elif Fpar == ProcessType.stripper:
            m.fs.unit.length_column.fix(12.1)
            # Lean loading at the bottom of column @ final time
            m.fs.unit.loading = Param(initialize=0.17982818165156983)
            for t in m.fs.time:
                # Vapor
                m.fs.unit.vapor_inlet.flow_mol[t].fix(17.496)
                m.fs.unit.vapor_inlet.temperature[t].fix(396.6)
                m.fs.unit.vapor_inlet.pressure[t].fix(col_pressure)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "CO2"].fix(0.0145)
                m.fs.unit.vapor_inlet.mole_frac_comp[t, "H2O"].fix(0.9855)
                # Liquid
                m.fs.unit.liquid_inlet.flow_mol[t].fix(84.48)
                m.fs.unit.liquid_inlet.temperature[t].fix(382.15)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "CO2"].fix(0.0331)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "H2O"].fix(0.8547)
                m.fs.unit.liquid_inlet.mole_frac_comp[t, "MEA"].fix(0.1122)

            # Initialize column
            m.fs.unit.initialize(
                homotopy_steps_h=[
                    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])

        res = solver.solve(m)

        # Solver status/condition
        assert res.solver.termination_condition == TerminationCondition.optimal
        assert res.solver.status == SolverStatus.ok

        # Performance Condition Testing at final time
        assert value(m.fs.unit.liquid_phase.properties[
                     m.fs.time.last(), 0].mole_frac_comp['CO2'] /
                     m.fs.unit.liquid_phase.properties[
                     m.fs.time.last(), 0].mole_frac_comp['MEA']) == \
            pytest.approx(value(m.fs.unit.loading), abs=1e-1)
