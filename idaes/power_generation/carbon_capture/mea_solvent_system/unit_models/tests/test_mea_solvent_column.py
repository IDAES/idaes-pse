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
from pyomo.environ import (ConcreteModel,
                           value,
                           Constraint,
                           SolverStatus,
                           TerminationCondition)
from pyomo.util.check_units import assert_units_consistent

# Import IDAES Libraries
from idaes.core import FlowsheetBlock
from idaes.power_generation.carbon_capture.mea_solvent_system.unit_models.mea_solvent_column \
    import MEAColumn
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.MEA_vapor \
    import flue_gas
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.MEA_solvent \
    import configuration as liquidconfig

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_solver, scaling as iscale

import idaes.logger as idaeslog


# -----------------------------------------------------------------------------
solver = get_solver()


class TestAbsorber:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up property package
        m.fs.vapor_properties = GenericParameterBlock(default=flue_gas)
        m.fs.liquid_properties = GenericParameterBlock(default=liquidconfig)

        # Create an instance of the column in the flowsheet
        m.fs.unit = MEAColumn(default={
            "finite_elements": 30,
            "has_pressure_change": False,
            "vapor_side": {
                "property_package": m.fs.vapor_properties},
            "liquid_side": {
                "property_package": m.fs.liquid_properties
            }})

        m.fs.unit.diameter_column.fix(0.64135)
        m.fs.unit.length_column.fix(18.15)

        # Fix operating conditions
        # Flue gas
        m.fs.unit.vapor_inlet.flow_mol.fix(21.48)
        m.fs.unit.vapor_inlet.temperature.fix(317.88)
        m.fs.unit.vapor_inlet.pressure.fix(107650)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "CO2"].fix(0.11453)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "H2O"].fix(0.08526)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "N2"].fix(0.73821)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "O2"].fix(0.06200)
        # Solvent liquid
        m.fs.unit.liquid_inlet.flow_mol.fix(37.55)
        m.fs.unit.liquid_inlet.temperature.fix(319.87)
        m.fs.unit.liquid_inlet.pressure.fix(107650)
        m.fs.unit.liquid_inlet.mole_frac_comp[0, "CO2"].fix(0.00963)
        m.fs.unit.liquid_inlet.mole_frac_comp[0, "H2O"].fix(0.87435)
        m.fs.unit.liquid_inlet.mole_frac_comp[0, "MEA"].fix(0.11602)

        m.fs.unit.mass_transfer_coeff_vap_comp.fix(3e-5)

        return m

    @pytest.mark.unit
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model) == 0

    # @pytest.mark.component
    # def test_unit_consistency(self, model):
    #     # TODO : Units are still a mess and need to be fixed
    #     assert_units_consistent(model.fs.unit)

    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model, outlvl=idaeslog.DEBUG)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        res = solver.solve(model)

        # Solver status/condition
        assert res.solver.termination_condition == TerminationCondition.optimal
        assert res.solver.status == SolverStatus.ok

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(22.9449, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.flow_mol[0])
        assert pytest.approx(0.0279819, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "CO2"])
        assert pytest.approx(0.222898, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "H2O"])
        assert pytest.approx(0.691079, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "N2"])
        assert pytest.approx(0.0580416, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "O2"])
        assert pytest.approx(107650, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.pressure[0])
        assert pytest.approx(333.578, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.temperature[0])

        assert pytest.approx(36.0851, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.flow_mol[0])
        assert pytest.approx(0.0604036, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "CO2"])
        assert pytest.approx(0.818866, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "H2O"])
        assert pytest.approx(0.120730, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "MEA"])
        assert pytest.approx(107650, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.pressure[0])
        assert pytest.approx(320.760, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.temperature[0])

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        vap_in = model.fs.unit.vapor_phase.properties[0, 0]
        vap_out = model.fs.unit.vapor_phase.properties[0, 1]
        liq_in = model.fs.unit.liquid_phase.properties[0, 1]
        liq_out = model.fs.unit.liquid_phase.properties[0, 0]

        # Material conservation
        for j in ["CO2", "H2O"]:
            assert 1e-6 >= abs(value(
                vap_in.get_material_flow_terms("Vap", j) +
                liq_in.get_material_flow_terms("Liq", j) -
                vap_out.get_material_flow_terms("Vap", j) -
                liq_out.get_material_flow_terms("Liq", j)))
        for j in ["N2", "O2"]:
            assert 1e-6 >= abs(value(
                vap_in.get_material_flow_terms("Vap", j) -
                vap_out.get_material_flow_terms("Vap", j)))
        for j in ["MEA"]:
            assert 1e-6 >= abs(value(
                liq_in.get_material_flow_terms("Liq", j) -
                liq_out.get_material_flow_terms("Liq", j)))

        # Energy conservation
        assert 1e-6 >= abs(value(
            vap_in.get_enthalpy_flow_terms("Vap") +
            liq_in.get_enthalpy_flow_terms("Liq") -
            vap_out.get_enthalpy_flow_terms("Vap") -
            liq_out.get_enthalpy_flow_terms("Liq")))
