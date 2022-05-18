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
Author: Paul Akula, Anuja Deshpande, Andrew Lee
"""
# Import Python libraries
import pytest

# Import Pyomo libraries
from pyomo.environ import ConcreteModel, value, SolverStatus, TerminationCondition
from pyomo.util.check_units import assert_units_consistent

# Import IDAES Libraries
import idaes
from idaes.core import FlowsheetBlock
from idaes.models_extra.column_models.solvent_column import PackedColumn
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.column_models.properties.MEA_vapor import flue_gas, wet_co2
from idaes.models_extra.column_models.properties.MEA_solvent import (
    configuration as liquid_config,
)

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver


# -----------------------------------------------------------------------------
solver = get_solver()


class TestAbsorberColumn:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up property package
        m.fs.vapor_properties = GenericParameterBlock(default=flue_gas)
        m.fs.liquid_properties = GenericParameterBlock(default=liquid_config)

        # Create an instance of the column in the flowsheet
        m.fs.unit = PackedColumn(
            default={
                "finite_elements": 10,
                "has_pressure_change": False,
                "vapor_phase": {"property_package": m.fs.vapor_properties},
                "liquid_phase": {"property_package": m.fs.liquid_properties},
            }
        )

        # Fix column design variables
        m.fs.unit.diameter_column.fix(0.65)
        m.fs.unit.length_column.fix(15)

        # Fix operating conditions
        # Flue gas
        m.fs.unit.vapor_inlet.flow_mol.fix(22)
        m.fs.unit.vapor_inlet.temperature.fix(318)
        m.fs.unit.vapor_inlet.pressure.fix(107650)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "CO2"].fix(0.12)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "H2O"].fix(0.09)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "N2"].fix(0.74)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "O2"].fix(0.06)
        # Solvent liquid
        m.fs.unit.liquid_inlet.flow_mol.fix(37)
        m.fs.unit.liquid_inlet.temperature.fix(320)
        m.fs.unit.liquid_inlet.pressure.fix(107650)
        m.fs.unit.liquid_inlet.mole_frac_comp[0, "CO2"].fix(0.01)
        m.fs.unit.liquid_inlet.mole_frac_comp[0, "H2O"].fix(0.87)
        m.fs.unit.liquid_inlet.mole_frac_comp[0, "MEA"].fix(0.12)

        m.fs.unit.holdup_liq.fix(1e-2)
        m.fs.unit.mass_transfer_coeff_vap[0, :, "CO2"].fix(3e-7)
        m.fs.unit.mass_transfer_coeff_vap[0, :, "H2O"].fix(4e-7)
        m.fs.unit.heat_transfer_coeff.fix(7100)
        m.fs.unit.area_interfacial.fix(200)

        # Apply scaling
        iscale.calculate_scaling_factors(m.fs.unit)

        return m

    @pytest.mark.unit
    def test_build(self, model):

        assert model.fs.unit.config.dynamic is False

        assert hasattr(model.fs.unit, "vapor_inlet")
        assert hasattr(model.fs.unit, "vapor_outlet")
        assert hasattr(model.fs.unit, "liquid_inlet")
        assert hasattr(model.fs.unit, "liquid_outlet")

    @pytest.mark.unit
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Solver status and condition
        assert results.solver.status == SolverStatus.ok
        assert results.solver.termination_condition == TerminationCondition.optimal

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(22.0991, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.flow_mol[0]
        )
        assert pytest.approx(0.0436641, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.159923, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.73668, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "N2"]
        )
        assert pytest.approx(0.0597309, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "O2"]
        )
        assert pytest.approx(107650, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.pressure[0]
        )
        assert pytest.approx(327.433, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.temperature[0]
        )

        assert pytest.approx(37.1209, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.flow_mol[0]
        )
        assert pytest.approx(0.0550919, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.825299, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.119609, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "MEA"]
        )
        assert pytest.approx(107650, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.pressure[0]
        )
        assert pytest.approx(335.497, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.temperature[0]
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        vap_comp = model.fs.unit.config.vapor_phase.property_package.component_list
        liq_comp = (
            model.fs.unit.config.liquid_phase.property_package.apparent_species_set
        )

        # Mass conservation test
        vap_in = model.fs.unit.vapor_phase.properties[0, 0]
        vap_out = model.fs.unit.vapor_phase.properties[0, 1]
        liq_in = model.fs.unit.liquid_phase.properties[0, 1]
        liq_out = model.fs.unit.liquid_phase.properties[0, 0]

        # Material conservation
        for j in liq_comp:
            if j in vap_comp:
                assert 1e-6 >= abs(
                    value(
                        vap_in.get_material_flow_terms("Vap", j)
                        + liq_in.get_material_flow_terms("Liq", j)
                        - vap_out.get_material_flow_terms("Vap", j)
                        - liq_out.get_material_flow_terms("Liq", j)
                    )
                )
            else:
                assert 1e-6 >= abs(
                    value(
                        liq_in.get_material_flow_terms("Liq", j)
                        - liq_out.get_material_flow_terms("Liq", j)
                    )
                )

        for j in vap_comp:
            if j not in liq_comp:
                assert 1e-6 >= abs(
                    value(
                        vap_in.get_material_flow_terms("Vap", j)
                        - vap_out.get_material_flow_terms("Vap", j)
                    )
                )

        # Energy conservation
        assert 1e-5 >= abs(
            value(
                vap_in.get_enthalpy_flow_terms("Vap")
                + liq_in.get_enthalpy_flow_terms("Liq")
                - vap_out.get_enthalpy_flow_terms("Vap")
                - liq_out.get_enthalpy_flow_terms("Liq")
            )
        )


class TestStripperColumn:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up property package
        m.fs.vapor_properties = GenericParameterBlock(default=wet_co2)
        m.fs.liquid_properties = GenericParameterBlock(default=liquid_config)

        # Create an instance of the column in the flowsheet
        m.fs.unit = PackedColumn(
            default={
                "finite_elements": 10,
                "has_pressure_change": False,
                "vapor_phase": {"property_package": m.fs.vapor_properties},
                "liquid_phase": {"property_package": m.fs.liquid_properties},
            }
        )

        # Fix column design variables
        m.fs.unit.diameter_column.fix(0.64135)
        m.fs.unit.length_column.fix(12.1)

        # Fix operating conditions
        # Flue gas
        m.fs.unit.vapor_inlet.flow_mol.fix(17.496)
        m.fs.unit.vapor_inlet.temperature.fix(396.6)
        m.fs.unit.vapor_inlet.pressure.fix(183430)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "CO2"].fix(0.0145)
        m.fs.unit.vapor_inlet.mole_frac_comp[0, "H2O"].fix(0.9855)
        # Solvent liquid
        m.fs.unit.liquid_inlet.flow_mol.fix(84.48)
        m.fs.unit.liquid_inlet.temperature.fix(382.15)
        m.fs.unit.liquid_inlet.pressure.fix(183430)
        m.fs.unit.liquid_inlet.mole_frac_comp[0, "CO2"].fix(0.0331)
        m.fs.unit.liquid_inlet.mole_frac_comp[0, "H2O"].fix(0.8547)
        m.fs.unit.liquid_inlet.mole_frac_comp[0, "MEA"].fix(0.1122)

        m.fs.unit.holdup_liq.fix(1e-2)
        m.fs.unit.mass_transfer_coeff_vap[0, :, "CO2"].fix(3e-5)
        m.fs.unit.mass_transfer_coeff_vap[0, :, "H2O"].fix(4e-5)
        m.fs.unit.heat_transfer_coeff.fix(7100)
        m.fs.unit.area_interfacial.fix(200)

        # Apply scaling
        iscale.calculate_scaling_factors(m.fs.unit)

        return m

    @pytest.mark.unit
    def test_build(self, model):

        assert model.fs.unit.config.dynamic is False

        assert hasattr(model.fs.unit, "vapor_inlet")
        assert hasattr(model.fs.unit, "vapor_outlet")
        assert hasattr(model.fs.unit, "liquid_inlet")
        assert hasattr(model.fs.unit, "liquid_outlet")

    @pytest.mark.unit
    def test_degrees_of_freedom(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Solver status and condition
        assert results.solver.status == SolverStatus.ok
        assert results.solver.termination_condition == TerminationCondition.optimal

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(11.3723, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.flow_mol[0]
        )
        assert pytest.approx(0.12857, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.87143, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(183430, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.pressure[0]
        )
        assert pytest.approx(403.912, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.temperature[0]
        )

        assert pytest.approx(90.6037, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.flow_mol[0]
        )
        assert pytest.approx(0.0175250, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.877858, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.104616, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "MEA"]
        )
        assert pytest.approx(183430, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.pressure[0]
        )
        assert pytest.approx(394.014, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.temperature[0]
        )

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        vap_comp = model.fs.unit.config.vapor_phase.property_package.component_list
        liq_comp = (
            model.fs.unit.config.liquid_phase.property_package.apparent_species_set
        )

        vap_in = model.fs.unit.vapor_phase.properties[0, 0]
        vap_out = model.fs.unit.vapor_phase.properties[0, 1]
        liq_in = model.fs.unit.liquid_phase.properties[0, 1]
        liq_out = model.fs.unit.liquid_phase.properties[0, 0]

        # Material conservation
        for j in liq_comp:
            if j in vap_comp:
                assert 1e-6 >= abs(
                    value(
                        vap_in.get_material_flow_terms("Vap", j)
                        + liq_in.get_material_flow_terms("Liq", j)
                        - vap_out.get_material_flow_terms("Vap", j)
                        - liq_out.get_material_flow_terms("Liq", j)
                    )
                )
            else:
                assert 1e-6 >= abs(
                    value(
                        liq_in.get_material_flow_terms("Liq", j)
                        - liq_out.get_material_flow_terms("Liq", j)
                    )
                )

        for j in vap_comp:
            if j not in liq_comp:
                assert 1e-6 >= abs(
                    value(
                        vap_in.get_material_flow_terms("Vap", j)
                        - vap_out.get_material_flow_terms("Vap", j)
                    )
                )

        # Energy conservation
        assert 2e-5 >= abs(
            value(
                vap_in.get_enthalpy_flow_terms("Vap")
                + liq_in.get_enthalpy_flow_terms("Liq")
                - vap_out.get_enthalpy_flow_terms("Vap")
                - liq_out.get_enthalpy_flow_terms("Liq")
            )
        )
