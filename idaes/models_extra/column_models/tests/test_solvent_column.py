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
from pyomo.environ import (
    ConcreteModel,
    value,
    SolverStatus,
    TerminationCondition,
    TransformationFactory,
)
from pyomo.util.check_units import assert_units_consistent

# Import IDAES Libraries
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
        m.fs = FlowsheetBlock(dynamic=False)

        # Set up property package
        m.fs.vapor_properties = GenericParameterBlock(**flue_gas)
        m.fs.liquid_properties = GenericParameterBlock(**liquid_config)

        # Create an instance of the column in the flowsheet
        m.fs.unit = PackedColumn(
            finite_elements=10,
            has_pressure_change=False,
            vapor_phase={"property_package": m.fs.vapor_properties},
            liquid_phase={"property_package": m.fs.liquid_properties},
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
        assert pytest.approx(22.649833, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.flow_mol[0]
        )
        assert pytest.approx(0.04261887, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.180333, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.718769, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "N2"]
        )
        assert pytest.approx(0.058278, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "O2"]
        )
        assert pytest.approx(107650, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.pressure[0]
        )
        assert pytest.approx(328.456, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.temperature[0]
        )

        assert pytest.approx(36.5701, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.flow_mol[0]
        )
        assert pytest.approx(0.0559114, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.8226781, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.121410, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "MEA"]
        )
        assert pytest.approx(107650, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.pressure[0]
        )
        assert pytest.approx(334.574, rel=1e-5) == value(
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
        m.fs = FlowsheetBlock(dynamic=False)

        # Set up property package
        m.fs.vapor_properties = GenericParameterBlock(**wet_co2)
        m.fs.liquid_properties_stripper = GenericParameterBlock(**liquid_config)

        # Set the heat of absorption value in stripper
        m.fs.liquid_properties_stripper.CO2.dh_abs_co2.unfix()
        m.fs.liquid_properties_stripper.CO2.dh_abs_co2.fix(-97000)

        # Create an instance of the column in the flowsheet
        m.fs.unit = PackedColumn(
            finite_elements=10,
            has_pressure_change=False,
            vapor_phase={"property_package": m.fs.vapor_properties},
            liquid_phase={"property_package": m.fs.liquid_properties_stripper},
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
        for x in m.fs.unit.liquid_phase.length_domain:
            iscale.set_scaling_factor(
                m.fs.unit.liquid_phase.properties[0, x].mole_frac_phase_comp_true[
                    "Liq", "CO2"
                ],
                1e4,
            )

        for (t, x, j), v in m.fs.unit.pressure_equil.items():
            if x != 0:
                iscale.set_scaling_factor(
                    v,
                    1
                    / value(
                        m.fs.unit.liquid_phase.properties[t, x].fug_phase_comp["Liq", j]
                    ),
                )
            else:
                iscale.set_scaling_factor(v, 1)

        for (t, x, j), v in m.fs.unit.interphase_mass_transfer.items():
            if x != 0:
                iscale.set_scaling_factor(
                    v, 1 / value(m.fs.unit.interphase_mass_transfer[t, x, j])
                )
            else:
                iscale.set_scaling_factor(v, 1)

        for x in m.fs.unit.vapor_phase.length_domain:
            iscale.set_scaling_factor(m.fs.unit.vapor_phase.heat[0, x], 1e-2)

            iscale.set_scaling_factor(
                m.fs.unit.vapor_phase.enthalpy_transfer[0, x], 0.1
            )

            iscale.set_scaling_factor(
                m.fs.unit.vapor_phase._enthalpy_flow[0, x, "Vap"], 1e-4
            )

            iscale.set_scaling_factor(
                m.fs.unit.vapor_phase.enthalpy_flow_dx[0, x, "Vap"], 1e-3
            )

        for x in m.fs.unit.liquid_phase.length_domain:
            iscale.set_scaling_factor(
                m.fs.unit.liquid_phase._enthalpy_flow[0, x, "Liq"], 1e-6
            )

            iscale.set_scaling_factor(
                m.fs.unit.liquid_phase.enthalpy_flow_dx[0, x, "Liq"], 1e-3
            )

            iscale.set_scaling_factor(
                m.fs.unit.liquid_phase.enthalpy_transfer[0, x], 0.1
            )

            iscale.set_scaling_factor(m.fs.unit.liquid_phase.heat[0, x], 1e-2)

        iscale.calculate_scaling_factors(m.fs.unit)

        xfrm = TransformationFactory("contrib.strip_var_bounds")
        xfrm.apply_to(m, reversible=True)

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
        assert pytest.approx(13.50864, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.flow_mol[0]
        )
        assert pytest.approx(0.109306, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.890693, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(183430, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.pressure[0]
        )
        assert pytest.approx(391.686, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.temperature[0]
        )

        assert pytest.approx(88.4673, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.flow_mol[0]
        )
        assert pytest.approx(0.0177850, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.8750719, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.107143, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.mole_frac_comp[0, "MEA"]
        )
        assert pytest.approx(183430, rel=1e-5) == value(
            model.fs.unit.liquid_outlet.pressure[0]
        )
        assert pytest.approx(394.134, rel=1e-5) == value(
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
