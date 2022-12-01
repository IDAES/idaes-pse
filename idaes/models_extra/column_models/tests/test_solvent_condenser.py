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
Tests for solvent condenser unit model.
Authors: Andrew Lee
"""

import pytest
from pyomo.environ import (
    check_optimal_termination,
    TransformationFactory,
    ConcreteModel,
    Constraint,
    Param,
    units,
    value,
)
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver

from idaes.models_extra.column_models.solvent_condenser import SolventCondenser
from idaes.models_extra.column_models.properties.MEA_solvent import (
    configuration as aqueous_mea,
)
from idaes.models_extra.column_models.properties.MEA_vapor import wet_co2
from idaes.core.util.exceptions import InitializationError


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
class TestStripperVaporFlow(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.liquid_properties = GenericParameterBlock(**aqueous_mea)
        m.fs.vapor_properties = GenericParameterBlock(**wet_co2)

        m.fs.unit = SolventCondenser(
            liquid_property_package=m.fs.liquid_properties,
            vapor_property_package=m.fs.vapor_properties,
        )

        m.fs.unit.inlet.flow_mol[0].fix(1.1117)
        m.fs.unit.inlet.temperature[0].fix(339.33)
        m.fs.unit.inlet.pressure[0].fix(184360)
        m.fs.unit.inlet.mole_frac_comp[0, "CO2"].fix(0.8817)
        m.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix(0.1183)

        m.fs.unit.reflux.flow_mol[0].fix(0.1083)

        iscale.set_scaling_factor(
            m.fs.unit.vapor_phase.mass_transfer_term[0, "Vap", "CO2"], 1e5
        )
        iscale.set_scaling_factor(
            m.fs.unit.vapor_phase.mass_transfer_term[0, "Vap", "H2O"], 100
        )

        iscale.set_scaling_factor(
            m.fs.unit.vapor_phase.properties_out[0].pressure, 1e-5
        )
        iscale.set_scaling_factor(
            m.fs.unit.vapor_phase.properties_out[0].fug_phase_comp["Vap", "CO2"], 1e-5
        )
        iscale.set_scaling_factor(
            m.fs.unit.vapor_phase.properties_out[0].fug_phase_comp["Vap", "H2O"], 1e-4
        )
        iscale.set_scaling_factor(
            m.fs.unit.vapor_phase.properties_out[0].temperature, 1e-2
        )
        iscale.set_scaling_factor(
            m.fs.unit.vapor_phase.properties_out[0].enth_mol_phase["Vap"], 1e-3
        )

        iscale.set_scaling_factor(m.fs.unit.vapor_phase.enthalpy_transfer[0], 1e-3)

        iscale.calculate_scaling_factors(m.fs.unit)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, model):

        assert hasattr(model.fs.unit, "inlet")
        assert len(model.fs.unit.inlet.vars) == 4
        assert hasattr(model.fs.unit.inlet, "flow_mol")
        assert hasattr(model.fs.unit.inlet, "mole_frac_comp")
        assert hasattr(model.fs.unit.inlet, "temperature")
        assert hasattr(model.fs.unit.inlet, "pressure")

        assert hasattr(model.fs.unit, "reflux")
        assert len(model.fs.unit.reflux.vars) == 4
        assert hasattr(model.fs.unit.reflux, "flow_mol")
        assert hasattr(model.fs.unit.reflux, "mole_frac_comp")
        assert hasattr(model.fs.unit.reflux, "temperature")
        assert hasattr(model.fs.unit.reflux, "pressure")

        assert hasattr(model.fs.unit, "vapor_outlet")
        assert len(model.fs.unit.vapor_outlet.vars) == 4
        assert hasattr(model.fs.unit.vapor_outlet, "flow_mol")
        assert hasattr(model.fs.unit.vapor_outlet, "mole_frac_comp")
        assert hasattr(model.fs.unit.vapor_outlet, "temperature")
        assert hasattr(model.fs.unit.vapor_outlet, "pressure")

        assert isinstance(model.fs.unit.unit_material_balance, Constraint)
        assert isinstance(model.fs.unit.unit_enthalpy_balance, Constraint)
        assert isinstance(model.fs.unit.unit_temperature_equality, Constraint)
        assert isinstance(model.fs.unit.unit_pressure_balance, Constraint)
        assert isinstance(model.fs.unit.zero_flow_param, Param)

        assert number_variables(model.fs.unit) == 55
        assert number_total_constraints(model.fs.unit) == 49
        assert number_unused_variables(model.fs.unit) == 0

    @pytest.mark.component
    def test_units(self, model):
        assert_units_consistent(model)
        assert_units_equivalent(model.fs.unit.heat_duty[0], units.W)

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

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

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(0.1083, rel=1e-5) == value(
            model.fs.unit.reflux.flow_mol[0]
        )
        assert pytest.approx(0, abs=1e-3) == value(
            model.fs.unit.reflux.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0, abs=1e-3) == value(
            model.fs.unit.reflux.mole_frac_comp[0, "MEA"]
        )
        assert pytest.approx(1, rel=1e-3) == value(
            model.fs.unit.reflux.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(184360, rel=1e-5) == value(
            model.fs.unit.reflux.pressure[0]
        )
        assert pytest.approx(303.244, rel=1e-5) == value(
            model.fs.unit.reflux.temperature[0]
        )

        assert pytest.approx(1.0034, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.flow_mol[0]
        )
        assert pytest.approx(0.976758, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.0232423, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(184360, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.pressure[0]
        )
        assert pytest.approx(303.244, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.temperature[0]
        )

        assert pytest.approx(-6278.10, rel=1e-5) == value(model.fs.unit.heat_duty[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        assert (
            abs(
                value(
                    model.fs.unit.inlet.flow_mol[0]
                    - model.fs.unit.reflux.flow_mol[0]
                    - model.fs.unit.vapor_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    model.fs.unit.inlet.flow_mol[0]
                    * model.fs.unit.inlet.mole_frac_comp[0, "CO2"]
                    - model.fs.unit.reflux.flow_mol[0]
                    * model.fs.unit.reflux.mole_frac_comp[0, "CO2"]
                    - model.fs.unit.vapor_outlet.flow_mol[0]
                    * model.fs.unit.vapor_outlet.mole_frac_comp[0, "CO2"]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    model.fs.unit.inlet.flow_mol[0]
                    * model.fs.unit.inlet.mole_frac_comp[0, "H2O"]
                    - model.fs.unit.reflux.flow_mol[0]
                    * model.fs.unit.reflux.mole_frac_comp[0, "H2O"]
                    - model.fs.unit.vapor_outlet.flow_mol[0]
                    * model.fs.unit.vapor_outlet.mole_frac_comp[0, "H2O"]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    model.fs.unit.reflux.flow_mol[0]
                    * model.fs.unit.reflux.mole_frac_comp[0, "MEA"]
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    model.fs.unit.vapor_phase.properties_in[0]._enthalpy_flow_term[
                        "Vap"
                    ]
                    - model.fs.unit.vapor_phase.properties_out[0]._enthalpy_flow_term[
                        "Vap"
                    ]
                    - model.fs.unit.liquid_phase[0]._enthalpy_flow_term["Liq"]
                    + model.fs.unit.heat_duty[0]
                )
            )
            <= 1e-6
        )


# -----------------------------------------------------------------------------
class TestStripperHeatDuty(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.liquid_properties = GenericParameterBlock(**aqueous_mea)
        m.fs.vapor_properties = GenericParameterBlock(**wet_co2)

        m.fs.unit = SolventCondenser(
            liquid_property_package=m.fs.liquid_properties,
            vapor_property_package=m.fs.vapor_properties,
        )

        m.fs.unit.inlet.flow_mol[0].fix(1.1117)
        m.fs.unit.inlet.temperature[0].fix(339.33)
        m.fs.unit.inlet.pressure[0].fix(184360)
        m.fs.unit.inlet.mole_frac_comp[0, "CO2"].fix(0.8817)
        m.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix(0.1183)

        m.fs.unit.heat_duty.fix(-6264)

        xfrm = TransformationFactory("contrib.strip_var_bounds")
        xfrm.apply_to(m, reversible=True)

        # Note: The strip_var_bounds transformation is used to remove the
        # bounds from all the variables in the model. This helps the model
        # converge if the variables are very close to one of their bounds.
        # However, an alternate method would be to specify an appropriate
        # value for the solver option bound_push during model initialization
        # such that the minimum absolute distance from the initial point to
        # the bounds is reduced beyond the default value of 0.01.

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, model):

        assert hasattr(model.fs.unit, "inlet")
        assert len(model.fs.unit.inlet.vars) == 4
        assert hasattr(model.fs.unit.inlet, "flow_mol")
        assert hasattr(model.fs.unit.inlet, "mole_frac_comp")
        assert hasattr(model.fs.unit.inlet, "temperature")
        assert hasattr(model.fs.unit.inlet, "pressure")

        assert hasattr(model.fs.unit, "reflux")
        assert len(model.fs.unit.reflux.vars) == 4
        assert hasattr(model.fs.unit.reflux, "flow_mol")
        assert hasattr(model.fs.unit.reflux, "mole_frac_comp")
        assert hasattr(model.fs.unit.reflux, "temperature")
        assert hasattr(model.fs.unit.reflux, "pressure")

        assert hasattr(model.fs.unit, "vapor_outlet")
        assert len(model.fs.unit.vapor_outlet.vars) == 4
        assert hasattr(model.fs.unit.vapor_outlet, "flow_mol")
        assert hasattr(model.fs.unit.vapor_outlet, "mole_frac_comp")
        assert hasattr(model.fs.unit.vapor_outlet, "temperature")
        assert hasattr(model.fs.unit.vapor_outlet, "pressure")

        assert isinstance(model.fs.unit.unit_material_balance, Constraint)
        assert isinstance(model.fs.unit.unit_enthalpy_balance, Constraint)
        assert isinstance(model.fs.unit.unit_temperature_equality, Constraint)
        assert isinstance(model.fs.unit.unit_pressure_balance, Constraint)
        assert isinstance(model.fs.unit.zero_flow_param, Param)

        assert number_variables(model.fs.unit) == 55
        assert number_total_constraints(model.fs.unit) == 49
        assert number_unused_variables(model.fs.unit) == 0

    @pytest.mark.component
    def test_units(self, model):
        assert_units_consistent(model)
        assert_units_equivalent(model.fs.unit.heat_duty[0], units.W)

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model)

    # @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert pytest.approx(0.108117, rel=1e-5) == value(
            model.fs.unit.reflux.flow_mol[0]
        )
        assert pytest.approx(0, abs=1e-3) == value(
            model.fs.unit.reflux.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0, abs=1e-3) == value(
            model.fs.unit.reflux.mole_frac_comp[0, "MEA"]
        )
        assert pytest.approx(1, rel=1e-3) == value(
            model.fs.unit.reflux.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(184360, rel=1e-5) == value(
            model.fs.unit.reflux.pressure[0]
        )
        assert pytest.approx(303.376, rel=1e-5) == value(
            model.fs.unit.reflux.temperature[0]
        )

        assert pytest.approx(1.00358, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.flow_mol[0]
        )
        assert pytest.approx(0.976580, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.0234192, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(184360, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.pressure[0]
        )
        assert pytest.approx(303.376, rel=1e-5) == value(
            model.fs.unit.vapor_outlet.temperature[0]
        )

        assert pytest.approx(-6264, rel=1e-5) == value(model.fs.unit.heat_duty[0])

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, model):
        assert (
            abs(
                value(
                    model.fs.unit.inlet.flow_mol[0]
                    - model.fs.unit.reflux.flow_mol[0]
                    - model.fs.unit.vapor_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    model.fs.unit.inlet.flow_mol[0]
                    * model.fs.unit.inlet.mole_frac_comp[0, "CO2"]
                    - model.fs.unit.reflux.flow_mol[0]
                    * model.fs.unit.reflux.mole_frac_comp[0, "CO2"]
                    - model.fs.unit.vapor_outlet.flow_mol[0]
                    * model.fs.unit.vapor_outlet.mole_frac_comp[0, "CO2"]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    model.fs.unit.inlet.flow_mol[0]
                    * model.fs.unit.inlet.mole_frac_comp[0, "H2O"]
                    - model.fs.unit.reflux.flow_mol[0]
                    * model.fs.unit.reflux.mole_frac_comp[0, "H2O"]
                    - model.fs.unit.vapor_outlet.flow_mol[0]
                    * model.fs.unit.vapor_outlet.mole_frac_comp[0, "H2O"]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    model.fs.unit.reflux.flow_mol[0]
                    * model.fs.unit.reflux.mole_frac_comp[0, "MEA"]
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    model.fs.unit.vapor_phase.properties_in[0]._enthalpy_flow_term[
                        "Vap"
                    ]
                    - model.fs.unit.vapor_phase.properties_out[0]._enthalpy_flow_term[
                        "Vap"
                    ]
                    - model.fs.unit.liquid_phase[0]._enthalpy_flow_term["Liq"]
                    + model.fs.unit.heat_duty[0]
                )
            )
            <= 1e-6
        )

    @pytest.mark.component
    def test_scaling(self, model):
        iscale.set_scaling_factor(
            model.fs.unit.vapor_phase.properties_out[0].fug_phase_comp["Vap", "CO2"],
            1e-5,
        )
        iscale.set_scaling_factor(
            model.fs.unit.vapor_phase.properties_out[0].fug_phase_comp["Vap", "H2O"],
            1e-3,
        )

        iscale.calculate_scaling_factors(model.fs.unit)

        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.unit_material_balance[0, "CO2"]
            )
            == 1
        )
        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.unit_material_balance[0, "H2O"]
            )
            == 1
        )
        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.unit_material_balance[0, "MEA"]
            )
            == 1e8
        )

        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.unit_phase_equilibrium[0, "CO2"]
            )
            == 1e-5
        )
        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.unit_phase_equilibrium[0, "H2O"]
            )
            == 1e-3
        )

        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.unit_temperature_equality[0]
            )
            == 1e-2
        )

        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.unit_enthalpy_balance[0]
            )
            == 1
        )

        assert (
            iscale.get_constraint_transform_applied_scaling_factor(
                model.fs.unit.unit_pressure_balance[0]
            )
            == 1e-5
        )

    @pytest.mark.component
    def test_initialization_error_dof(self, model):
        model.fs.unit.reflux.flow_mol[0].fix(100)

        with pytest.raises(InitializationError):
            model.fs.unit.initialize()
