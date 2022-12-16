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
Tests for ControlVolumeBlockData.

Author: Andrew Lee
"""
import pytest

from pyomo.environ import check_optimal_termination, ConcreteModel, value, Var, units
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.models.unit_models.plug_flow_reactor import PFR

from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)
from idaes.models.properties.examples.saponification_reactions import (
    SaponificationReactionParameterBlock,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
    number_derivative_variables,
)
from idaes.core.util.testing import (
    PhysicalParameterTestBlock,
    ReactionParameterTestBlock,
    initialization_tester,
)
from idaes.core.solvers import get_solver


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.component
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()
    m.fs.reactions = ReactionParameterTestBlock(property_package=m.fs.properties)

    m.fs.unit = PFR(property_package=m.fs.properties, reaction_package=m.fs.reactions)

    # Check unit config arguments
    assert len(m.fs.unit.config) == 19

    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_heat_transfer
    assert not m.fs.unit.config.has_pressure_change
    assert not m.fs.unit.config.has_equilibrium_reactions
    assert not m.fs.unit.config.has_phase_equilibrium
    assert not m.fs.unit.config.has_heat_of_reaction
    assert m.fs.unit.config.property_package is m.fs.properties
    assert m.fs.unit.config.reaction_package is m.fs.reactions

    assert m.fs.unit.config.length_domain_set == [0.0, 1.0]
    assert m.fs.unit.config.transformation_method == "dae.finite_difference"
    assert m.fs.unit.config.transformation_scheme == "BACKWARD"
    assert m.fs.unit.config.finite_elements == 20
    assert m.fs.unit.config.collocation_points == 3


# -----------------------------------------------------------------------------
class TestSaponification(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()
        m.fs.reactions = SaponificationReactionParameterBlock(
            property_package=m.fs.properties
        )

        m.fs.unit = PFR(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_equilibrium_reactions=False,
            has_heat_transfer=True,
            has_heat_of_reaction=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_vol.fix(1.0)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.inlet.temperature.fix(303.15)
        m.fs.unit.inlet.pressure.fix(101325.0)

        m.fs.unit.length.fix(0.5)
        m.fs.unit.area.fix(0.1)

        m.fs.unit.heat_duty.fix(0)
        m.fs.unit.deltaP.fix(0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, sapon):
        assert hasattr(sapon.fs.unit, "inlet")
        assert len(sapon.fs.unit.inlet.vars) == 4
        assert hasattr(sapon.fs.unit.inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.inlet, "temperature")
        assert hasattr(sapon.fs.unit.inlet, "pressure")

        assert hasattr(sapon.fs.unit, "outlet")
        assert len(sapon.fs.unit.outlet.vars) == 4
        assert hasattr(sapon.fs.unit.outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.outlet, "temperature")
        assert hasattr(sapon.fs.unit.outlet, "pressure")

        assert isinstance(sapon.fs.unit.area, Var)
        assert isinstance(sapon.fs.unit.length, Var)
        assert isinstance(sapon.fs.unit.volume, Var)
        assert hasattr(sapon.fs.unit, "performance_eqn")
        assert hasattr(sapon.fs.unit.control_volume, "heat")
        assert hasattr(sapon.fs.unit, "heat_duty")
        assert hasattr(sapon.fs.unit, "deltaP")

        assert number_variables(sapon) == 654
        assert number_total_constraints(sapon) == 590
        assert number_unused_variables(sapon) == 14
        assert number_derivative_variables(sapon) == 0

    @pytest.mark.component
    def test_units(self, sapon):
        assert_units_consistent(sapon)
        assert_units_equivalent(sapon.fs.unit.volume, units.m**3)
        assert_units_equivalent(sapon.fs.unit.length, units.m)
        assert_units_equivalent(sapon.fs.unit.area, units.m**2)
        assert_units_equivalent(sapon.fs.unit.heat_duty[0, 0], units.W / units.m)
        assert_units_equivalent(sapon.fs.unit.deltaP[0, 0], units.Pa / units.m)

    @pytest.mark.unit
    def test_dof(self, sapon):
        assert degrees_of_freedom(sapon) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, sapon):
        initialization_tester(sapon)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, sapon):
        results = solver.solve(sapon)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, sapon):
        assert (
            pytest.approx(101325.0, abs=1e-2) == sapon.fs.unit.outlet.pressure[0].value
        )
        assert (
            pytest.approx(303.6, abs=1e-2) == sapon.fs.unit.outlet.temperature[0].value
        )
        assert pytest.approx(62.29, abs=1e-2) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, sapon):
        assert (
            abs(
                value(
                    sapon.fs.unit.inlet.flow_vol[0] - sapon.fs.unit.outlet.flow_vol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    sapon.fs.unit.inlet.flow_vol[0]
                    * sum(
                        sapon.fs.unit.inlet.conc_mol_comp[0, j]
                        for j in sapon.fs.properties.component_list
                    )
                    - sapon.fs.unit.outlet.flow_vol[0]
                    * sum(
                        sapon.fs.unit.outlet.conc_mol_comp[0, j]
                        for j in sapon.fs.properties.component_list
                    )
                )
            )
            <= 1e-6
        )

        hrxn = 0
        for x in sapon.fs.unit.control_volume.length_domain:
            if x != 0:
                hrxn += value(
                    sapon.fs.unit.control_volume.heat_of_reaction[0, x]
                    * (x - sapon.fs.unit.control_volume.length_domain.prev(x))
                    * sapon.fs.unit.control_volume.length
                )
        assert pytest.approx(1847000, abs=1e3) == hrxn
        assert (
            abs(
                value(
                    (
                        sapon.fs.unit.inlet.flow_vol[0]
                        * sapon.fs.properties.dens_mol
                        * sapon.fs.properties.cp_mol
                        * (
                            sapon.fs.unit.inlet.temperature[0]
                            - sapon.fs.properties.temperature_ref
                        )
                    )
                    - (
                        sapon.fs.unit.outlet.flow_vol[0]
                        * sapon.fs.properties.dens_mol
                        * sapon.fs.properties.cp_mol
                        * (
                            sapon.fs.unit.outlet.temperature[0]
                            - sapon.fs.properties.temperature_ref
                        )
                    )
                    + hrxn
                )
            )
            <= 1e-6
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, sapon):
        perf_dict = sapon.fs.unit._get_performance_contents()

        assert perf_dict == {"vars": {"Area": sapon.fs.unit.area}}
