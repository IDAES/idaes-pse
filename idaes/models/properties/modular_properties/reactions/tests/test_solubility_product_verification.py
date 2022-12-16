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
Verification tests for solubility product
"""

import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    units as pyunits,
    value,
)

# Import IDAES cores
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
    ConcentrationForm,
)
from idaes.models.unit_models import EquilibriumReactor
from idaes.core import (
    Component,
    EnergyBalanceType,
    FlowsheetBlock,
    LiquidPhase,
    SolidPhase,
)
from idaes.core.solvers import get_solver

from idaes.models.properties.modular_properties.state_definitions import FpcTP
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.reactions.equilibrium_forms import (
    solubility_product,
)
from idaes.models.properties.modular_properties.reactions.equilibrium_constant import (
    ConstantKeq,
)


solver = get_solver()


def dummy_h(b, *args, **kwargs):
    return b.temperature


thermo_config = {
    # Specifying components
    "components": {
        "H2O": {"type": Component, "enth_mol_liq_comp": dummy_h},
        "Na+": {"type": Component, "enth_mol_liq_comp": dummy_h},
        "Cl-": {"type": Component, "enth_mol_liq_comp": dummy_h},
        "NaCl": {"type": Component, "enth_mol_sol_comp": dummy_h},
    },
    # Specifying phases
    "phases": {
        "Liq": {
            "type": LiquidPhase,
            "equation_of_state": Ideal,
            "component_list": ["H2O", "Na+", "Cl-"],
        },
        "Sol": {
            "type": SolidPhase,
            "equation_of_state": Ideal,
            "component_list": ["NaCl"],
        },
    },
    # Set base units of measurement
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    # Specifying state definition
    "state_definition": FpcTP,
    "state_bounds": {
        "flow_mol_phase_comp": (0, 100, 1000, pyunits.mol / pyunits.s),
        "temperature": (273.15, 300, 450, pyunits.K),
        "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
    },
    "pressure_ref": (1e5, pyunits.Pa),
    "temperature_ref": (300, pyunits.K),
}

# Solubility of NaCl in H2O is 360 g/L, or 6.16 mol/L
# 1L of H2O is 55.56 mol, results in Ksp = 8.235e-3 on mole frac basis
rxn_config = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "equilibrium_reactions": {
        "S1": {
            "stoichiometry": {
                ("Liq", "Na+"): 1,
                ("Liq", "Cl-"): 1,
                ("Sol", "NaCl"): -1,
            },
            "equilibrium_constant": ConstantKeq,
            "equilibrium_form": solubility_product,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {"k_eq_ref": (8.235e-3, None)},
        }
    },
}


class TestSingleState(object):
    @pytest.fixture
    def model(self):
        m = ConcreteModel()

        # Add a test thermo package for validation
        m.pparams = GenericParameterBlock(**thermo_config)
        m.rparams = GenericReactionParameterBlock(
            property_package=m.pparams, **rxn_config
        )

        m.state = m.pparams.build_state_block([0], defined_state=False)

        m.rxn = m.rparams.build_reaction_block(
            [0], state_block=m.state, has_equilibrium=True
        )

        return m

    @pytest.mark.component
    def test_saturated(self, model):
        assert model.state[0].phase_component_set == [
            ("Liq", "H2O"),
            ("Liq", "Na+"),
            ("Liq", "Cl-"),
            ("Sol", "NaCl"),
        ]

        model.state[0].temperature.fix(298.15)
        model.state[0].pressure.fix(101325)

        # Solve for saturated states (i.e. solubility product applied) by
        # fixing the flowrate of solids to a positive value.
        # Then start with a large amount of Na+ and solve for the flowrate of
        # Cl- which satisfies the solubility product (ignore electroneutrality)
        model.state[0].flow_mol_phase_comp["Liq", "H2O"].fix(55.56)
        model.state[0].flow_mol_phase_comp["Liq", "Na+"].fix(20)
        model.state[0].flow_mol_phase_comp["Liq", "Cl-"].set_value(2.5)
        model.state[0].flow_mol_phase_comp["Sol", "NaCl"].fix(1)

        for i in range(20, 2, -1):
            model.state[0].flow_mol_phase_comp["Liq", "Na+"].fix(i)

            results = solver.solve(model)

            assert check_optimal_termination(results)

            assert pytest.approx(8.235e-3, abs=1e-8) == value(
                model.state[0].mole_frac_phase_comp["Liq", "Na+"]
                * model.state[0].mole_frac_phase_comp["Liq", "Cl-"]
            )

    @pytest.mark.component
    def test_subsaturated(self, model):
        assert model.state[0].phase_component_set == [
            ("Liq", "H2O"),
            ("Liq", "Na+"),
            ("Liq", "Cl-"),
            ("Sol", "NaCl"),
        ]

        model.state[0].temperature.fix(298.15)
        model.state[0].pressure.fix(101325)

        # For subsaturated systems, fix the Na and Cl flows and check that no
        # solid is formed.
        model.state[0].flow_mol_phase_comp["Liq", "H2O"].fix(55.56)
        model.state[0].flow_mol_phase_comp["Liq", "Na+"].fix(0)
        model.state[0].flow_mol_phase_comp["Liq", "Cl-"].fix(0)
        model.state[0].flow_mol_phase_comp["Sol", "NaCl"].set_value(0)

        for i in range(11):
            for j in range(11):
                if i * j < 40:
                    # Zeroes result in evaluation errors, so use a small number
                    if i == 0:
                        i = 1e-8
                    if j == 0:
                        j = 1e-8

                    model.state[0].flow_mol_phase_comp["Liq", "Na+"].fix(i)
                    model.state[0].flow_mol_phase_comp["Liq", "Cl-"].fix(j)

                    results = solver.solve(model)

                    assert check_optimal_termination(results)

                    assert pytest.approx(0, abs=1e-5) == value(
                        model.state[0].flow_mol_phase_comp["Sol", "NaCl"]
                    )


class TestUnit(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # Add a test thermo package for validation
        m.fs.pparams = GenericParameterBlock(**thermo_config)
        m.fs.rparams = GenericReactionParameterBlock(
            property_package=m.fs.pparams, **rxn_config
        )

        # Don't include energy balances, as the test doesn't have a proper
        # enthalpy model. Fix outlet T instead.
        m.fs.R101 = EquilibriumReactor(
            property_package=m.fs.pparams,
            reaction_package=m.fs.rparams,
            has_equilibrium_reactions=True,
            has_rate_reactions=False,
            energy_balance_type=EnergyBalanceType.none,
        )

        m.fs.R101.inlet.flow_mol_phase_comp[0, "Liq", "H2O"].fix(55.56)
        m.fs.R101.inlet.flow_mol_phase_comp[0, "Liq", "Na+"].fix(1e-8)
        m.fs.R101.inlet.flow_mol_phase_comp[0, "Liq", "Cl-"].fix(1e-8)
        m.fs.R101.inlet.flow_mol_phase_comp[0, "Sol", "NaCl"].fix(1e-8)
        m.fs.R101.inlet.temperature[0].fix(298.15)
        m.fs.R101.inlet.pressure[0].fix(101325)

        m.fs.R101.outlet.temperature[0].fix(298.15)

        return m

    @pytest.mark.integration
    def test_subsaturated(self, model):
        for i in range(0, 6):
            # Values of zero cause probelms, use a small number instead
            if i == 0:
                i == 1e-8

            model.fs.R101.inlet.flow_mol_phase_comp[0, "Sol", "NaCl"].fix(i)

            results = solver.solve(model)

            assert check_optimal_termination(results)

            assert pytest.approx(0, abs=1.1e-6) == value(
                model.fs.R101.outlet.flow_mol_phase_comp[0, "Sol", "NaCl"]
            )
            assert pytest.approx(i, rel=1e-6, abs=1e-6) == value(
                model.fs.R101.outlet.flow_mol_phase_comp[0, "Liq", "Na+"]
            )
            assert pytest.approx(i, rel=1e-6, abs=1e-6) == value(
                model.fs.R101.outlet.flow_mol_phase_comp[0, "Liq", "Cl-"]
            )

    @pytest.mark.integration
    def test_saturated(self, model):
        model.fs.R101.inlet.flow_mol_phase_comp[0, "Sol", "NaCl"].fix(6.16)

        results = solver.solve(model)

        assert check_optimal_termination(results)

        assert pytest.approx(0, abs=1.1e-3) == value(
            model.fs.R101.outlet.flow_mol_phase_comp[0, "Sol", "NaCl"]
        )
        assert pytest.approx(6.158968, rel=1e-5) == value(
            model.fs.R101.outlet.flow_mol_phase_comp[0, "Liq", "Na+"]
        )
        assert pytest.approx(6.158968, rel=1e-5) == value(
            model.fs.R101.outlet.flow_mol_phase_comp[0, "Liq", "Cl-"]
        )

    @pytest.mark.integration
    def test_supersaturated(self, model):
        for i in range(7, 21):
            model.fs.R101.inlet.flow_mol_phase_comp[0, "Sol", "NaCl"].fix(i)

            results = solver.solve(model)

            assert check_optimal_termination(results)

            assert pytest.approx(i, rel=1e-5) == value(
                model.fs.R101.outlet.flow_mol_phase_comp[0, "Sol", "NaCl"]
                + model.fs.R101.outlet.flow_mol_phase_comp[0, "Liq", "Na+"]
            )
            assert pytest.approx(6.159876, rel=1e-5) == value(
                model.fs.R101.outlet.flow_mol_phase_comp[0, "Liq", "Na+"]
            )
            assert pytest.approx(6.159876, rel=1e-5) == value(
                model.fs.R101.outlet.flow_mol_phase_comp[0, "Liq", "Cl-"]
            )
