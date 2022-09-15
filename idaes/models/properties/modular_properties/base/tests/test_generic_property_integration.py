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
Tests for generic property package core code

Author: Andrew Lee
"""
import pytest

# Import Pyomo units
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Expression,
    units as pyunits,
    value,
)

# Import IDAES cores
from idaes.core import (
    FlowsheetBlock,
    VaporPhase,
    Component,
    MaterialBalanceType,
    ControlVolume1DBlock,
)

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
    StateIndex,
)
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.pure import RPP4
from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn
from idaes.models.properties.modular_properties.reactions.equilibrium_constant import (
    van_t_hoff,
)
from idaes.models.properties.modular_properties.reactions.equilibrium_forms import (
    power_law_equil,
)
from idaes.models.properties.modular_properties.base.utility import ConcentrationForm

from idaes.models.unit_models import Heater
from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.core.solvers import get_solver


configuration = {
    # Specifying components
    "components": {
        "a": {
            "type": Component,
            "enth_mol_ig_comp": RPP4,
            "parameter_data": {
                "mw": (78.1136e-3, pyunits.kg / pyunits.mol),
                "pressure_crit": (48.9e5, pyunits.Pa),
                "temperature_crit": (562.2, pyunits.K),
                "cp_mol_ig_comp_coeff": {
                    "A": (-40, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (0, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (0, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (0, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "enth_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),
            },
        },
        "b": {
            "type": Component,
            "enth_mol_ig_comp": RPP4,
            "parameter_data": {
                "mw": (78.1136e-3, pyunits.kg / pyunits.mol),
                "pressure_crit": (48.9e5, pyunits.Pa),
                "temperature_crit": (562.2, pyunits.K),
                "cp_mol_ig_comp_coeff": {
                    "A": (-40, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (0, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (0, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (0, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "enth_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),
            },
        },
        "c": {
            "type": Component,
            "enth_mol_ig_comp": RPP4,
            "parameter_data": {
                "mw": (78.1136e-3, pyunits.kg / pyunits.mol),
                "pressure_crit": (48.9e5, pyunits.Pa),
                "temperature_crit": (562.2, pyunits.K),
                "cp_mol_ig_comp_coeff": {
                    "A": (-40, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (0, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (0, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (0, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "enth_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),
            },
        },
        "d": {
            "type": Component,
            "enth_mol_ig_comp": RPP4,
            "parameter_data": {
                "mw": (78.1136e-3, pyunits.kg / pyunits.mol),
                "pressure_crit": (48.9e5, pyunits.Pa),
                "temperature_crit": (562.2, pyunits.K),
                "cp_mol_ig_comp_coeff": {
                    "A": (-40, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (0, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (0, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (0, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "enth_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),
            },
        },
    },
    # Specifying phases
    "phases": {"Vap": {"type": VaporPhase, "equation_of_state": Ideal}},
    # Set base units of measurement
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {
        "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
        "temperature": (273.15, 300, 450, pyunits.K),
        "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
    },
    "state_components": StateIndex.true,  # need this for inherent rxn terms in material balances
    "pressure_ref": (1e5, pyunits.Pa),
    "temperature_ref": (300, pyunits.K),
    "inherent_reactions": {
        "e1": {
            "stoichiometry": {("Vap", "a"): -1, ("Vap", "b"): 1},
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {"dh_rxn_ref": -20000, "k_eq_ref": 2, "T_eq_ref": 350},
        }
    },
}


class TestInherentReactions(object):
    @pytest.fixture()
    def frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.params = GenericParameterBlock(**configuration)

        return m

    @pytest.mark.unit
    def test_build(self, frame):
        frame.fs.props = frame.fs.params.build_state_block([1], defined_state=False)

        assert isinstance(frame.fs.props[1].inherent_equilibrium_constraint, Constraint)
        assert len(frame.fs.props[1].inherent_equilibrium_constraint) == 1

        assert isinstance(frame.fs.props[1].dh_rxn, Expression)
        assert len(frame.fs.props[1].dh_rxn) == 1

    @pytest.mark.component
    def test_heater_w_inherent_rxns_comp_phase(self, frame):
        frame.fs.H101 = Heater(
            property_package=frame.fs.params,
            material_balance_type=MaterialBalanceType.componentPhase,
        )

        frame.fs.H101.inlet.flow_mol.fix(100)
        frame.fs.H101.inlet.mole_frac_comp.fix(0.25)
        frame.fs.H101.inlet.temperature.fix(350)
        frame.fs.H101.inlet.pressure.fix(101325)

        frame.fs.H101.heat_duty.fix(0)

        assert (degrees_of_freedom(frame)) == 0

        frame.fs.H101.initialize()

        solver = get_solver()

        results = solver.solve(frame)

        assert check_optimal_termination(results)

        assert value(frame.fs.H101.control_volume.properties_out[0].k_eq["e1"]) == 2

        assert value(
            frame.fs.H101.control_volume.properties_out[0].k_eq["e1"]
        ) == pytest.approx(
            value(
                frame.fs.H101.outlet.mole_frac_comp[0, "b"]
                / frame.fs.H101.outlet.mole_frac_comp[0, "a"]
            ),
            rel=1e-5,
        )

        assert value(frame.fs.H101.outlet.mole_frac_comp[0, "a"]) == pytest.approx(
            1 / 6, rel=1e-5
        )
        assert value(frame.fs.H101.outlet.mole_frac_comp[0, "b"]) == pytest.approx(
            1 / 3, rel=1e-5
        )
        assert value(frame.fs.H101.outlet.mole_frac_comp[0, "c"]) == pytest.approx(
            1 / 4, rel=1e-5
        )
        assert value(frame.fs.H101.outlet.mole_frac_comp[0, "d"]) == pytest.approx(
            1 / 4, rel=1e-5
        )

        assert value(frame.fs.H101.outlet.flow_mol[0]) == pytest.approx(100, rel=1e-5)
        assert value(frame.fs.H101.outlet.temperature[0]) == pytest.approx(
            350, rel=1e-5
        )
        assert value(frame.fs.H101.outlet.pressure[0]) == pytest.approx(
            101325, rel=1e-5
        )

    @pytest.mark.component
    def test_heater_w_inherent_rxns_comp_total(self, frame):
        frame.fs.H101 = Heater(
            property_package=frame.fs.params,
            material_balance_type=MaterialBalanceType.componentTotal,
        )

        frame.fs.H101.inlet.flow_mol.fix(100)
        frame.fs.H101.inlet.mole_frac_comp.fix(0.25)
        frame.fs.H101.inlet.temperature.fix(350)
        frame.fs.H101.inlet.pressure.fix(101325)

        frame.fs.H101.heat_duty.fix(0)

        assert (degrees_of_freedom(frame)) == 0

        frame.fs.H101.initialize()

        solver = get_solver()

        results = solver.solve(frame)

        assert check_optimal_termination(results)

        assert value(frame.fs.H101.control_volume.properties_out[0].k_eq["e1"]) == 2

        assert value(
            frame.fs.H101.control_volume.properties_out[0].k_eq["e1"]
        ) == pytest.approx(
            value(
                frame.fs.H101.outlet.mole_frac_comp[0, "b"]
                / frame.fs.H101.outlet.mole_frac_comp[0, "a"]
            ),
            rel=1e-5,
        )

        assert value(frame.fs.H101.outlet.mole_frac_comp[0, "a"]) == pytest.approx(
            1 / 6, rel=1e-5
        )
        assert value(frame.fs.H101.outlet.mole_frac_comp[0, "b"]) == pytest.approx(
            1 / 3, rel=1e-5
        )
        assert value(frame.fs.H101.outlet.mole_frac_comp[0, "c"]) == pytest.approx(
            1 / 4, rel=1e-5
        )
        assert value(frame.fs.H101.outlet.mole_frac_comp[0, "d"]) == pytest.approx(
            1 / 4, rel=1e-5
        )

        assert value(frame.fs.H101.outlet.flow_mol[0]) == pytest.approx(100, rel=1e-5)
        assert value(frame.fs.H101.outlet.temperature[0]) == pytest.approx(
            350, rel=1e-5
        )
        assert value(frame.fs.H101.outlet.pressure[0]) == pytest.approx(
            101325, rel=1e-5
        )

    @pytest.mark.component
    def test_CV1D_w_inherent_rxns_comp_phase(self, frame):
        frame.fs.cv = ControlVolume1DBlock(
            property_package=frame.fs.params,
            transformation_method="dae.finite_difference",
            transformation_scheme="BACKWARD",
            finite_elements=2,
        )

        frame.fs.cv.add_geometry()

        frame.fs.cv.add_state_blocks(has_phase_equilibrium=False)

        frame.fs.cv.add_material_balances(
            balance_type=MaterialBalanceType.componentPhase
        )

        frame.fs.cv.add_energy_balances()
        frame.fs.cv.add_momentum_balances()
        frame.fs.cv.apply_transformation()

        frame.fs.cv.properties[0, 0].flow_mol.fix(100)
        frame.fs.cv.properties[0, 0].mole_frac_comp.fix(0.25)
        frame.fs.cv.properties[0, 0].temperature.fix(350)
        frame.fs.cv.properties[0, 0].pressure.fix(101325)

        frame.fs.cv.area.fix(1)
        frame.fs.cv.length.fix(1)

        assert (degrees_of_freedom(frame)) == 0

        frame.fs.cv.initialize()

        solver = get_solver()

        results = solver.solve(frame)

        assert check_optimal_termination(results)

        assert pytest.approx(2, rel=1e-8) == value(
            frame.fs.cv.properties[0, 1].k_eq["e1"]
        )

        assert value(frame.fs.cv.properties[0, 1].k_eq["e1"]) == pytest.approx(
            value(
                frame.fs.cv.properties[0, 1].mole_frac_comp["b"]
                / frame.fs.cv.properties[0, 1].mole_frac_comp["a"]
            ),
            rel=1e-5,
        )

        assert value(frame.fs.cv.properties[0, 1].mole_frac_comp["a"]) == pytest.approx(
            1 / 6, rel=1e-5
        )
        assert value(frame.fs.cv.properties[0, 1].mole_frac_comp["b"]) == pytest.approx(
            1 / 3, rel=1e-5
        )
        assert value(frame.fs.cv.properties[0, 1].mole_frac_comp["c"]) == pytest.approx(
            1 / 4, rel=1e-5
        )
        assert value(frame.fs.cv.properties[0, 1].mole_frac_comp["d"]) == pytest.approx(
            1 / 4, rel=1e-5
        )

        assert value(frame.fs.cv.properties[0, 1].flow_mol) == pytest.approx(
            100, rel=1e-5
        )
        assert value(frame.fs.cv.properties[0, 1].temperature) == pytest.approx(
            350, rel=1e-5
        )
        assert value(frame.fs.cv.properties[0, 1].pressure) == pytest.approx(
            101325, rel=1e-5
        )

    @pytest.mark.component
    def test_CV1D_w_inherent_rxns_comp_total(self, frame):
        frame.fs.cv = ControlVolume1DBlock(
            property_package=frame.fs.params,
            transformation_method="dae.finite_difference",
            transformation_scheme="BACKWARD",
            finite_elements=2,
        )

        frame.fs.cv.add_geometry()

        frame.fs.cv.add_state_blocks(has_phase_equilibrium=False)

        frame.fs.cv.add_material_balances(
            balance_type=MaterialBalanceType.componentTotal
        )

        frame.fs.cv.add_energy_balances()
        frame.fs.cv.add_momentum_balances()
        frame.fs.cv.apply_transformation()

        frame.fs.cv.properties[0, 0].flow_mol.fix(100)
        frame.fs.cv.properties[0, 0].mole_frac_comp.fix(0.25)
        frame.fs.cv.properties[0, 0].temperature.fix(350)
        frame.fs.cv.properties[0, 0].pressure.fix(101325)

        frame.fs.cv.area.fix(1)
        frame.fs.cv.length.fix(1)

        assert (degrees_of_freedom(frame)) == 0

        frame.fs.cv.initialize()

        solver = get_solver()

        results = solver.solve(frame)

        assert check_optimal_termination(results)

        assert pytest.approx(2, rel=1e-8) == value(
            frame.fs.cv.properties[0, 1].k_eq["e1"]
        )

        assert value(frame.fs.cv.properties[0, 1].k_eq["e1"]) == pytest.approx(
            value(
                frame.fs.cv.properties[0, 1].mole_frac_comp["b"]
                / frame.fs.cv.properties[0, 1].mole_frac_comp["a"]
            ),
            rel=1e-5,
        )

        assert value(frame.fs.cv.properties[0, 1].mole_frac_comp["a"]) == pytest.approx(
            1 / 6, rel=1e-5
        )
        assert value(frame.fs.cv.properties[0, 1].mole_frac_comp["b"]) == pytest.approx(
            1 / 3, rel=1e-5
        )
        assert value(frame.fs.cv.properties[0, 1].mole_frac_comp["c"]) == pytest.approx(
            1 / 4, rel=1e-5
        )
        assert value(frame.fs.cv.properties[0, 1].mole_frac_comp["d"]) == pytest.approx(
            1 / 4, rel=1e-5
        )

        assert value(frame.fs.cv.properties[0, 1].flow_mol) == pytest.approx(
            100, rel=1e-5
        )
        assert value(frame.fs.cv.properties[0, 1].temperature) == pytest.approx(
            350, rel=1e-5
        )
        assert value(frame.fs.cv.properties[0, 1].pressure) == pytest.approx(
            101325, rel=1e-5
        )
