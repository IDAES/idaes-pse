##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Tests for rate forms
"""

import pytest

from pyomo.environ import \
    ConcreteModel, units as pyunits, value, SolverStatus, TerminationCondition

from idaes.generic_models.properties.core.generic.generic_property import \
    GenericParameterBlock
from idaes.generic_models.properties.core.generic.generic_reaction import \
    GenericReactionParameterBlock, ConcentrationForm

# Import IDAES cores
from idaes.core import LiquidPhase, SolidPhase, Component
from idaes.core.util import get_solver

from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.generic_models.properties.core.reactions.equilibrium_forms import \
    solubility_product
from idaes.generic_models.properties.core.reactions.equilibrium_constant import \
    ConstantKeq


solver = get_solver()


def dummy_h(b, **kwargs):
    return b.temperature


thermo_config = {
    # Specifying components
    "components": {
        'A': {"type": Component,
              "enth_mol_liq_comp": dummy_h},
        'B': {"type": Component,
              "enth_mol_liq_comp": dummy_h},
        'C': {"type": Component,
              "enth_mol_sol_comp": dummy_h}},

    # Specifying phases
    "phases":  {'Liq': {"type": LiquidPhase,
                        "equation_of_state": Ideal,
                        "component_list": ["A", "B"]},
                'Sol': {"type": SolidPhase,
                        "equation_of_state": Ideal,
                        "component_list": ["C"]}},

    # Set base units of measurement
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},

    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {"flow_mol": (0, 100, 1000, pyunits.mol/pyunits.s),
                     "temperature": (273.15, 300, 450, pyunits.K),
                     "pressure": (5e4, 1e5, 1e6, pyunits.Pa)},
    "pressure_ref": (1e5, pyunits.Pa),
    "temperature_ref": (300, pyunits.K)}


rxn_config = {
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "equilibrium_reactions": {
        "S2": {"stoichiometry": {("Liq", "A"): -1,
                                 ("Liq", "B"): -1,
                                 ("Sol", "C"): 1},
               "equilibrium_constant": ConstantKeq,
               "equilibrium_form": solubility_product,
               "concentration_form": ConcentrationForm.moleFraction,
               "parameter_data": {
                   "k_eq_ref": (0.1, None)}}}}


@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()

    # Add a test thermo package for validation
    m.pparams = GenericParameterBlock(default=thermo_config)
    m.rparams = GenericReactionParameterBlock(default={
        "property_package": m.pparams, **rxn_config})

    m.state = m.pparams.build_state_block(
        [0], default={"defined_state": False})

    m.rxn = m.rparams.build_reaction_block(
        [0], default={"state_block": m.state, "has_equilibrium": True})

    return m


@pytest.mark.integration
def test_saturated(model):
    assert model.state[0].phase_component_set == [
        ("Liq", "A"), ("Liq", "B"), ("Sol", "C")]

    model.state[0].mole_frac_comp["C"].set_value(0.18)

    # Set liquid phase compositions
    for i in range(15):
        A = 0.15+0.05*i
        B = 0.1/A
        C = 1-A-B

        model.state[0].mole_frac_comp["A"].fix(A)
        model.state[0].mole_frac_comp["B"].fix(B)

        results = solver.solve(model)

        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

        assert pytest.approx(C, rel=1e-5) == value(
            model.state[0].mole_frac_comp["C"])


@pytest.mark.integration
def test_subsaturated(model):
    # Set liquid phase compositions
    for i in range(12):
        A = 1-0.01*i
        B = 1-A

        model.state[0].mole_frac_comp["A"].fix(A)
        model.state[0].mole_frac_comp["B"].fix(B)

        results = solver.solve(model)

        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

        assert pytest.approx(0, abs=1e-5) == value(
            model.state[0].mole_frac_comp["C"])

    for i in range(12):
        B = 1-0.01*i
        A = 1-B

        model.state[0].mole_frac_comp["A"].fix(A)
        model.state[0].mole_frac_comp["B"].fix(B)

        results = solver.solve(model, tee=True)

        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

        assert pytest.approx(0, abs=1e-5) == value(
            model.state[0].mole_frac_comp["C"])
