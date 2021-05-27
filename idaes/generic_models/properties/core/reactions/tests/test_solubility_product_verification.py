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
    Block, ConcreteModel, Param, Var, units as pyunits, value

from idaes.generic_models.properties.core.generic.generic_property import \
    GenericParameterBlock
from idaes.generic_models.properties.core.generic.generic_reaction import \
    GenericReactionParameterBlock, ConcentrationForm

# Import IDAES cores
from idaes.core import LiquidPhase, SolidPhase, Component

from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.generic_models.properties.core.reactions.equilibrium_forms import \
    solubility_product


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


@pytest.fixture
def model():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = GenericParameterBlock(default=thermo_config)

    m.state = m.pparams.build_state_block([0])

    return m


def test_1(model):
    assert model.state[0].temperature.fixed
