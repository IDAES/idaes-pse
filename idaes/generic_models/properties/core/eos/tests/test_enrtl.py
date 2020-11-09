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
Benzene-Toluene phase equilibrium package using ideal liquid and vapor.

Example property package using the Generic Property Package Framework.
This exmample shows how to set up a property package to do benzene-toluene
phase equilibrium in the generic framework using ideal liquid and vapor
assumptions along with methods drawn from the pre-built IDAES property
libraries.
"""
# Import Python libraries
import logging
import pytest

# Import Pyomo units
from pyomo.environ import ConcreteModel, units as pyunits

# Import IDAES cores
from idaes.core import AqueousPhase
from idaes.core.components import *

from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.enrtl import ENRTL

from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)


# Set up logger
_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Configuration dictionary for an ideal Benzene-Toluene system

# Data Sources:
# [1] The Properties of Gases and Liquids (1987)
#     4th edition, Chemical Engineering Series - Robert C. Reid
# [3] Engineering Toolbox, https://www.engineeringtoolbox.com
#     Retrieved 1st December, 2019

configuration = {
    # Specifying components
    "components": {
        'H2O': {"type": Solvent,
                "parameter_data": {
                    "mw": (18E-3, pyunits.kg/pyunits.mol)}},
        'CO2': {"type": Solute,
                "parameter_data": {
                    "mw": (44E-3, pyunits.kg/pyunits.mol)}},
        'KHCO3': {"type": Apparent,
                  "parameter_data": {
                      "mw": (100.1E-3, pyunits.kg/pyunits.mol)}},
        'K+': {"type": Cation,
               "charge": +1,
               "parameter_data": {
                    "mw": (39.1E-3, pyunits.kg/pyunits.mol)}},
        'HCO3-': {"type": Anion,
                  "charge": -1,
                  "parameter_data": {
                       "mw": (61E-3, pyunits.kg/pyunits.mol)}},
        'N2': {"type": Component,
               "parameter_data": {
                    "mw": (28E-3, pyunits.kg/pyunits.mol)}}},

    # Specifying phases
    "phases":  {'Liq': {"type": AqueousPhase,
                        "equation_of_state": ENRTL}},

    # Set base units of measurement
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},

    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {"flow_mol": (0, 100, 1000, pyunits.mol/pyunits.s),
                     "temperature": (273.15, 300, 500, pyunits.K),
                     "pressure": (5e4, 1e5, 1e6, pyunits.Pa)},
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),

    # Defining phase equilibria
    }


@pytest.mark.unit
def test_component_lists():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={'dynamic': False})

    m.fs.props = GenericParameterBlock(default=configuration)

    m.fs.state_1 = m.fs.props.build_state_block(
            [1],
            default={"defined_state": True,
                     "species_basis": "true"})

    m.fs.state_2 = m.fs.props.build_state_block(
            [1],
            default={"defined_state": True,
                     "species_basis": "apparent"})

    assert m.fs.props._electrolyte

    assert m.fs.props.anion_set == ["HCO3-"]
    assert m.fs.props.cation_set == ["K+"]
    assert m.fs.props.solvent_set == ["H2O"]
    assert m.fs.props.solute_set == ["CO2"]
    assert m.fs.props._apparent_set == ["KHCO3"]
    assert m.fs.props._non_aqueous_set == ["N2"]

    assert m.fs.props.true_species_set == [
        "HCO3-", "K+", "H2O", "CO2", "N2"]
    assert m.fs.props.apparent_species_set == [
        "H2O", "CO2", "KHCO3", "N2"]
    assert m.fs.props.component_list == [
        "HCO3-", "K+", "H2O", "CO2", "KHCO3", "N2"]

    assert m.fs.props.true_phase_component_set == [
        ("Liq", "HCO3-"),  ("Liq", "K+"),  ("Liq", "H2O"),
        ("Liq", "CO2"),  ("Liq", "N2")]
    assert m.fs.props.apparent_phase_component_set == [
        ("Liq", "H2O"), ("Liq", "CO2"),  ("Liq", "KHCO3"),  ("Liq", "N2")]

    assert m.fs.state_1[1].component_list is m.fs.props.true_species_set
    assert m.fs.state_2[1].component_list is m.fs.props.apparent_species_set

    assert m.fs.state_1[1].phase_component_set is \
        m.fs.props.true_phase_component_set
    assert m.fs.state_2[1].phase_component_set is \
        m.fs.props.apparent_phase_component_set
