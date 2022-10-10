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

__author__ = "Douglas Allan"

import pytest
import numpy as np

import pyomo.environ as pyo
from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties import (
    GenericParameterBlock,
    GenericStateBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    get_prop,
    EosType,
)
import idaes.models_extra.power_generation.unit_models.soc_submodels.common as common

from idaes.models.properties.modular_properties.base.utility import (
    get_method,
    get_component_object as cobj,
)


def approx(x):
    return pytest.approx(x, rel=1e-6)


@pytest.fixture
def model_ideal():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=[0], time_units=pyo.units.s)
    comps = list(common.h_params.keys())
    comps.remove("Vac")
    comps.remove("O^2-")
    comps.remove("e^-")
    m.fs.paramsIdeal = GenericParameterBlock(
        **get_prop(comps, {"Vap"}, eos=EosType.IDEAL), doc="Air property parameters"
    )
    m.fs.propsIdeal = GenericStateBlock(parameters=m.fs.paramsIdeal, defined_state=True)
    m.fs.propsIdeal.flow_mol.fix(1000)
    # Technically the NIST standard state is 1 bar, but this is what Andrew has at the moment
    m.fs.propsIdeal.pressure.fix(1.01325e5)
    return m


@pytest.mark.component
def test_thermo(model_ideal):
    m = model_ideal
    for comp in m.fs.propsIdeal.component_list:
        obj = cobj(m.fs.propsIdeal, comp)
        # Avoid having to solve the block by getting the ideal gas entropy method directly
        entr_mol_expr = get_method(m.fs.propsIdeal, "entr_mol_ig_comp", comp)
        m.fs.propsIdeal.mole_frac_comp.fix(1e-19)
        m.fs.propsIdeal.mole_frac_comp[comp].fix(1)
        for T in np.linspace(500, 1200, 71):
            m.fs.propsIdeal.temperature.fix(T)
            assert pyo.value(common._comp_int_energy_expr(T, comp)) == approx(
                pyo.value(m.fs.propsIdeal.energy_internal_mol_phase_comp["Vap", comp])
            )
            assert pyo.value(common._comp_enthalpy_expr(T, comp)) == approx(
                pyo.value(m.fs.propsIdeal.enth_mol_phase_comp["Vap", comp])
            )
            assert pyo.value(common._comp_entropy_expr(T, comp)) == approx(
                pyo.value(entr_mol_expr(m.fs.propsIdeal, obj, T * pyo.units.K))
            )
