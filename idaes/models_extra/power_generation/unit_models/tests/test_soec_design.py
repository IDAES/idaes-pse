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

import pyomo.environ as pyo

from idaes.core import FlowsheetBlock
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
import idaes.core.util.scaling as iscale
from idaes.models_extra.power_generation.properties.natural_gas_PR import get_prop
from idaes.models_extra.power_generation.unit_models.soec_design import (
    SoecDesign,
    EosType,
)
import pytest


def flowsheet(eos=EosType.PR):
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    sweep_comp = {
        "O2": 0.2074,
        "H2O": 0.0099,
        "CO2": 0.0003,
        "N2": 0.7732,
        "Ar": 0.0092,
    }

    feed_comp = {
        "H2": 0.01,
        "H2O": 0.99,
    }

    m.fs.o2_side_prop_params = GenericParameterBlock(
        **get_prop(sweep_comp, {"Vap"}, eos=eos), doc="Air property parameters"
    )
    m.fs.h2_side_prop_params = GenericParameterBlock(
        **get_prop(feed_comp, {"Vap"}, eos=eos), doc="Flue gas property parameters"
    )
    m.fs.soec = SoecDesign(
        oxygen_side_property_package=m.fs.o2_side_prop_params,
        hydrogen_side_property_package=m.fs.h2_side_prop_params,
        reaction_eos=eos,
    )

    m.fs.soec.hydrogen_side_inlet.temperature.fix(1023)
    m.fs.soec.hydrogen_side_inlet.pressure.fix(20e5)
    m.fs.soec.hydrogen_side_inlet.flow_mol.fix(2)
    for c, v in feed_comp.items():
        m.fs.soec.hydrogen_side_inlet.mole_frac_comp[:, c].fix(v)
    m.fs.soec.oxygen_side_inlet.temperature.fix(1023)
    m.fs.soec.oxygen_side_inlet.pressure.fix(20e5)
    m.fs.soec.oxygen_side_inlet.flow_mol.fix(2)
    for c, v in sweep_comp.items():
        m.fs.soec.oxygen_side_inlet.mole_frac_comp[:, c].fix(v)
    m.fs.soec.hydrogen_side_outlet_temperature.fix(1023)
    m.fs.soec.oxygen_side_outlet_temperature.fix(1023)
    m.fs.soec.water_utilization.fix(0.7)
    iscale.calculate_scaling_factors(m)
    m.fs.soec.initialize(optarg={"max_iter": 30})
    return m


@pytest.mark.component
def test_soec_design_ideal():
    m = flowsheet(eos=EosType.IDEAL)
    solver = pyo.SolverFactory("ipopt")
    res = solver.solve(m)
    # Make sure it converged
    assert pyo.check_optimal_termination(res)
    # Based on the inlet composition and the water utilization I know this
    assert pytest.approx(0.703) == pyo.value(
        m.fs.soec.hydrogen_side_outlet.mole_frac_comp[0, "H2"]
    )
    # should be at the thermoneutral volatage and I know about what that is
    assert pytest.approx(1.287, abs=0.02) == pyo.value(m.fs.soec.cell_potential[0])


@pytest.mark.component
def test_soec_design_pr():
    m = flowsheet(eos=EosType.IDEAL)
    solver = pyo.SolverFactory("ipopt")
    res = solver.solve(m)
    # Make sure it converged
    assert pyo.check_optimal_termination(res)
    # Based on the inlet composition and the water utilization I know this
    assert pytest.approx(0.703) == pyo.value(
        m.fs.soec.hydrogen_side_outlet.mole_frac_comp[0, "H2"]
    )
    # should be at the thermoneutral volatage and I know about what that is
    assert pytest.approx(1.287, abs=0.02) == pyo.value(m.fs.soec.cell_potential[0])
