##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Tests for CSTR unit model.
Authors: Andrew Lee, Vibhav Dabadghao
"""

import pytest
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from idaes.unit_models.cstr import CSTR
from idaes.property_models.examples.saponification_thermo import (
    SaponificationParameterBlock)
from idaes.property_models.examples.saponification_reactions import (
    SaponificationReactionParameterBlock)
from idaes.ui.report import (degrees_of_freedom,
                             fixed_variables,
                             stale_variables,
                             unfixed_variables,
                             active_equalities,
                             inactive_equalities)
from idaes.core.util.testing import (count_constraints,
                                     count_variables,
                                     get_default_solver)


# Get default solver from test utils
solver = get_default_solver()


# -----------------------------------------------------------------------------
# create basic model for testing
@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = SaponificationParameterBlock()
    m.fs.reactions = SaponificationReactionParameterBlock(default={
                            "property_package": m.fs.properties})

    m.fs.cstr = CSTR(default={"property_package": m.fs.properties,
                              "reaction_package": m.fs.reactions,
                              "has_equilibrium_reactions": False,
                              "has_heat_transfer": True,
                              "has_pressure_change": False})

    return m

# -----------------------------------------------------------------------------
@pytest.mark.model_build
def test_build(model):
    assert hasattr(model.fs.cstr, "inlet")
    assert len(model.fs.cstr.inlet.vars) == 4
    assert hasattr(model.fs.cstr.inlet, "flow_vol")
    assert hasattr(model.fs.cstr.inlet, "conc_mol_comp")
    assert hasattr(model.fs.cstr.inlet, "temperature")
    assert hasattr(model.fs.cstr.inlet, "pressure")

    assert hasattr(model.fs.cstr, "outlet")
    assert len(model.fs.cstr.outlet.vars) == 4
    assert hasattr(model.fs.cstr.outlet, "flow_vol")
    assert hasattr(model.fs.cstr.outlet, "conc_mol_comp")
    assert hasattr(model.fs.cstr.outlet, "temperature")
    assert hasattr(model.fs.cstr.outlet, "pressure")

    assert hasattr(model.fs.cstr, "cstr_performance_eqn")
    assert hasattr(model.fs.cstr.control_volume, "heat")
    assert hasattr(model.fs.cstr, "heat_duty")


@pytest.mark.degrees_of_freedom
def test_dof(model):
    model.fs.cstr.inlet.flow_vol.fix(1.0e-03)
    model.fs.cstr.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
    model.fs.cstr.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
    model.fs.cstr.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
    model.fs.cstr.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
    model.fs.cstr.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

    model.fs.cstr.inlet.temperature.fix(303.15)
    model.fs.cstr.inlet.pressure.fix(101325.0)

    model.fs.cstr.volume.fix(1.5e-03)
    model.fs.cstr.heat_duty.fix(0)

    assert degrees_of_freedom(model) == 0

    assert count_variables(model) == 42
    assert count_constraints(model) == 16


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.initialization
def test_initialize(model):
    assert degrees_of_freedom(model) == 0

    f_vars_1 = fixed_variables(model)
    u_vars_1 = unfixed_variables(model)
    a_cons_1 = active_equalities(model)
    i_cons_1 = inactive_equalities(model)

    model.fs.cstr.initialize(outlvl=5,
                             optarg={'tol': 1e-6})

    f_vars_2 = fixed_variables(model)
    u_vars_2 = unfixed_variables(model)
    a_cons_2 = active_equalities(model)
    i_cons_2 = inactive_equalities(model)

    for v in f_vars_1:
        assert v in f_vars_2
    for v in u_vars_1:
        assert v in u_vars_2
    for c in a_cons_1:
        assert c in a_cons_2
    for c in i_cons_1:
        assert c in i_cons_2

    assert len(list(stale_variables(model))) == 0

    assert (pytest.approx(101325.0, abs=1e-2) ==
            model.fs.cstr.outlet.pressure[0].value)
    assert (pytest.approx(303.15, abs=1e-2) ==
            model.fs.cstr.outlet.temperature[0].value)
    assert (pytest.approx(20.80, abs=1e-2) ==
            model.fs.cstr.outlet.conc_mol_comp[0, "EthylAcetate"].value)
