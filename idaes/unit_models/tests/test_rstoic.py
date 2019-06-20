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
Tests for IDAES Stoichiometric reactor.

Author: Chinedu Okoli, Andrew Lee
"""
import pytest

from pyomo.environ import ConcreteModel, SolverFactory

from idaes.core import FlowsheetBlock
from idaes.unit_models.stoichiometric_reactor import StoichiometricReactor
from idaes.property_models.examples.saponification_thermo import (
    SaponificationParameterBlock)
from idaes.property_models.examples.saponification_reactions import (
    SaponificationReactionParameterBlock)
from idaes.core.util.model_statistics import degrees_of_freedom


# -----------------------------------------------------------------------------
# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None


# -----------------------------------------------------------------------------
def test_build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = SaponificationParameterBlock()
    m.fs.reactions = SaponificationReactionParameterBlock(default={
                            "property_package": m.fs.properties})

    m.fs.rstoic = StoichiometricReactor(default={
                            "property_package": m.fs.properties,
                            "reaction_package": m.fs.reactions,
                            "has_heat_transfer": True,
                            "has_pressure_change": True})

    assert hasattr(m.fs.rstoic, "inlet")
    assert len(m.fs.rstoic.inlet.vars) == 4
    assert hasattr(m.fs.rstoic.inlet, "flow_vol")
    assert hasattr(m.fs.rstoic.inlet, "conc_mol_comp")
    assert hasattr(m.fs.rstoic.inlet, "temperature")
    assert hasattr(m.fs.rstoic.inlet, "pressure")

    assert hasattr(m.fs.rstoic, "outlet")
    assert len(m.fs.rstoic.outlet.vars) == 4
    assert hasattr(m.fs.rstoic.outlet, "flow_vol")
    assert hasattr(m.fs.rstoic.outlet, "conc_mol_comp")
    assert hasattr(m.fs.rstoic.outlet, "temperature")
    assert hasattr(m.fs.rstoic.outlet, "pressure")

    assert hasattr(m.fs.rstoic.control_volume, "heat")
    assert hasattr(m.fs.rstoic, "heat_duty")
    assert hasattr(m.fs.rstoic, "deltaP")


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = SaponificationParameterBlock()
    m.fs.reactions = SaponificationReactionParameterBlock(default={
                            "property_package": m.fs.properties})

    m.fs.rstoic = StoichiometricReactor(default={
                            "property_package": m.fs.properties,
                            "reaction_package": m.fs.reactions,
                            "has_heat_transfer": False,
                            "has_pressure_change": False})

    m.fs.rstoic.inlet.flow_vol.fix(1.0)
    m.fs.rstoic.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
    m.fs.rstoic.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
    m.fs.rstoic.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
    m.fs.rstoic.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
    m.fs.rstoic.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

    m.fs.rstoic.inlet.temperature.fix(303.15)
    m.fs.rstoic.inlet.pressure.fix(101325.0)

    m.fs.rstoic.rate_reaction_extent[:, 'R1'].fix(
            0.9*m.fs.rstoic.inlet.conc_mol_comp[0, "NaOH"].value)

    assert degrees_of_freedom(m) == 0

    m.fs.rstoic.initialize(outlvl=5,
                           optarg={'tol': 1e-6})

    assert (pytest.approx(101325.0, abs=1e-2) ==
            m.fs.rstoic.outlet.pressure[0].value)
    assert (pytest.approx(303.15, abs=1e-2) ==
            m.fs.rstoic.outlet.temperature[0].value)
    assert (pytest.approx(90, abs=1e-2) ==
            m.fs.rstoic.outlet.conc_mol_comp[0, "Ethanol"].value)
