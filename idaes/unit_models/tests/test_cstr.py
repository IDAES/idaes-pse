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
from pyomo.environ import ConcreteModel, SolverFactory
from idaes.core import FlowsheetBlock
from idaes.unit_models.cstr import CSTR
from idaes.property_models.examples.saponification_thermo import (
    SaponificationParameterBlock)
from idaes.property_models.examples.saponification_reactions import (
    SaponificationReactionParameterBlock)
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_fixed_variables,
                                              number_activated_constraints,
                                              number_unused_variables)


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
class TestSaponification(object):
    @pytest.fixture(scope="class")
    def cstr_sap(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = SaponificationParameterBlock()
        m.fs.reactions = SaponificationReactionParameterBlock(default={
                                "property_package": m.fs.properties})

        m.fs.cstr = CSTR(default={"property_package": m.fs.properties,
                                  "reaction_package": m.fs.reactions,
                                  "has_equilibrium_reactions": False,
                                  "has_heat_transfer": True,
                                  "has_heat_of_reaction": True,
                                  "has_pressure_change": True})

        return m

    @pytest.mark.build
    def test_build(self, cstr_sap):

        assert hasattr(cstr_sap.fs.cstr, "inlet")
        assert len(cstr_sap.fs.cstr.inlet.vars) == 4
        assert hasattr(cstr_sap.fs.cstr.inlet, "flow_vol")
        assert hasattr(cstr_sap.fs.cstr.inlet, "conc_mol_comp")
        assert hasattr(cstr_sap.fs.cstr.inlet, "temperature")
        assert hasattr(cstr_sap.fs.cstr.inlet, "pressure")

        assert hasattr(cstr_sap.fs.cstr, "outlet")
        assert len(cstr_sap.fs.cstr.outlet.vars) == 4
        assert hasattr(cstr_sap.fs.cstr.outlet, "flow_vol")
        assert hasattr(cstr_sap.fs.cstr.outlet, "conc_mol_comp")
        assert hasattr(cstr_sap.fs.cstr.outlet, "temperature")
        assert hasattr(cstr_sap.fs.cstr.outlet, "pressure")

        assert hasattr(cstr_sap.fs.cstr, "cstr_performance_eqn")
        assert hasattr(cstr_sap.fs.cstr, "volume")
        assert hasattr(cstr_sap.fs.cstr, "heat_duty")
        assert hasattr(cstr_sap.fs.cstr, "deltaP")

        assert number_variables(cstr_sap) == 27
        assert number_total_constraints(cstr_sap) == 16
        assert number_unused_variables(cstr_sap) == 0

    def test_dof(self, cstr_sap):
        cstr_sap.fs.cstr.inlet.flow_vol.fix(1.0e-03)
        cstr_sap.fs.cstr.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        cstr_sap.fs.cstr.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        cstr_sap.fs.cstr.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        cstr_sap.fs.cstr.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        cstr_sap.fs.cstr.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        cstr_sap.fs.cstr.inlet.temperature.fix(303.15)
        cstr_sap.fs.cstr.inlet.pressure.fix(101325.0)

        cstr_sap.fs.cstr.volume.fix(1.5e-03)
        cstr_sap.fs.cstr.heat_duty.fix(0)
        cstr_sap.fs.cstr.deltaP.fix(0)

        assert degrees_of_freedom(cstr_sap) == 0

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_initialize(self, cstr_sap):
        cstr_sap.fs.cstr.initialize(optarg={'tol': 1e-6})

        assert degrees_of_freedom(cstr_sap) == 0
        assert number_activated_constraints(cstr_sap) == 16
        assert number_fixed_variables(cstr_sap) == 11

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solution(self, cstr_sap):
        assert (pytest.approx(101325.0, abs=1e-2) ==
                cstr_sap.fs.cstr.outlet.pressure[0].value)
        assert (pytest.approx(304.09, abs=1e-2) ==
                cstr_sap.fs.cstr.outlet.temperature[0].value)
        assert (pytest.approx(20.32, abs=1e-2) ==
                cstr_sap.fs.cstr.outlet.conc_mol_comp[0, "EthylAcetate"].value)

    @pytest.mark.ui
    def test_report(self, cstr_sap):
        cstr_sap.fs.cstr.report()
