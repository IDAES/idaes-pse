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
Tests for translator block.
Authors: Andrew Lee
"""

import pytest

from pyomo.environ import ConcreteModel
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models.translator import Translator

from idaes.generic_models.properties.activity_coeff_models.BTX_activity_coeff_VLE \
    import BTXParameterBlock
from idaes.generic_models.properties.examples.saponification_thermo import \
    SaponificationParameterBlock

from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              number_unused_variables)
from idaes.core.util.testing import (PhysicalParameterTestBlock,
                                     initialization_tester)
from idaes.core.util import get_solver


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = Translator(default={
            "inlet_property_package": m.fs.properties,
            "outlet_property_package": m.fs.properties})

    # Check unit config arguments
    assert len(m.fs.unit.config) == 8

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.inlet_property_package is m.fs.properties
    assert m.fs.unit.config.outlet_property_package is m.fs.properties
    assert m.fs.unit.config.outlet_state_defined
    assert not m.fs.unit.config.has_phase_equilibrium


# -----------------------------------------------------------------------------
class TestTranslate(object):
    @pytest.fixture(scope="class")
    def trans(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties1 = SaponificationParameterBlock()
        m.fs.properties2 = BTXParameterBlock()

        m.fs.unit = Translator(
                default={"inlet_property_package": m.fs.properties1,
                         "outlet_property_package": m.fs.properties2})

        m.fs.unit.inlet.flow_vol.fix(1.0e-03)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, trans):
        assert hasattr(trans.fs.unit, "properties_in")
        assert hasattr(trans.fs.unit, "properties_out")

        assert hasattr(trans.fs.unit, "inlet")
        assert len(trans.fs.unit.inlet.vars) == 4
        assert hasattr(trans.fs.unit.inlet, "flow_vol")
        assert hasattr(trans.fs.unit.inlet, "conc_mol_comp")
        assert hasattr(trans.fs.unit.inlet, "temperature")
        assert hasattr(trans.fs.unit.inlet, "pressure")

        assert hasattr(trans.fs.unit, "outlet")
        assert len(trans.fs.unit.outlet.vars) == 4
        assert hasattr(trans.fs.unit.outlet, "flow_mol")
        assert hasattr(trans.fs.unit.outlet, "mole_frac_comp")
        assert hasattr(trans.fs.unit.outlet, "temperature")
        assert hasattr(trans.fs.unit.outlet, "pressure")

        assert number_variables(trans) == 25
        assert number_total_constraints(trans) == 12
        assert number_unused_variables(trans) == 8

    @pytest.mark.component
    def test_units(self, trans):
        assert_units_consistent(trans)

    @pytest.mark.unit
    def test_dof(self, trans):
        assert degrees_of_freedom(trans) == 5

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, trans):
        initialization_tester(trans, dof=5)

    # No solve, as problem has is missing linking constraints

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, trans):
        trans.fs.unit.report()
