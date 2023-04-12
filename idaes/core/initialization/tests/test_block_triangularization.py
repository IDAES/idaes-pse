#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for Block Triangularization initialization
"""
import pytest
import types

from pyomo.environ import ConcreteModel, Constraint, units, value, Var

from idaes.core import FlowsheetBlock
from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)
from idaes.core.initialization.initializer_base import InitializationStatus
from idaes.models.unit_models.pressure_changer import (
    Turbine,
)

from idaes.models.properties import iapws95

__author__ = "Andrew Lee"


class TestBTSubMethods:
    @pytest.mark.unit
    def test_config(self):
        initializer = BlockTriangularizationInitializer()

        assert hasattr(initializer, "config")

        assert "constraint_tolerance" in initializer.config

        assert "block_solver" in initializer.config
        assert "block_solver_options" in initializer.config
        assert "calculate_variable_options" in initializer.config

    # TODO: Tests for prechecks and initialization_routine stand alone

    @pytest.fixture
    def model(self):
        m = ConcreteModel()

        m.v1 = Var()
        m.v2 = Var()
        m.v3 = Var()
        m.v4 = Var()

        m.c1 = Constraint(expr=m.v1 == m.v2)
        m.c2 = Constraint(expr=2 * m.v2 == m.v3 + m.v4)
        m.c3 = Constraint(expr=m.v3 - m.v4 == 0)

        # Add a dummy method for fixing initialization states
        def fix_initialization_states(blk):
            blk.v1.fix(4)

        m.fix_initialization_states = types.MethodType(fix_initialization_states, m)

        return m

    @pytest.mark.component
    def test_workflow(self, model):
        initializer = BlockTriangularizationInitializer()

        status = initializer.initialize(model)

        assert model.v1.value == 4
        assert model.v2.value == 4
        assert model.v3.value == 4
        assert model.v4.value == 4

        assert not model.v1.fixed

        assert status == InitializationStatus.Ok


@pytest.mark.component
def test_sympy_differentiate_failure():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()

    def perf_callback(b):
        unit_hd = units.J / units.kg
        unit_vflw = units.m**3 / units.s

        @b.Constraint(m.fs.time)
        def pc_isen_eff_eqn(b, t):
            prnt = b.parent_block()
            vflw = prnt.control_volume.properties_in[t].flow_vol
            return prnt.efficiency_isentropic[t] == 0.9 / 3.975 * vflw / unit_vflw

        @b.Constraint(m.fs.time)
        def pc_isen_head_eqn(b, t):
            prnt = b.parent_block()
            vflw = prnt.control_volume.properties_in[t].flow_vol
            return (
                b.head_isentropic[t] / 1000
                == -75530.8 / 3.975 / 1000 * vflw / unit_vflw * unit_hd
            )

    m.fs.unit = Turbine(
        property_package=m.fs.properties,
        support_isentropic_performance_curves=True,
        isentropic_performance_curves={"build_callback": perf_callback},
    )

    # set inputs
    m.fs.unit.inlet.flow_mol[0].fix(1000)  # mol/s
    Tin = 500  # K
    Pin = 1000000  # Pa
    hin = value(iapws95.htpx(Tin * units.K, Pin * units.Pa))
    m.fs.unit.inlet.enth_mol[0].fix(hin)
    m.fs.unit.inlet.pressure[0].fix(Pin)

    initializer = BlockTriangularizationInitializer()

    with pytest.raises(
        TypeError,
        match="A TypeError was encountered in the scc solver. This can occur "
        "in models that involve ExternalFunctions. We suggest you try "
        "setting calculate_variable_options=\u007b'diff_mode': "
        "pyomo.core.expr.calculus.differentiate.Modes.reverse_numeric\u007d. ",
    ):
        initializer.initialize(m)
