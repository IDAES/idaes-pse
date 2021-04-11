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
Test categorize_dae_variables_and_constraints function.
"""

import pytest
import pyomo.environ as pyo
import pyomo.dae as dae
from pyomo.common.collections import ComponentSet
from pyomo.dae.flatten import flatten_dae_components

from idaes.apps.caprese.categorize import (
        categorize_dae_variables_and_constraints,
        )
from idaes.apps.caprese.common.config import VariableCategory as VC
from idaes.apps.caprese.common.config import ConstraintCategory as CC

__author__ = "Robert Parker"

@pytest.mark.unit
def test_categorize_deriv():
    """ The simplest test. Identify a differential and a derivative var.
    """
    m = pyo.ConcreteModel()
    m.time = dae.ContinuousSet(initialize=[0, 1])
    m.v = pyo.Var(m.time, initialize=0)
    m.dv = dae.DerivativeVar(m.v, wrt=m.time)
    m.diff_eqn = pyo.Constraint(
            m.time,
            rule={t: m.dv[t] == -m.v[t]**2 for t in m.time}
            )
    with pytest.raises(TypeError):
        # If we find a derivative var, we will try to access the disc eq.
        var_partition, con_partition = categorize_dae_variables_and_constraints(
                m,
                [m.v, m.dv],
                [m.diff_eqn],
                m.time,
                )
    disc = pyo.TransformationFactory('dae.finite_difference')
    disc.apply_to(m, wrt=m.time, nfe=1, scheme='BACKWARD')
    var_partition, con_partition = categorize_dae_variables_and_constraints(
            m,
            [m.v, m.dv],
            [m.diff_eqn, m.dv_disc_eq],
            m.time,
            )
    assert len(var_partition[VC.DIFFERENTIAL]) == 1
    assert var_partition[VC.DIFFERENTIAL][0] is m.v
    assert len(var_partition[VC.DERIVATIVE]) == 1
    assert var_partition[VC.DERIVATIVE][0] is m.dv
    assert len(con_partition[CC.DIFFERENTIAL]) == 1
    assert con_partition[CC.DIFFERENTIAL][0] is m.diff_eqn
    assert len(con_partition[CC.DISCRETIZATION]) == 1
    assert con_partition[CC.DISCRETIZATION][0] is m.dv_disc_eq

    for categ in VC:
        if (categ is not VC.DIFFERENTIAL and categ is not VC.DERIVATIVE
                and categ in var_partition):
            assert len(var_partition[categ]) == 0
    for categ in CC:
        if (categ is not CC.DIFFERENTIAL and categ is not CC.DISCRETIZATION
                and categ in con_partition):
            assert len(con_partition[categ]) == 0


@pytest.mark.unit
def test_categorize_deriv_fixed():
    """ If one of the derivative or diff var are fixed, the other
    should be categorized as algebraic.
    """
    m = pyo.ConcreteModel()
    m.time = dae.ContinuousSet(initialize=[0, 1])
    m.v = pyo.Var(m.time, initialize=0)
    m.dv = dae.DerivativeVar(m.v, wrt=m.time)
    m.diff_eqn = pyo.Constraint(
            m.time,
            rule={t: m.dv[t] == -m.v[t]**2 for t in m.time}
            )
    disc = pyo.TransformationFactory('dae.finite_difference')
    disc.apply_to(m, wrt=m.time, nfe=1, scheme='BACKWARD')
    
    #
    # Fix differential variable, e.g. it is an input
    #
    m.v.fix()
    m.diff_eqn.deactivate()

    var_partition, con_partition = categorize_dae_variables_and_constraints(
            m,
            [m.v, m.dv],
            [m.diff_eqn, m.dv_disc_eq],
            m.time,
            )
    # Expected categories have expected variables
    assert len(var_partition[VC.ALGEBRAIC]) == 1
    assert var_partition[VC.ALGEBRAIC][0] is m.dv
    assert len(con_partition[CC.ALGEBRAIC]) == 1
    assert con_partition[CC.ALGEBRAIC][0] is m.dv_disc_eq

    # Unexpected categories are empty
    for categ in VC:
        if (categ is not VC.ALGEBRAIC and categ is not VC.UNUSED
                and categ in var_partition):
            assert len(var_partition[categ]) == 0
    for categ in CC:
        if (categ is not CC.ALGEBRAIC and categ is not CC.UNUSED
                and categ in con_partition):
            assert len(con_partition[categ]) == 0

    #
    # We can accomplish something similar by making m.v an input
    #
    m.v.unfix()
    var_partition, con_partition = categorize_dae_variables_and_constraints(
            m,
            [m.v, m.dv],
            [m.diff_eqn, m.dv_disc_eq],
            m.time,
            input_vars=[m.v],
            )
    # Expected categories have expected variables
    assert len(var_partition[VC.ALGEBRAIC]) == 1
    assert var_partition[VC.ALGEBRAIC][0] is m.dv
    assert len(var_partition[VC.INPUT]) == 1
    assert var_partition[VC.INPUT][0] is m.v
    assert len(con_partition[CC.ALGEBRAIC]) == 1
    assert con_partition[CC.ALGEBRAIC][0] is m.dv_disc_eq

    #
    # Fix derivative var, e.g. pseudo-steady state
    #
    for var in m.dv.values():
        var.fix(0.0)
    m.diff_eqn.activate()
    m.dv_disc_eq.deactivate()
    var_partition, con_partition = categorize_dae_variables_and_constraints(
            m,
            [m.v, m.dv],
            [m.diff_eqn, m.dv_disc_eq],
            m.time,
            )
    # Expected categories have expected variables
    assert len(var_partition[VC.ALGEBRAIC]) == 1
    assert var_partition[VC.ALGEBRAIC][0] is m.v
    assert len(con_partition[CC.ALGEBRAIC]) == 1
    assert con_partition[CC.ALGEBRAIC][0] is m.diff_eqn


@pytest.mark.unit
def test_categorize_simple_model():
    """ Categorize variables and equations in the "simple model" used
    for the base class unit tests.
    """
    pass
