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

import pyomo.environ as aml
import pyomo.dae as dae
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core.util.model_statistics import degrees_of_freedom
import pytest

__author__ = "Robert Parker"


def make_model(horizon=5, nfe=10, ncp=2):
    """The simplest DAE I can think of. An isothermal CSTR with one
    reaction. No nlp solver should have trouble with this one.
    Other than the flow_out*conc_out term in the material_balances,
    this model is entirely linear.

    """
    m = aml.ConcreteModel()
    m.time = dae.ContinuousSet(bounds=(0, horizon))
    m.components = aml.Set(initialize=["A", "B"])

    m.volume = aml.Param(initialize=1.0)
    m.max_height = aml.Var(initialize=1.0)
    m.max_height.fix()

    stoich_dict = {
        "A": -1,
        "B": 1,
    }
    m.stoich = aml.Param(m.components, initialize=stoich_dict)
    m.k_rxn = aml.Param(initialize=1.0)
    m.rxn_order = aml.Param(initialize=1, domain=aml.Integers)

    # Variables
    m.flow_in = aml.Var(m.time)
    m.flow_out = aml.Var(m.time)

    m.conc_in = aml.Var(m.time, m.components)
    m.conc = aml.Var(m.time, m.components)

    m.rate = aml.Var(m.time, m.components)

    m.dcdt = dae.DerivativeVar(m.conc, wrt=m.time)

    # Equations
    def mb_rule(m, t, j):
        return (
            m.volume * m.dcdt[t, j]
            == m.flow_in[t] * m.conc_in[t, j]
            - m.flow_out[t] * m.conc[t, j]
            + m.volume * m.rate[t, j]
        )

    m.material_balance = aml.Constraint(m.time, m.components, rule=mb_rule)

    def flow_rule(m, t):
        return m.flow_in[t] == m.flow_out[t]

    m.flow_eqn = aml.Constraint(m.time, rule=flow_rule)

    def rate_rule(m, t, j):
        return m.rate[t, j] == m.stoich[j] * m.k_rxn * m.conc[t, "A"] ** m.rxn_order

    m.rate_eqn = aml.Constraint(m.time, m.components, rule=rate_rule)

    disc = aml.TransformationFactory("dae.collocation")
    disc.apply_to(m, wrt=m.time, nfe=nfe, ncp=2, scheme="LAGRANGE-RADAU")

    # Degrees of freedom:
    m.conc_in[:, "A"].fix(5.0)
    m.conc_in[:, "B"].fix(0.0)

    m.conc[0, "A"].fix(0.0)
    m.conc[0, "B"].fix(0.0)

    m.flow_in.fix(1.0)

    return m


def make_small_model():
    return make_model(horizon=1, nfe=2)


def initialize_t0(model):
    time = model.time
    t0 = time.first()
    calculate_variable_from_constraint(
        model.flow_out[t0],
        model.flow_eqn[t0],
    )
    for j in model.components:
        calculate_variable_from_constraint(
            model.rate[t0, j],
            model.rate_eqn[t0, j],
        )
        calculate_variable_from_constraint(
            model.dcdt[t0, j],
            model.material_balance[t0, j],
        )


def copy_values_forward(model):
    time = model.time
    t0 = time.first()
    for t in time:
        for j in model.components:
            model.conc[t, j].set_value(model.conc[t0, j].value)
            model.rate[t, j].set_value(model.rate[t0, j].value)
            model.dcdt[t, j].set_value(model.dcdt[t0, j].value)
        model.flow_out[t].set_value(model.flow_out[t0].value)


@pytest.mark.unit
def test_model():
    m = make_model()
    assert degrees_of_freedom(m) == 0
