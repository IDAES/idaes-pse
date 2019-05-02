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

import pandas as pd
import pyomo.environ as pyo
from pyomo.common.config import ConfigBlock
from idaes.core import StateBlock, StateBlockData, declare_process_block_class
from idaes.core.util.tables import stream_table, state_table
###
# Create a dummy StateBlock class
##

@declare_process_block_class("TestStateBlock", block_class=StateBlock)
class StateTestBlockData(StateBlockData):
    def build(self):
        self.x = pyo.Var(initialize=0)
        self.pressure = pyo.Var()
        self.enth_mol = pyo.Var()
        self.temperature = pyo.Expression(expr=self.enth_mol*2.0)
        self.div0 = pyo.Expression(expr=4.0/self.x)
        self.flow_mol = pyo.Var()


def make_model():
    m = pyo.ConcreteModel()
    m.state_a = TestStateBlock()
    m.state_b = TestStateBlock([1,2,3])

    m.state_a.pressure = 11000
    m.state_a.enth_mol = 1100
    m.state_a.flow_mol = 110

    m.state_b[1].pressure = 10000
    m.state_b[1].enth_mol = 1000
    m.state_b[1].flow_mol = 100

    m.state_b[2].pressure = 20000
    m.state_b[2].enth_mol = 2000
    m.state_b[2].flow_mol = 200

    m.state_b[3].pressure = 30000
    m.state_b[3].enth_mol = 3000
    m.state_b[3].flow_mol = 300

    return m


def test_stream_table():
    m = make_model()

    sd = {"a":m.state_a, "b1": m.state_b[1]}
    st = stream_table(sd,
        attributes=["pressure", "temperature", "div0", "flow_mol", "not_there"],
        heading=["P", "T", "ERR", "F", "Miss"])

    assert(st.loc["a"]["P"]==11000)
    assert(abs(st.loc["a"]["T"] - 1100*2) < 0.001)
    assert(pd.isna(st.loc["a"]["ERR"]))
    assert(pd.isna(st.loc["a"]["Miss"]))

    assert(st.loc["b1"]["P"]==10000)
    assert(abs(st.loc["b1"]["T"] - 1000*2) < 0.001)
    assert(pd.isna(st.loc["b1"]["ERR"]))
    assert(pd.isna(st.loc["b1"]["Miss"]))

def test_stream_table():
    m = make_model()

    st = state_table(m,
        attributes=["pressure", "temperature", "div0", "flow_mol", "not_there"],
        heading=["P", "T", "ERR", "F", "Miss"])
    print(st)
    assert(st.loc["state_a"]["P"]==11000)
    assert(abs(st.loc["state_a"]["T"] - 1100*2) < 0.001)
    assert(pd.isna(st.loc["state_a"]["ERR"]))
    assert(pd.isna(st.loc["state_a"]["Miss"]))

    assert(st.loc["state_b[1]"]["P"]==10000)
    assert(abs(st.loc["state_b[1]"]["T"] - 1000*2) < 0.001)
    assert(pd.isna(st.loc["state_b[1]"]["ERR"]))
    assert(pd.isna(st.loc["state_b[1]"]["Miss"]))
