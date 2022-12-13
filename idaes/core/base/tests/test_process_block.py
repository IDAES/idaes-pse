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
"""
Tests ProcessBlock and ProcessBlockData.

Author: John Eslick
"""
import pytest
from pyomo.environ import ConcreteModel, Var, value
from pyomo.common.config import ConfigValue
from idaes.core import ProcessBlockData, declare_process_block_class


@declare_process_block_class("MyBlock")
class MyBlockData(ProcessBlockData):
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare("xinit", ConfigValue(default=1001, domain=float))
    CONFIG.declare("yinit", ConfigValue(default=1002, domain=float))

    def build(self):
        super(MyBlockData, self).build()
        self.x = Var(initialize=self.config.xinit)
        self.y = Var(initialize=self.config.yinit)


class TestProcessBlock(object):
    @pytest.mark.unit
    def test_scalar_noargs(self):
        m = ConcreteModel()
        m.b = MyBlock()
        assert isinstance(m.b.x, Var)
        assert isinstance(m.b.y, Var)
        assert value(m.b.x) == 1001
        assert value(m.b.y) == 1002

    @pytest.mark.unit
    def test_vec_noargs(self):
        m = ConcreteModel()
        m.b = MyBlock([1, 2, 3])
        assert isinstance(m.b[1].x, Var)
        assert isinstance(m.b[1].y, Var)
        assert isinstance(m.b[2].x, Var)
        assert isinstance(m.b[2].y, Var)
        assert isinstance(m.b[3].x, Var)
        assert isinstance(m.b[3].y, Var)
        assert value(m.b[1].x) == 1001
        assert value(m.b[1].y) == 1002
        assert value(m.b[2].x) == 1001
        assert value(m.b[2].y) == 1002
        assert value(m.b[3].x) == 1001
        assert value(m.b[3].y) == 1002

    @pytest.mark.unit
    def test_scalar_args1(self):
        m = ConcreteModel()
        m.b = MyBlock(xinit=1, yinit=2)
        assert isinstance(m.b.x, Var)
        assert isinstance(m.b.y, Var)
        assert value(m.b.x) == 1
        assert value(m.b.y) == 2

    @pytest.mark.unit
    def test_scalar_args2(self):
        m = ConcreteModel()
        m.b = MyBlock(initialize={None: {"xinit": 1, "yinit": 2}})
        assert isinstance(m.b.x, Var)
        assert isinstance(m.b.y, Var)
        assert value(m.b.x) == 1
        assert value(m.b.y) == 2

    @pytest.mark.unit
    def test_vec_args(self):
        m = ConcreteModel()
        m.b = MyBlock(
            [1, 2, 3], xinit=1, yinit=2, initialize={2: {"xinit": 2001, "yinit": 2002}}
        )
        assert isinstance(m.b[1].x, Var)
        assert isinstance(m.b[1].y, Var)
        assert isinstance(m.b[2].x, Var)
        assert isinstance(m.b[2].y, Var)
        assert isinstance(m.b[3].x, Var)
        assert isinstance(m.b[3].y, Var)
        assert m.b[1].index() == 1
        assert value(m.b[1].x) == 1
        assert value(m.b[1].y) == 2
        assert value(m.b[2].x) == 2001
        assert value(m.b[2].y) == 2002
        assert value(m.b[3].x) == 1
        assert value(m.b[3].y) == 2

    @pytest.mark.unit
    def test_user_map(self):
        m = ConcreteModel()

        def new_imap(i):
            if i > 0 and i < 4:
                return 1
            else:
                return i

        m.b = MyBlock(
            [0, 1, 2, 3, 4],
            idx_map=new_imap,
            initialize={
                0: {"xinit": 2001, "yinit": 2002},
                1: {"xinit": 5001, "yinit": 5002},
                2: {"xinit": 6001, "yinit": 6002},
                4: {"xinit": 7001, "yinit": 7002},
            },
        )
        assert value(m.b[0].x) == 2001
        assert value(m.b[0].y) == 2002
        assert value(m.b[1].x) == 5001
        assert value(m.b[1].y) == 5002
        assert value(m.b[2].x) == 5001  # although this index (2) is in initialize
        assert value(m.b[2].y) == 5002  # the idx_map function maps it to 1
        assert value(m.b[3].x) == 5001
        assert value(m.b[3].y) == 5002
        assert value(m.b[4].x) == 7001
        assert value(m.b[4].y) == 7002
