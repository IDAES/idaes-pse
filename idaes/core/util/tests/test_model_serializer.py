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
Test for functions to save/load Pyomo model state to a dict or json
"""

import unittest
import os

from pyomo.environ import *
from idaes.core.util import to_json, from_json, StoreSpec
from idaes.util.system import mkdtemp
import shutil
import pytest

__author__ = "John Eslick"


class TestModelSerialize(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dirname = mkdtemp()
        cls.fname = os.path.join(cls.dirname, "crAzYStuff1010202030.json")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dirname)

    def tearDown(self):
        try:
            os.remove(self.fname)
        except:
            pass

    def setup_model01(self):
        model = ConcreteModel()
        model.b = Block([1,2,3])
        a = model.b[1].a = Var(bounds=(-100, 100), initialize=2)
        b = model.b[1].b = Var(bounds=(-100, 100), initialize=20)
        model.b[1].c = Constraint(expr=b==10*a)
        a.fix(2)
        return model

    def setup_model01b(self):
        model = ConcreteModel()
        model.b = Block(["1","2","3"])
        a = model.b["1"].a = Var(bounds=(-100, 100), initialize=2)
        b = model.b["1"].b = Var(bounds=(-100, 100), initialize=20)
        model.b["1"].c = Constraint(expr=b==10*a)
        a.fix(2)
        return model

    def setup_model02(self):
        model = ConcreteModel()
        a = model.a = Param(default=1, mutable=True)
        b = model.b = Param(default=2, mutable=True)
        c = model.c = Param(initialize=4)
        x = model.x = Var([1,2], initialize={1:1.5, 2:2.5}, bounds=(-10,10))
        model.f = Objective(expr=(x[1] - a)**2 + (x[2] - b)**2)
        model.g = Constraint(expr=x[1] + x[2] - c >= 0)
        model.dual = Suffix(direction=Suffix.IMPORT)
        model.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
        model.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
        return model

    def setup_model02b(self):
        model = ConcreteModel()
        model.p = Param(["a", "b"], default={"a":1, "b":2}, mutable=True)
        a = model.p["a"]
        b = model.p["b"]
        c = model.c = Param(initialize=4)
        x = model.x = Var(["1","2"], initialize={"1":1.5, "2":2.5}, bounds=(-10,10))
        model.f = Objective(expr=(x["1"] - a)**2 + (x["2"] - b)**2)
        model.g = Constraint(expr=x["1"] + x["2"] - c >= 0)
        model.dual = Suffix(direction=Suffix.IMPORT)
        model.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
        model.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
        return model


    @pytest.mark.unit
    def test01(self):
        """
        Simple test of load save json
        """
        model = self.setup_model01()
        a = model.b[1].a
        b = model.b[1].b
        to_json(model, fname=self.fname, human_read=True)
        # change variable values
        a.value = 0.11
        b.value = 0.11
        a.unfix()
        model.b[1].deactivate()
        b.setlb(2)
        b.setub(4)
        # reload values
        from_json(model, fname=self.fname)
        #make sure they are right
        assert(a.fixed)
        assert(model.b[1].active)
        assert(abs(value(b) - 20) < 1e-4)
        assert(abs(value(a) - 2) < 1e-4)
        assert(abs(b.lb - -100) < 1e-4)
        assert(abs(b.ub - 100) < 1e-4)

    @pytest.mark.unit
    def test01b(self):
        """
        Simple test of load save json
        """
        model = self.setup_model01b()
        a = model.b["1"].a
        b = model.b["1"].b
        to_json(model, fname=self.fname, human_read=True)
        # change variable values
        a.value = 0.11
        b.value = 0.11
        a.unfix()
        model.b["1"].deactivate()
        b.setlb(2)
        b.setub(4)
        # reload values
        from_json(model, fname=self.fname)
        #make sure they are right
        assert(a.fixed)
        assert(model.b["1"].active)
        assert(abs(value(b) - 20) < 1e-4)
        assert(abs(value(a) - 2) < 1e-4)
        assert(abs(b.lb - -100) < 1e-4)
        assert(abs(b.ub - 100) < 1e-4)

    @pytest.mark.unit
    def test02(self):
        """Test with suffixes"""
        model = self.setup_model02()
        x = model.x
        model.dual[model.g] = 1
        model.ipopt_zL_out[x[1]] = 0
        model.ipopt_zL_out[x[2]] = 0
        model.ipopt_zU_out[x[1]] = 0
        model.ipopt_zU_out[x[2]] = 0
        to_json(model, fname=self.fname, human_read=True)
        model.x[1].value = 10
        model.x[2].value = 10
        model.dual[model.g] = 10
        model.ipopt_zL_out[x[1]] = 10
        model.ipopt_zL_out[x[2]] = 10
        model.ipopt_zU_out[x[1]] = 10
        model.ipopt_zU_out[x[2]] = 10
        from_json(model, fname=self.fname)
        assert(abs(value(x[1]) - 1.5) < 1e-5)
        assert(abs(value(x[2]) - 2.5) < 1e-5)
        assert(abs(model.dual[model.g] - 1) < 1e-5)
        assert(abs(model.ipopt_zL_out[x[1]]) < 1e-5)
        assert(abs(model.ipopt_zL_out[x[2]]) < 1e-5)
        assert(abs(model.ipopt_zU_out[x[1]]) < 1e-5)
        assert(abs(model.ipopt_zU_out[x[2]]) < 1e-5)

    @pytest.mark.unit
    def test02b(self):
        """Test with suffixes"""
        model = self.setup_model02b()
        x = model.x
        model.dual[model.g] = 1
        model.ipopt_zL_out[x["1"]] = 0
        model.ipopt_zL_out[x["2"]] = 0
        model.ipopt_zU_out[x["1"]] = 0
        model.ipopt_zU_out[x["2"]] = 0
        to_json(model, fname=self.fname, human_read=True)
        model.x["1"].value = 10
        model.x["2"].value = 10
        model.dual[model.g] = 10
        model.ipopt_zL_out[x["1"]] = 10
        model.ipopt_zL_out[x["2"]] = 10
        model.ipopt_zU_out[x["1"]] = 10
        model.ipopt_zU_out[x["2"]] = 10
        model.p["a"] = 10
        model.p["b"] = 10
        from_json(model, fname=self.fname)
        assert value(model.x["1"]) == pytest.approx(1.5)
        assert value(model.x["2"]) == pytest.approx(2.5)
        assert value(model.p["a"]) == pytest.approx(1)
        assert value(model.p["b"]) == pytest.approx(2)
        assert model.dual[model.g] == pytest.approx(1)
        assert model.ipopt_zL_out[x["1"]] == pytest.approx(0)
        assert model.ipopt_zL_out[x["2"]] == pytest.approx(0)
        assert model.ipopt_zU_out[x["1"]] == pytest.approx(0)
        assert model.ipopt_zU_out[x["2"]] == pytest.approx(0)

    @pytest.mark.unit
    def test03(self):
        """
        This tests a StoreSpec object meant for initialization.  It reloads
        the saved state of variable values and whether they are fixed or
        unfixed but it only loads values for variables that were originally
        fixed. It's use would be in doing initialization that changes which
        variables are fixed and does whatever, but when the original state is
        reloaded only the original values for fixed variables are reloaded.  You
        basically end up with the same problem, just a different initial guess.
        """
        model = self.setup_model02()
        x = model.x
        x[1].fix(1)
        wts = StoreSpec.value_isfixed(only_fixed=True)
        to_json(model, fname=self.fname, human_read=True, wts=wts)
        x[1].unfix()
        x[1].value = 2
        x[2].value = 10
        from_json(model, fname=self.fname, wts=wts)
        assert(x[1].fixed)
        assert(abs(value(x[1]) - 1) < 1e-5)
        assert(abs(value(x[2]) - 10) < 1e-5)

    @pytest.mark.unit
    def test04(self):
        """
        Like test03, but this StoreSpec also saves/loads active/deactivated
        component attribute and parameter values.
        """
        model = self.setup_model02()
        x = model.x
        x[1].fix(1)
        wts = StoreSpec.value_isfixed_isactive(only_fixed=True)
        to_json(model, fname=self.fname, human_read=True, wts=wts)
        x[1].unfix()
        x[1].value = 2
        x[2].value = 10
        model.g.deactivate()
        from_json(model, fname=self.fname, wts=wts)
        assert(x[1].fixed)
        assert(abs(value(x[1]) - 1) < 1e-5)
        assert(abs(value(x[2]) - 10) < 1e-5)
        assert(model.g.active)

    @pytest.mark.unit
    def test05(self):
        """Try just saving values"""
        model = self.setup_model02()
        model.x[1].value = 1
        wts = StoreSpec.value()
        to_json(model, fname=self.fname, human_read=True, wts=wts)
        model.g.deactivate()
        model.x[1].setlb(-4)
        model.x[1].value = 3
        model.x[2].value = 6
        from_json(model, fname=self.fname, wts=wts)
        assert(abs(value(model.x[1]) - 1) < 1e-5)
        assert(abs(model.x[1].lb + 4) < 1e-5)
        assert(abs(value(model.x[2]) - 2.5) < 1e-5)
        assert(not model.g.active)

    @pytest.mark.unit
    def test06(self):
        """Try just saving bounds"""
        model = self.setup_model02()
        model.x[1].value = 1
        wts = StoreSpec.bound()
        to_json(model, fname=self.fname, human_read=True, wts=wts)
        model.g.deactivate()
        model.x[1].setlb(-4)
        model.x[1].value = 3
        model.x[2].value = 6
        from_json(model, fname=self.fname, wts=wts)
        assert(abs(value(model.x[1]) - 3) < 1e-5)
        assert(abs(model.x[1].lb + 10) < 1e-5)
        assert(abs(value(model.x[2]) - 6) < 1e-5)
        assert(not model.g.active)

    @pytest.mark.unit
    def test07(self):
        """Try just saving just if fixed"""
        model = self.setup_model02()
        model.x[1].fix(1)
        wts = StoreSpec.isfixed()
        to_json(model, fname=self.fname, human_read=True, wts=wts)
        model.g.deactivate()
        model.x[1].setlb(-4)
        model.x[1].unfix()
        model.x[2].fix(6)
        from_json(model, fname=self.fname, wts=wts)
        assert(abs(value(model.x[1]) - 1) < 1e-5)
        assert(abs(model.x[1].lb + 4) < 1e-5)
        assert(abs(value(model.x[2]) - 6) < 1e-5)
        assert(model.x[1].fixed)
        assert(not model.x[2].fixed)
        assert(not model.g.active)

    @pytest.mark.unit
    def test08(self):
        """Try just saving suffixes"""
        model = self.setup_model02()
        model.x[1].fix(1)

        model.dual[model.g] = 1
        model.ipopt_zL_out[model.x[1]] = 1
        model.ipopt_zL_out[model.x[2]] = 1
        model.ipopt_zU_out[model.x[1]] = 1
        model.ipopt_zU_out[model.x[2]] = 1

        to_json(model, fname=self.fname, wts=StoreSpec.suffix())

        model.dual[model.g] = 10
        model.ipopt_zL_out[model.x[1]] = 10
        model.ipopt_zL_out[model.x[2]] = 10
        model.ipopt_zU_out[model.x[1]] = 10
        model.ipopt_zU_out[model.x[2]] = 10
        model.g.deactivate()
        model.x[1].setlb(-4)
        model.x[1].unfix()
        model.x[2].fix(6)

        from_json(model, fname=self.fname, wts=StoreSpec.suffix())
        assert(abs(value(model.x[1]) - 1) < 1e-5)
        assert(abs(value(model.x[2]) - 6) < 1e-5)
        assert(not model.x[1].fixed)
        assert(model.x[2].fixed)
        assert(not model.g.active)
        assert(abs(model.dual[model.g] - 1) < 1e-5)
        assert(abs(model.ipopt_zL_out[model.x[1]] - 1) < 1e-5)
        assert(abs(model.ipopt_zL_out[model.x[2]] - 1) < 1e-5)
        assert(abs(model.ipopt_zU_out[model.x[1]] - 1) < 1e-5)
        assert(abs(model.ipopt_zU_out[model.x[2]] - 1) < 1e-5)
        assert(abs(model.x[1].lb + 4) < 1e-5)

    @pytest.mark.unit
    def test09(self):
        """Try just saving suffixes, and suffix filter"""
        model = self.setup_model02()

        model.dual[model.g] = 1
        model.ipopt_zL_out[model.x[1]] = 1
        model.ipopt_zL_out[model.x[2]] = 1
        model.ipopt_zU_out[model.x[1]] = 1
        model.ipopt_zU_out[model.x[2]] = 1

        wts = StoreSpec.suffix(suffix_filter=("dual",))
        to_json(model, fname=self.fname, wts=wts)

        model.dual[model.g] = 10
        model.ipopt_zL_out[model.x[1]] = 10
        model.ipopt_zL_out[model.x[2]] = 10
        model.ipopt_zU_out[model.x[1]] = 10
        model.ipopt_zU_out[model.x[2]] = 10

        from_json(model, fname=self.fname, wts=wts)
        assert(abs(model.dual[model.g] - 1) < 1e-5)
        assert(abs(model.ipopt_zL_out[model.x[1]] - 10) < 1e-5)
        assert(abs(model.ipopt_zL_out[model.x[2]] - 10) < 1e-5)
        assert(abs(model.ipopt_zU_out[model.x[1]] - 10) < 1e-5)
        assert(abs(model.ipopt_zU_out[model.x[2]] - 10) < 1e-5)

    @pytest.mark.unit
    def test10(self):
        """Try just saving suffixes, and suffix filter only on write"""
        model = self.setup_model02()

        model.dual[model.g] = 1
        model.ipopt_zL_out[model.x[1]] = 1
        model.ipopt_zL_out[model.x[2]] = 1
        model.ipopt_zU_out[model.x[1]] = 1
        model.ipopt_zU_out[model.x[2]] = 1

        wts = StoreSpec.suffix(suffix_filter=("dual",))
        to_json(model, fname=self.fname, wts=wts)

        model.dual[model.g] = 10
        model.ipopt_zL_out[model.x[1]] = 10
        model.ipopt_zL_out[model.x[2]] = 10
        model.ipopt_zU_out[model.x[1]] = 10
        model.ipopt_zU_out[model.x[2]] = 10

        from_json(model, fname=self.fname, wts=StoreSpec.suffix())
        assert(abs(model.dual[model.g] - 1) < 1e-5)
        assert(abs(model.ipopt_zL_out[model.x[1]] - 10) < 1e-5)
        assert(abs(model.ipopt_zL_out[model.x[2]] - 10) < 1e-5)
        assert(abs(model.ipopt_zU_out[model.x[1]] - 10) < 1e-5)
        assert(abs(model.ipopt_zU_out[model.x[2]] - 10) < 1e-5)

    @pytest.mark.unit
    def test11(self):
        """Try just saving suffixes, and suffix filter only on read"""
        model = self.setup_model02()

        model.dual[model.g] = 1
        model.ipopt_zL_out[model.x[1]] = 1
        model.ipopt_zL_out[model.x[2]] = 1
        model.ipopt_zU_out[model.x[1]] = 1
        model.ipopt_zU_out[model.x[2]] = 1


        to_json(model, fname=self.fname, wts=StoreSpec.suffix())

        model.dual[model.g] = 10
        model.ipopt_zL_out[model.x[1]] = 10
        model.ipopt_zL_out[model.x[2]] = 10
        model.ipopt_zU_out[model.x[1]] = 10
        model.ipopt_zU_out[model.x[2]] = 10

        wts = StoreSpec.suffix(suffix_filter=("dual",))
        from_json(model, fname=self.fname, wts=wts)
        assert(abs(model.dual[model.g] - 1) < 1e-5)
        assert(abs(model.ipopt_zL_out[model.x[1]] - 10) < 1e-5)
        assert(abs(model.ipopt_zL_out[model.x[2]] - 10) < 1e-5)
        assert(abs(model.ipopt_zU_out[model.x[1]] - 10) < 1e-5)
        assert(abs(model.ipopt_zU_out[model.x[2]] - 10) < 1e-5)

if __name__ == '__main__':
    unittest.main()
