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

__author__ = "John Eslick"

import pytest
from pyomo.environ import ConcreteModel, value, units as pyunits
from pyomo.common.fileutils import this_file_dir
from pyomo.core.base.external import AMPLExternalFunction
import idaes.generic_models.properties.swco2 as swco2
import csv
import os
import idaes


# Mark module as an integration test
pytestmark = pytest.mark.integration


def read_data(fname, mw):
    dfile = os.path.join(this_file_dir(), fname)
    data = {
        "T": [],  # T in K col 0
        "P": [],  # P in kPa col 1
        "rho": [],  # density kg/m3 col 2
        "U": [],  # internal energy kJ/kg col 4
        "H": [],  # enthalpy kJ/kg col 5
        "S": [],  # entropy kJ/kg/K col 6
        "cv": [],
        "cp": [],
        "w": [],
        "phase": [],  # liquid, vapor, or supercritical col 13
        "visc": [],
        "tc": [],
    }

    with open(dfile, 'r') as csvfile:
        dat = csv.reader(csvfile, delimiter='\t', quotechar='"')
        for i in range(7):
            next(dat)  # skip header
        for row in dat:
            data["T"].append(float(row[0]))
            data["P"].append(float(row[1])*1e6)
            data["rho"].append(float(row[2]))

            # Different references states
            data["U"].append((float(row[4])-506.7791289)*mw*1000)
            data["H"].append((float(row[5])*mw*1000) + -22303.24810876697)
            data["S"].append((float(row[6])-2.739003)*mw*1000)

            data["cv"].append(float(row[7]))
            data["cp"].append(float(row[8]))
            data["w"].append(float(row[9]))
            data["visc"].append(float(row[11]))
            data["tc"].append(float(row[12]))
            data["phase"].append(row[13])
    return data


def read_sat_data(fname, mw):
    dfile = os.path.join(this_file_dir(), fname)
    data = {}
    data["T"] = []  # T in K col 0
    data["P"] = []  # P in kPa col 1
    data["rhol"] = []  # density kg/m3 col 2
    data["rhov"] = []  # density kg/m3 col 15

    with open(dfile, 'r') as csvfile:
        dat = csv.reader(csvfile, delimiter='\t', quotechar='"')
        for i in range(7):
            next(dat)  # skip header
        for row in dat:
            data["T"].append(float(row[0]))
            data["P"].append(float(row[1])*1e6)
            data["rhol"].append(float(row[2]))
            data["rhov"].append(float(row[14]))
    return data


def between(y, x0, x1):
    return 0 > (y-x0)*(y-x1)


def unary_derivative_test(f, x0, d=1e-5, tol=0.02):
    """Test derivatives for function f(x) against f.d. approx (with assert)

    Args:
        f: ExternalFunction to test
        x: f argument 1
        d: f.d. step for grad
        tol: assert derivitive value tolerance
    """
    assert(isinstance(f, AMPLExternalFunction))
    y, g, h = f.evaluate_fgh(args=(x0,))
    yf, gf, hf = f.evaluate_fgh(args=(x0 + d,))
    yb, gb, hb = f.evaluate_fgh(args=(x0 - d,))
    gfdf = (yf - y)/d
    gfdb = -(yb - y)/d
    hfdf = (gf[0] - g[0])/d
    hfdb = -(gb[0] - g[0])/d

    zero_cut = 1e-9  # how close to zero before maybe it is zero?

    # check that the forward and backward FD approximations are close enough
    # that the accuracy is good enough for the test and that the detivative
    # is not ~ zero.  I know this rough but what can you do?
    if abs(g[0]) > zero_cut:
        assert (abs((gfdf - g[0])/g[0]) < tol or
                abs((gfdb - g[0])/g[0]) < tol or
                between(g[0], gfdf, gfdb))

    if abs(h[0]) > zero_cut:
        assert (abs((hfdf - h[0])/h[0]) < tol or
                abs((hfdb - h[0])/h[0]) < tol or
                between(h[0], hfdf, hfdb))


def binary_derivative_test(f, x0, x1, d0=1e-5, d1=1e-5, tol=0.02):
    """
    Test derivatives for function f(x0, x1) against f.d. approx (with assert)

    Args:
        f: ExternalFunction to test
        x0: f argument 1
        x1: f argument 2
        d0: f.d. step for grad[0] (x0)
        d1: f.d. step for grad[1] (x1)
        tol: assert derivitive value tolerance
    """
    assert(isinstance(f, AMPLExternalFunction))
    y, g, h = f.evaluate_fgh(args=(value(x0), value(x1)))
    yf0, gf0, hf0 = f.evaluate_fgh(args=(value(x0 + d0), value(x1)))
    yb0, gb0, hb0 = f.evaluate_fgh(args=(value(x0 - d0), value(x1)))
    yf1, gf1, hf1 = f.evaluate_fgh(args=(value(x0), value(x1 + d1)))
    yb1, gb1, hb1 = f.evaluate_fgh(args=(value(x0), value(x1 - d1)))
    gf = [(yf0 - y)/d0, (yf1 - y)/d1]
    gb = [-(yb0 - y)/d0, -(yb1 - y)/d1]
    hf = [(gf0[0] - g[0])/d0, (gf0[1] - g[1])/d0, (gf1[1] - g[1])/d1]
    hb = [-(gb0[0] - g[0])/d0, -(gb0[1] - g[1])/d0, -(gb1[1] - g[1])/d1]

    zero_cut = 1e-9  # how close to zero before maybe it is zero?
    # check that the forward and backward FD approximations are close enough
    # that the accuracy is good enough for the test and that the detivative
    # is not ~ zero.  I know this rough but what can you do?

    if abs(g[0]) > zero_cut:  # derivative is not 0
        assert(abs((gf[0] - g[0])/g[0]) < tol or
               abs((gb[0] - g[0])/g[0]) < tol or
               between(g[0], gf[0], gb[0]))

    if abs(g[1]) > zero_cut and abs((gf[1] - gb[1])/g[1]) < tol:
        assert (abs((gf[1] - g[1])/g[1]) < tol or
                abs((gb[1] - g[1])/g[1]) < tol or
                between(g[1], gf[1], gb[1]))

    if abs(h[0]) > zero_cut and abs((hf[0] - hb[0])/h[0]) < tol:
        assert (abs((hf[0] - h[0])/h[0]) < tol or
                abs((hb[0] - h[0])/h[0]) < tol or
                between(h[0], hf[0], hb[0]))

    if abs(h[1]) > zero_cut and abs((hf[1] - hb[1])/h[1]) < tol:
        assert (abs((hf[1] - h[1])/h[1]) < tol or
                abs((hb[1] - h[1])/h[1]) < tol or
                between(h[1], hf[1], hb[1]))

    if abs(h[2]) > zero_cut and abs((hf[2] - hb[2])/h[2]) < tol:
        assert (abs((hf[2] - h[2])/h[2]) < tol or
                abs((hb[2] - h[2])/h[2]) < tol or
                between(h[2], hf[2], hb[2]))


@pytest.mark.skipif(not swco2.swco2_available(),
                    reason="Property lib not available")
class TestHelm(object):
    mw = 0.0440098
    Tc = 304.128
    Pc = 7377300  # Pa
    rhoc = 467.6  # kg/m3
    Pmin = 1000  # Pa
    Pmax = 100*Pc  # Pa
    Tmax = 800  # K
    Tmin = 270
    pparam = swco2
    pparam_construct = swco2.SWCO2ParameterBlock
    pdata = "prop_swco2_nist_webbook.txt"
    pdata_sat = "sat_prop_swco2_nist_webbook.txt"

    @pytest.fixture(scope="class")
    def model_transport(self):
        # This model is used to test transport properties
        model = ConcreteModel()
        model.prop_param = swco2.SWCO2ParameterBlock()
        model.prop_in = swco2.SWCO2StateBlock(
            default={"parameters": model.prop_param}
        )
        return model

    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.prop = self.pparam_construct()
        model.te = self.pparam.HelmholtzThermoExpressions(
            model, parameters=model.prop)
        return model

    def tst_hs_state(self, model):
        """test development of a p(h, s) function """
        mw = self.mw
        data = read_data(self.pdata, self.mw)
        model.te.add_funcs(
            names=[
                "func_h", "func_s"
            ]
        )

        from scipy.optimize import fsolve

        data = read_data(self.pdata, self.mw)
        for i, T in enumerate(data["T"]):
            p = data["P"][i]
            h = data["H"][i]
            s = data["S"][i]
            rho = data["rho"][i]

            if p < self.Pmin or p > self.Pmax:
                continue
            if T < self.Tmin or T > self.Tmax:
                continue

            def f(a):
                _delta = a[0]
                _tau = a[1]
                fh = (value(model.func_h(_delta, _tau)) - h/mw/1000)
                fs = (value(model.func_s(_delta, _tau)) - s/mw/1000)
                return fh, fs

            x0 = [rho/self.rhoc, self.Tc/T]
            print("solve")
            print(f(x0))
            xs = fsolve(f, x0)
            print(xs)
            print(f(xs))
            print("T = {}, {}".format(self.Tc/xs[1], T))
            # assert sum(map(abs, f)) < 1e-5
        assert False

    def test_external_memo(self, model):
        """ This tests the memoization in the external functions.  There is a
        special set of functions that return the memoized value of a function.
        This test should catch things like functions that aren't memoized, or
        function with the wrong lookup tag getting or setting the wrong value.
        """
        mw = self.mw
        Tc = self.Tc

        model.te.add_funcs(
            names=[
                "memo_test_tau",
            ]
        )

        data = read_data(self.pdata, self.mw)
        il = None
        ig = None
        for i, phase in enumerate(data["phase"]):
            if phase == "vapor" and ig is None:
                ig = i
            elif phase == "liquid" and il is None:
                il = i

        for i in [il, ig]:
            T = data["T"][i]
            p = data["P"][i]
            h = data["H"][i]
            phase = data["phase"][i]

            assert (value(model.memo_test_tau(h/mw/1000, p/1000)) ==
                    pytest.approx(Tc/T, rel=0.10))
            binary_derivative_test(
                f=model.memo_test_tau, x0=h/mw/1000, x1=p/1000)

    def test_thermo_expression_writter(self, model):
        te = model.te
        mw = self.mw
        Tc = self.Tc

        model.te.add_funcs(
            names=[
                "memo_test_tau",
            ]
        )

        data = read_data(self.pdata, self.mw)
        for i, T in enumerate(data["T"]):
            p = data["P"][i]*pyunits.Pa
            h = data["H"][i]*pyunits.J/pyunits.mol
            s = data["S"][i]*pyunits.J/pyunits.mol/pyunits.K
            u = data["U"][i]*pyunits.J/pyunits.mol
            if data["phase"][i] == "vapor":
                x = 1
            else:
                x = 0

            if value(p) < self.Pmin or value(p) > self.Pmax:
                continue
            if T < self.Tmin or T > self.Tmax:
                continue

            # if it fails I want to know where so pring T and P
            print("T = {}, P = {}".format(T, p))

            # Test state variable with P in the set, these are pretty reliable
            assert value(te.s(h=h, p=p)) == pytest.approx(
                value(s), rel=0.05)
            assert value(te.s(u=u, p=p)) == pytest.approx(
                value(s), rel=0.05)
            assert value(te.h(T=T*pyunits.K, p=p, x=x)) == pytest.approx(
                value(h), rel=0.05)
            assert value(te.h(s=s, p=p)) == pytest.approx(
                value(h), rel=0.05)
            assert value(te.p(s=s, T=T*pyunits.K)) == pytest.approx(
                value(p), rel=0.05)

            # test the deriviatives that are critical to the thermo expressions
            binary_derivative_test(f=model.func_p_stau, x0=s/mw/1000, x1=Tc/T)
            binary_derivative_test(f=model.func_tau, x0=h/mw/1000, x1=p/1000)
            binary_derivative_test(
                f=model.func_tau_sp, x0=s/mw/1000, x1=p/1000)
            binary_derivative_test(
                f=model.func_tau_up, x0=u/mw/1000, x1=p/1000)
            binary_derivative_test(f=model.func_vf, x0=h/mw/1000, x1=p/1000)
            binary_derivative_test(f=model.func_vfs, x0=s/mw/1000, x1=p/1000)
            binary_derivative_test(f=model.func_vfu, x0=u/mw/1000, x1=p/1000)

    def test_solve_vapor_density(self, model):
        """ The density calculations should be tested by the thermo expression
        tests, but they are pretty fundimental to everything else, so test them
        a little more on their own.
        """
        te = model.te
        data = read_data(self.pdata, self.mw)
        for i, T in enumerate(data["T"]):
            if (data["phase"][i] == "vapor" or
                    data["phase"][i] == "supercritical"):
                rho = value(te.rho_vap(
                    p=data["P"][i]*pyunits.Pa, T=T*pyunits.K, x=1))
                assert rho == pytest.approx(data["rho"][i], rel=1e-2)

    def test_solve_liquid_density(self, model):
        """ The density calculations should be tested by the thermo expression
        tests, but they are pretty fundimental to everything else, so test them
        a little more on their own.
        """
        te = model.te
        data = read_data(self.pdata, self.mw)
        for i, T in enumerate(data["T"]):
            if (data["phase"][i] == "liquid" or
                    data["phase"][i] == "supercritical"):
                rho = value(te.rho_liq(
                    p=data["P"][i]*pyunits.Pa, T=T*pyunits.K, x=0))
                assert rho == pytest.approx(data["rho"][i], rel=1e-2)

    def test_solve_sat_density(self, model):
        """ The density calculations should be tested by the thermo expression
        tests, but they are pretty fundimental to everything else, so test them
        a little more on their own.
        """
        # test saturated liquid and vapor density solve (critical part of the
        # phase equlibrium calc)
        te = model.te
        data = read_sat_data(self.pdata_sat, self.mw)
        for i, T in enumerate(data["T"]):
            if T > 304.128:
                # if this goes over the critical temperature this makes no
                # sense while we're looking at the two phase region. (not sure
                # how it got in the data)
                pass
            else:
                tol = 1e-2
                # test p, x spec
                rhol = value(
                    te.rho_liq(p=data["P"][i]*pyunits.Pa, x=0))
                rhov = value(
                    te.rho_vap(p=data["P"][i]*pyunits.Pa, x=1))
                assert rhol == pytest.approx(data["rhol"][i], rel=tol)
                assert rhov == pytest.approx(data["rhov"][i], rel=tol)
                # test t, x spec
                rhol = value(te.rho_liq(T=T*pyunits.K, x=0))
                rhov = value(te.rho_vap(T=T*pyunits.K, x=1))
                assert rhol == pytest.approx(data["rhol"][i], rel=tol)
                assert rhov == pytest.approx(data["rhov"][i], rel=tol)

            # Ignore the phase equilibrium and use T,P data to calc densities
            if T > 296:
                tol = 1e-1  # data needs more sig fig
            rhol = value(te.rho_liq(
                p=data["P"][i]*pyunits.Pa, T=T*pyunits.K, x=0))
            rhov = value(te.rho_vap(
                p=data["P"][i]*pyunits.Pa, T=T*pyunits.K, x=1))
            assert rhol == pytest.approx(data["rhol"][i], rel=tol)
            assert rhov == pytest.approx(data["rhov"][i], rel=tol)

    def test_functions_of_delta_and_tau(self, model):
        """
        These are the bisic direct from density and temperature propery
        calculation tests.
        """
        def tau(_T):
            return self.Tc/_T

        def delta(_rho):
            return _rho/self.rhoc

        def check(_T, _rho, func, val, rel=1e-4):
            val = pytest.approx(val, rel=rel)
            assert value(func(delta(_rho), tau(_T))) == val

        model.te.add_funcs(
            names=[
                "func_p",
                "func_u",
                "func_h",
                "func_s",
                "func_cp",
                "func_cv",
                "func_w",
                "func_g",
                "func_f",
            ]
        )

        data = read_data(self.pdata, self.mw)
        mw = self.mw
        for i, T in enumerate(data["T"]):
            p = data["P"][i]/1000
            u = data["U"][i]/1000/mw
            s = data["S"][i]/1000/mw
            h = data["H"][i]/1000/mw
            cp = data["cp"][i]
            cv = data["cv"][i]
            w = data["w"][i]
            rho = data["rho"][i]

            check(T, rho, func=model.func_p, val=p, rel=1e-2)
            check(T, rho, func=model.func_u, val=u, rel=1e-2)
            check(T, rho, func=model.func_h, val=h, rel=1e-2)
            check(T, rho, func=model.func_s, val=s, rel=1e-2)
            check(T, rho, func=model.func_cp, val=cp, rel=1e-2)
            check(T, rho, func=model.func_cv, val=cv, rel=1e-2)
            check(T, rho, func=model.func_w, val=w, rel=1e-2)

            if p < 1000:
                continue

            binary_derivative_test(f=model.func_p, x0=delta(rho), x1=tau(T))
            binary_derivative_test(f=model.func_u, x0=delta(rho), x1=tau(T))
            binary_derivative_test(f=model.func_s, x0=delta(rho), x1=tau(T))
            binary_derivative_test(f=model.func_h, x0=delta(rho), x1=tau(T))
            binary_derivative_test(f=model.func_cp, x0=delta(rho), x1=tau(T))
            binary_derivative_test(f=model.func_cv, x0=delta(rho), x1=tau(T))
            binary_derivative_test(f=model.func_w, x0=delta(rho), x1=tau(T))
            binary_derivative_test(f=model.func_f, x0=delta(rho), x1=tau(T))
            binary_derivative_test(f=model.func_g, x0=delta(rho), x1=tau(T))

    def test_transport(self, model_transport):
        """ Test transport properties.  The tolerances are pretty forgiving here.
        The values are closer for the most part, but the estimation methods
        aren't the same or are super sensitive near the critical point.  Just
        want a sanity check not for high accuracy.
        """
        m = model_transport
        data = read_data(self.pdata, self.mw)
        for i, T in enumerate(data["T"]):
            m.prop_in.temperature.set_value(T)
            m.prop_in.pressure = data["P"][i]
            ph = {"vapor":"Vap", "liquid":"Liq", "supercritical":"Liq"}[data["phase"][i]]
            mu = value(m.prop_in.visc_d_phase[ph])
            tc = value(m.prop_in.therm_cond_phase[ph])
            assert tc == pytest.approx(data["tc"][i], rel=2.5e-1)
            assert mu == pytest.approx(data["visc"][i], rel=2.5e-1)
