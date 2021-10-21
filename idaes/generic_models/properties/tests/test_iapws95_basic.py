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
import idaes.generic_models.properties.iapws95 as iapws95
from idaes.generic_models.properties.iapws95 import \
    iapws95_available as prop_available
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
            data["U"].append(float(row[4])*mw*1000)
            data["H"].append(float(row[5])*mw*1000)
            data["S"].append(float(row[6])*mw*1000)
            #Some things are undefined at critical point.
            try:
                data["visc"].append(float(row[11]))
            except:
                data["visc"].append(None)
            try:
                data["tc"].append(float(row[12]))
            except:
                data["tc"].append(None)
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
    """Test derivatives for function f(x0, x1) against f.d. approx
    (with assert)

    Args:
        f: ExternalFunction to test
        x0: f argument 1
        x1: f argument 2
        d0: f.d. step for grad[0] (x0)
        d1: f.d. step for grad[1] (x1)
        tol: assert derivitive value tolerance
    """
    assert(isinstance(f, AMPLExternalFunction))
    y, g, h = f.evaluate_fgh(args=(x0, x1))
    yf0, gf0, hf0 = f.evaluate_fgh(args=(x0 + d0, x1))
    yb0, gb0, hb0 = f.evaluate_fgh(args=(x0 - d0, x1))
    yf1, gf1, hf1 = f.evaluate_fgh(args=(x0, x1 + d1))
    yb1, gb1, hb1 = f.evaluate_fgh(args=(x0, x1 - d1))
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


class TestHelm(object):
    mw = 0.01801528
    Tc = 647.096
    Pc = 2.2064e7  # Pa
    rhoc = 322  # kg/m3
    Pmin = 1000  # Pa
    Pmax = 20*Pc  # Pa
    Tmax = 1000  # K
    Tmin = 274
    pparam = iapws95
    pparam_construct = iapws95.Iapws95ParameterBlock
    pdata = "prop_iapws95_nist_webbook.txt"
    pdata_sat = "sat_prop_iapws95_nist_webbook.txt"
    onein = 10

    @pytest.fixture(scope="class")
    def model_transport(self):
        # This model is used to test transport properties
        model = ConcreteModel()
        model.prop_param = iapws95.Iapws95ParameterBlock()
        model.prop_in = iapws95.Iapws95StateBlock(
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

    def test_thermo_expression_writter(self, model):
        te = model.te

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
            if i % self.onein:
                continue

            # Test state variable with P in the set, these are pretty reliable
            assert value(te.s(h=h, p=p)) == pytest.approx(
                value(s), rel=0.05)
            assert value(te.h(u=u, p=p)) == pytest.approx(
                value(h), rel=0.05)
            assert value(te.h(T=T*pyunits.K, p=p, x=x)) == pytest.approx(
                value(h), rel=0.05)
            assert value(te.h(s=s, p=p)) == pytest.approx(
                value(h), rel=0.05)
            if value(p) < 2e6 or T < 290:
                # need data with more significant figure to test here
                # generally p(s, T) prbably isn't that useful for liquids
                # so I'll come back to it later with high precision data.
                continue
            assert value(te.p(s=s, T=T*pyunits.K)) == pytest.approx(
                value(p), rel=0.1)

            # Commenting out the deriative test for now.  The derivatives are
            # tested in the CO2 tests and are the same for all Helmholtz EOSs
            # running all these is pretty time consuming since there are so
            # many data points for water.

            # test the deriviatives that are critical to the thermo expressions
            #binary_derivative_test(f=model.func_p_stau, x0=s/mw/1000, x1=Tc/T)
            #binary_derivative_test(f=model.func_tau, x0=h/mw/1000, x1=p/1000)
            #binary_derivative_test(f=model.func_tau_sp, x0=s/mw/1000, x1=p/1000)
            #binary_derivative_test(f=model.func_tau_up, x0=u/mw/1000, x1=p/1000)
            #binary_derivative_test(f=model.func_vf, x0=h/mw/1000, x1=p/1000)
            #binary_derivative_test(f=model.func_vfs, x0=s/mw/1000, x1=p/1000)
            #binary_derivative_test(f=model.func_vfu, x0=u/mw/1000, x1=p/1000)

    def test_solve_vapor_density(self, model):
        """ The density calculations should be tested by the thermo expression
        tests, but they are pretty fundimental to everything else, so test them
        a little more on their own.
        """
        te = model.te
        data = read_data(self.pdata, self.mw)
        for i, T in enumerate(data["T"]):
            if data["phase"][i] == "vapor":
                rho = value(te.rho_vap(
                    p=data["P"][i]*pyunits.Pa, T=T*pyunits.K, x=1))
                assert rho == pytest.approx(value(data["rho"][i]), rel=1e-2)

    def test_solve_liquid_density(self, model):
        """ The density calculations should be tested by the thermo expression
        tests, but they are pretty fundimental to everything else, so test them
        a little more on their own.
        """
        te = model.te
        data = read_data(self.pdata, self.mw)
        for i, T in enumerate(data["T"]):
            if data["phase"][i] == "liquid":
                rho = value(te.rho_liq(
                    p=data["P"][i]*pyunits.Pa, T=T*pyunits.K, x=0))
                print("T {}, P {}, rho dat {}, rho {}".format(
                    T, data["P"][i], data["rho"][i], rho))
                assert rho == pytest.approx(value(data["rho"][i]), rel=1e-1)

    def test_solve_supercritical_density(self, model):
        """ The density calculations should be tested by the thermo expression
        tests, but they are pretty fundimental to everything else, so test them
        a little more on their own.
        """
        te = model.te
        data = read_data(self.pdata, self.mw)
        for i, T in enumerate(data["T"]):
            if data["phase"][i] == "supercritical":
                rhol = value(
                    te.rho_liq(p=data["P"][i]*pyunits.Pa, T=T*pyunits.K, x=0))
                rhov = value(
                    te.rho_vap(p=data["P"][i]*pyunits.Pa, T=T*pyunits.K, x=0))
                assert rhol == pytest.approx(value(data["rho"][i]), rel=0.5e-1)
                assert rhov == pytest.approx(value(data["rho"][i]), rel=0.5e-1)

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
            if T > self.Tc - 0.1:
                # if this goes over the critical temperature this makes no
                # sense
                # also really close to the critical temperature for various
                # reasons we'll diverge a bit from what we are comparing to
                pass
            else:
                tol = 1e-2
                # test p, x spec
                rhol = value(te.rho_liq(
                    p=data["P"][i]*pyunits.Pa, x=0))
                rhov = value(te.rho_vap(
                    p=data["P"][i]*pyunits.Pa, x=1))
                assert rhol == pytest.approx(value(data["rhol"][i]), rel=tol)
                assert rhov == pytest.approx(value(data["rhov"][i]), rel=tol)
                # test t, x spec
                rhol = value(te.rho_liq(T=T*pyunits.K, x=0))
                rhov = value(te.rho_vap(T=T*pyunits.K, x=1))
                assert rhol == pytest.approx(value(data["rhol"][i]), rel=tol)
                assert rhov == pytest.approx(value(data["rhov"][i]), rel=tol)

            # Ignore the phase equilibrium and use T,P data to calc densities
            if T > 296:
                tol = 1e-1  # data needs more sig fig
            rhol = value(
                te.rho_liq(p=data["P"][i]*pyunits.Pa, T=T*pyunits.K, x=0))
            rhov = value(
                te.rho_vap(p=data["P"][i]*pyunits.Pa, T=T*pyunits.K, x=1))
            assert rhol == pytest.approx(value(data["rhol"][i]), rel=tol)
            assert rhov == pytest.approx(value(data["rhov"][i]), rel=tol)

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
        mw = 0.01801528
        for i, T in enumerate(data["T"]):
            p = data["P"][i]/1000
            u = data["U"][i]/1000/mw
            s = data["S"][i]/1000/mw
            h = data["H"][i]/1000/mw
            rho = data["rho"][i]

            if p < self.Pmin or p > self.Pmax:
                continue
            if T < self.Tmin or T > self.Tmax:
                continue

            check(T, rho, func=model.func_p, val=p, rel=1e-01)
            check(T, rho, func=model.func_u, val=u, rel=1e-01)
            check(T, rho, func=model.func_h, val=h, rel=1e-01)
            check(T, rho, func=model.func_s, val=s, rel=1e-01)

            binary_derivative_test(f=model.func_p, x0=delta(rho), x1=tau(T))
            binary_derivative_test(f=model.func_u, x0=delta(rho), x1=tau(T))
            binary_derivative_test(f=model.func_s, x0=delta(rho), x1=tau(T))
            binary_derivative_test(f=model.func_h, x0=delta(rho), x1=tau(T))
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
            if abs(self.Pc - data["P"][i]) < 3e6 and abs(self.Tc - T) < 25:
                #undefined at critical point and vers sensitive close to it.
                continue
            print(f"T = {T}, P = {data['P'][i]}, TC = {data['tc'][i]}, Visc = {data['visc'][i]}")
            assert tc == pytest.approx(data["tc"][i], rel=5e-1)
            assert mu == pytest.approx(data["visc"][i], rel=5e-1)
