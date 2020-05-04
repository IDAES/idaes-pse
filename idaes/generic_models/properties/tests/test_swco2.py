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

__author__ = "John Eslick"

import pytest
from pyomo.environ import ConcreteModel, value, Var, SolverFactory
from pyomo.environ import ExternalFunction as EF
from pyomo.common.fileutils import this_file_dir
import idaes.generic_models.properties.swco2 as swco2
from idaes.generic_models.unit_models import Compressor
from idaes.core import FlowsheetBlock
import csv
import os
import idaes

_so = os.path.join(idaes.lib_directory, "swco2_external.so")

def swco2_available():
    """Make sure the compiled IAPWS-95 functions are available. Yes, in Windows
    the .so extension is still used.
    """
    return os.path.isfile(_so)

if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6}
else:
    solver = None

def read_data(fname):
    dfile = os.path.join(this_file_dir(), fname)
    data = {
        "T": [], # T in K col 0
        "P": [], # P in kPa col 1
        "rho": [], # density kg/m3 col 2
        "U": [], # internal energy kJ/kg col 4
        "H": [], # enthalpy kJ/kg col 5
        "S": [], # entropy kJ/kg/K col 6
        "cv": [],
        "cp": [],
        "w": [],
        "phase": [], # liquid, vapor, or supercritical col 13
        "visc": [],
        "tc": [],
    }

    with open(dfile, 'r') as csvfile:
        dat = csv.reader(csvfile, delimiter='\t', quotechar='"')
        for i in range(7):
            next(dat) # skip header
        for row in dat:
            data["T"].append(float(row[0]))
            data["P"].append(float(row[1])*1000)
            data["rho"].append(float(row[2]))
            data["U"].append(float(row[4])-506.7791289) # differnt reference state
            data["H"].append(float(row[5])-506.7791289) #different reference state
            data["S"].append(float(row[6])-2.739003) #differnt reference state
            data["cv"].append(float(row[7]))
            data["cp"].append(float(row[8]))
            data["w"].append(float(row[9]))
            data["visc"].append(float(row[11]))
            data["tc"].append(float(row[12]))
            data["phase"].append(row[13])
    return data

def read_sat_data(fname):
    dfile = os.path.join(this_file_dir(), fname)
    data = {}
    data["T"] = [] # T in K col 0
    data["P"] = [] # P in kPa col 1
    data["rhol"] = [] # density kg/m3 col 2
    data["rhov"] = [] # density kg/m3 col 15

    with open(dfile, 'r') as csvfile:
        dat = csv.reader(csvfile, delimiter='\t', quotechar='"')
        for i in range(7):
            next(dat) # skip header
        for row in dat:
            data["T"].append(float(row[0]))
            data["P"].append(float(row[1])*1000)
            data["rhol"].append(float(row[2]))
            data["rhov"].append(float(row[14]))
    return data


@pytest.mark.skipif(not swco2_available(), reason="Span-Wagner lib not available")
class TestSWCO2(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        plib = _so
        model.func_p = EF(library=plib, function="p")
        model.func_u = EF(library=plib, function="u")
        model.func_s = EF(library=plib, function="s")
        model.func_h = EF(library=plib, function="h")
        model.func_hvpt = EF(library=plib, function="hvpt")
        model.func_hlpt = EF(library=plib, function="hlpt")
        model.func_tau = EF(library=plib, function="tau")
        model.func_tau_sp = EF(library=plib, function="tau_sp")
        model.func_vf = EF(library=plib, function="vf")
        model.func_g = EF(library=plib, function="g")
        model.func_f = EF(library=plib, function="f")
        model.func_cv = EF(library=plib, function="cv")
        model.func_cp = EF(library=plib, function="cp")
        model.func_w = EF(library=plib, function="w")
        model.func_delta_liq = EF(library=plib, function="delta_liq")
        model.func_delta_vap = EF(library=plib, function="delta_vap")
        model.func_delta_sat_l = EF(library=plib, function="delta_sat_l")
        model.func_delta_sat_v = EF(library=plib, function="delta_sat_v")
        model.func_p_sat = EF(library=plib, function="p_sat")
        model.func_tau_sat = EF(library=plib, function="tau_sat")
        model.func_phi0 = EF(library=plib, function="phi0")
        model.func_phi0_delta = EF(library=plib, function="phi0_delta")
        model.func_phi0_delta2 = EF(library=plib, function="phi0_delta2")
        model.func_phi0_tau = EF(library=plib, function="phi0_tau")
        model.func_phi0_tau2 = EF(library=plib, function="phi0_tau2")
        model.func_phir = EF(library=plib, function="phir")
        model.func_phir_delta = EF(library=plib, function="phir_delta")
        model.func_phir_delta2 = EF(library=plib, function="phir_delta2")
        model.func_phir_tau = EF(library=plib, function="phir_tau")
        model.func_phir_tau2 = EF(library=plib, function="phir_tau2")
        model.func_phir_delta_tau = EF(library=plib, function="phir_delta_tau")
        return model

    def test_solve_sat_density(self, model):
        # test saturated liquid and vapor density solve (critical part of the
        # phase equlibrium calc)
        data = read_sat_data("sat_prop_swco2_nist_webbook.txt")
        for i, T in enumerate(data["T"]):
            rhol = value(467.6*model.func_delta_sat_l(304.128/T))
            rhov = value(467.6*model.func_delta_sat_v(304.128/T))
            if T > 296:
                tol = 1e-2 #data needs more sig fig
            else:
                tol = 1e-3
            assert rhol == pytest.approx(data["rhol"][i], rel=tol)
            assert rhov == pytest.approx(data["rhov"][i], rel=tol)

            # Use the data to check the regular delta(P, T) functions.
            # need not pressure sig figs for liquid, so have to reduce tol more
            if T > 296:
                tol = 1e-1 # data needs more sig fig
            rhol = value(467.6*model.func_delta_liq(data["P"][i], 304.128/T))
            rhov = value(467.6*model.func_delta_vap(data["P"][i], 304.128/T))
            assert rhol == pytest.approx(data["rhol"][i], rel=tol)
            assert rhov == pytest.approx(data["rhov"][i], rel=tol)

    def test_solve_tau_hp(self, model):
        data = read_data("prop_swco2_nist_webbook.txt")
        for i, T in enumerate(data["T"]):
            temp_sat = value(304.128/model.func_tau_sat(data["P"][i]))
            temp = value(304.128/model.func_tau(data["H"][i], data["P"][i]))
            print("h = {}, p = {}, T = {}, T_sat = {}".format(data["H"][i], data["P"][i], T, temp_sat))
            assert temp == pytest.approx(T, rel=0.1)

    def test_solve_tau_sp(self, model):
        data = read_data("prop_swco2_nist_webbook.txt")
        for i, T in enumerate(data["T"]):
            if data["phase"][i] not in ["vapor"]:
                continue
            print("s = {}, p = {}, T = {}".format(data["S"][i], data["P"][i], T))
            temp = value(304.128/model.func_tau_sp(data["S"][i], data["P"][i]))
            assert temp == pytest.approx(T, rel=1e-3)

    def test_solve_vapor_density(self, model):
        data = read_data("prop_swco2_nist_webbook.txt")
        for i, T in enumerate(data["T"]):
            if data["phase"][i] == "vapor" or data["phase"][i] == "supercritical":
                rho = value(467.6*model.func_delta_vap(data["P"][i], 304.128/T))
                assert rho == pytest.approx(data["rho"][i], rel=1e-2)

    def test_solve_liquid_density(self, model):
        data = read_data("prop_swco2_nist_webbook.txt")
        for i, T in enumerate(data["T"]):
            if data["phase"][i] == "liquid" or data["phase"][i] == "supercritical":
                rho = value(467.6*model.func_delta_liq(data["P"][i], 304.128/T))
                print(T, data["P"][i], data["phase"][i], rho)
                assert rho == pytest.approx(data["rho"][i], rel=1e-2)

    def test_functions_of_delta_and_tau(self, model):
        def tau(_T):
            return 304.128/_T
        def delta(_rho):
            return _rho/467.6

        def check(_T, _rho, func, val, rel=1e-4):
            print(func.name)
            val = pytest.approx(val, rel=rel)
            assert value(func(delta(_rho), tau(_T))) == val

        data = read_data("prop_swco2_nist_webbook.txt")
        for i, T in enumerate(data["T"]):
            check(T, data["rho"][i], func=model.func_p, val=data["P"][i], rel=1e-2)
            check(T, data["rho"][i], func=model.func_u, val=data["U"][i], rel=1e-2)
            check(T, data["rho"][i], func=model.func_h, val=data["H"][i], rel=1e-2)
            check(T, data["rho"][i], func=model.func_s, val=data["S"][i], rel=1e-2)
            check(T, data["rho"][i], func=model.func_cp, val=data["cp"][i], rel=1e-2)
            check(T, data["rho"][i], func=model.func_cv, val=data["cv"][i], rel=1e-2)
            check(T, data["rho"][i], func=model.func_w, val=data["w"][i], rel=1e-2)

    @pytest.fixture(scope="class")
    def model2(self):
        model = ConcreteModel()
        model.prop_param = swco2.SWCO2ParameterBlock()
        model.prop_in = swco2.SWCO2StateBlock(
            default={"parameters": model.prop_param}
        )
        return model

    def test_transport(self, model2):
        data = (
           #(T,      P,        mu,          tc,         tol_mu,  tol_tc)
            (200,    0.1e6,    10.06e-6,    9.63e-3,    0.02,    0.02),
            (240,    0.1e6,    12.07e-6,    12.23e-3,   0.02,    0.02),
            (240,    10.0e6,   188.91e-6,   159.24e-3,  0.02,    0.02),
            (240,    25e6,     214.16e-6,   171.59e-3,  0.02,    0.02),
            (280,    0.1e6,    14.05e-6,    15.19e-3,   0.02,    0.02),
            (280,    2.5e6,    14.51e-6,    17.46e-3,   0.02,    0.04),
            (280,    30.0e6,   134.98e-6,   136.63e-3,  0.02,    0.02),
            (280,    50.0e6,   160.73e-6,   152.67e-3,  0.02,    0.02),
            (306,    7.5e6,    24.06e-6,    23.93e-3,   0.02,    0.40), # supercritcal near critical
            (400,    5.0e6,    20.37e-6,    27.46e-3,   0.02,    0.02),
            (400,    50.0e6,   65.67e-6,    86.48e-3,   0.02,    0.02),
            (400,    100e6,    101.41e-6,   121.41e-3,  0.02,    0.02),
            (600,    50e6,     42.20e-6,    64.11e-3,   0.02,    0.02),
            (1000,   5e6,      41.42e-6,    71.18e-3,   0.02,    0.02),
            (1000,   50e6,     46.16e-6,    80.33e-3,   0.02,    0.02),
            (1000,   100e6,    54.79e-6,    92.16e-3,   0.02,    0.02),
        )

        for d in data:
            tol_mu = d[4]
            tol_tc = d[5]
            P = d[1]
            T = d[0]
            mu_data = d[2]
            tc_data = d[3]
            model2.prop_in.temperature.set_value(T)
            model2.prop_in.pressure = P
            Tsat = value(model2.prop_in.temperature_sat)
            if P >= 7.377e6 and T >= 304.1282:
                # super critical, which we clasify as liquid
                ph = "Liq"
            if P <= 0.5179e6:
                # below triple point
                ph = "Vap"
            elif T > Tsat + 0.5:
                ph = "Vap"
            elif T < Tsat - 0.5:
                ph = "Liq"
            else:
                # if too close to sat, don't want to get into a situation where
                # I don't know if it's liquid or vapor, so using extreme caution
                # If we're woried about it, we can add tests on the sat curve, and
                # test both liquid and vapor
                continue

            rho = value(model2.prop_in.dens_mass_phase[ph])
            mu = value(model2.prop_in.visc_d_phase[ph])
            tc = value(model2.prop_in.therm_cond_phase[ph])
            print("T = {}, P = {}, mu = {}, tc = {}, rho = {}, phase = {}".format(
                    T, P, mu, tc, rho, ph
                )
            )
            print("Data: mu = {}, tc = {}".format(mu_data, tc_data))
            assert tc == pytest.approx(tc_data, rel=tol_tc)
            assert mu == pytest.approx(mu_data, rel=tol_mu)


@pytest.mark.skipif(not swco2.swco2_available(), reason="Library not available")
class TestIntegration(object):
    @pytest.fixture(scope="class")
    def compressor_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = swco2.SWCO2ParameterBlock()
        m.fs.unit = Compressor(default={"property_package": m.fs.properties})
        return m

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_verify(self, compressor_model):
        model = compressor_model
        # Verify the turbine results against 3 known test cases

        # Case Data (90% isentropic efficency)
        cases = {
            "F": (1000, 1000), # mol/s
            "Tin": (500, 300), # K
            "Pin": (10, 100), # kPa
            "W": (3414.29266, 2796.30966), # kW
            "Tout": (574.744119, 372.6675676), # K
            "Pout": (20, 250), # kPa
            "xout": (1.0, 1.0), # vapor fraction
            "Tisen": (567.418852, 365.7680891),
        }

        for i , F in enumerate(cases["F"]):
            Tin = cases["Tin"][i]
            Tout = cases["Tout"][i]
            Pin = cases["Pin"][i]*1000
            Pout = cases["Pout"][i]*1000
            hin = swco2.htpx(T=Tin, P=Pin)
            W = cases["W"][i]*1000
            Tis = cases["Tisen"][i]
            xout = cases["xout"][i]

            model.fs.unit.inlet.flow_mol[0].fix(F)
            model.fs.unit.inlet.enth_mol[0].fix(hin)
            model.fs.unit.inlet.pressure[0].fix(Pin)
            model.fs.unit.deltaP.fix(Pout - Pin)
            model.fs.unit.efficiency_isentropic.fix(0.9)
            model.fs.unit.initialize(optarg={'tol': 1e-6})
            results = solver.solve(model)

            Tout = pytest.approx(cases["Tout"][i], rel=1e-2)
            Pout = pytest.approx(cases["Pout"][i]*1000, rel=1e-2)
            Pout = pytest.approx(cases["Pout"][i]*1000, rel=1e-2)
            W = pytest.approx(cases["W"][i]*1000, rel=1e-2)
            xout = pytest.approx(xout, rel=1e-2)
            prop_out = model.fs.unit.control_volume.properties_out[0]
            prop_in = model.fs.unit.control_volume.properties_in[0]
            prop_is = model.fs.unit.properties_isentropic[0]

            assert value(prop_in.temperature) == pytest.approx(Tin, rel=1e-3)
            assert value(prop_is.temperature) == pytest.approx(Tis, rel=1e-3)
            assert value(model.fs.unit.control_volume.work[0]) == W
            assert value(prop_out.pressure) == Pout
            assert value(prop_out.temperature) == Tout
            assert value(prop_out.vapor_frac) == xout
