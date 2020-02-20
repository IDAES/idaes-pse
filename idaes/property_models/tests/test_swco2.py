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
from pyomo.environ import ConcreteModel, value, Var
from pyomo.environ import ExternalFunction as EF
from pyomo.common.fileutils import this_file_dir
import csv
import os
import idaes

_so = os.path.join(idaes.lib_directory, "swco2_external.so")

def swco2_available():
    """Make sure the compiled IAPWS-95 functions are available. Yes, in Windows
    the .so extention is still used.
    """
    return os.path.isfile(_so)

def read_data(fname):
    dfile = os.path.join(this_file_dir(), fname)
    data = {}
    data["T"] = [] # T in K col 0
    data["P"] = [] # P in kPa col 1
    data["rho"] = [] # density kg/m3 col 2
    data["U"] = [] # internal energy kJ/kg col 4
    data["H"] = [] # enthalpy kJ/kg col 5
    data["S"] = [] # entropy kJ/kg/K col 6
    data["cv"] = []
    data["cp"] = []
    data["w"] = []
    data["phase"] = [] # liquid, vapor, or supercirtical col 13

    with open(dfile, 'r') as csvfile:
        dat = csv.reader(csvfile, delimiter='\t', quotechar='"')
        next(dat)  # skip header
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
        next(dat)  # skip header
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
        #test saturated liquid and vapor density solve (ciritical part of the
        #phase equlibrium calc)
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
            # unfotunatly the pressure reported in the data doesn't have enough
            # figures for a really accurate check
            rhol = value(467.6*model.func_delta_liq(data["P"][i], 304.128/T))
            rhov = value(467.6*model.func_delta_vap(data["P"][i], 304.128/T))
            if T > 285:
                tol = 1e-1 #data needs more sig fig
            print(T, rhol, rhov)
            assert rhol == pytest.approx(data["rhol"][i], rel=tol)
            assert rhov == pytest.approx(data["rhov"][i], rel=tol)


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
            #print(func.name)
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
