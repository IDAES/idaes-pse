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
from __future__ import division, print_function, absolute_import

import pytest
from pyomo.environ import *
from pyomo.opt import SolverFactory
from idaes.property_models import iapws95_ph as iapws95
from idaes.property_models.iapws95 import is_available as iapws_available
import csv
import math
import os


def read_data(fname, col):
    dfile = os.path.dirname(__file__)
    dfile = os.path.join(dfile, fname)
    cond = []  # Tuple (T [K],P [Pa], data) pressure in file is MPa
    with open(dfile, 'r') as csvfile:
        dat = csv.reader(csvfile, delimiter='\t', quotechar='"')
        next(dat)  # skip header
        for row in dat:
            try:
                x = float(row[col])
            except:
                x = row[col]
            cond.append((float(row[0]), float(row[1]) * 1e6, x))
    return cond

@pytest.mark.skipif(not iapws_available(), reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_limits():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(
        default={"parameters":model.prop_param})
    x = value(model.prop_in.func_delta_liq(1e5, -1000))
    assert(math.isnan(x))
    x = value(model.prop_in.func_delta_liq(1e5, 100))
    assert(math.isnan(x))
    x = value(model.prop_in.func_delta_liq(1e5, 1e7))
    assert(math.isnan(x))
    x = value(model.prop_in.func_delta_vap(1e5, -1000))
    assert(math.isnan(x))
    x = value(model.prop_in.func_delta_vap(1e5, 100))
    assert(math.isnan(x))
    x = value(model.prop_in.func_delta_vap(1e5, 1e7))
    assert(math.isnan(x))

@pytest.mark.skipif(not iapws_available(), reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_tau_sat():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("sat_prop.txt", col=2)
    for c in cond:
        tau = value(model.prop_in.func_tau_sat(c[1]/1000.0))
        T = 647.096/tau
        print("{}, {}, {}".format(c[1], c[0], T))
        assert(abs(T-c[0]) < 0.1)

@pytest.mark.skipif(not iapws_available(), reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_liquid_density_sat():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("sat_prop.txt", col=2)
    for c in cond:
        model.prop_in.temperature.set_value(c[0])
        psat = value(model.prop_in.pressure_sat)
        model.prop_in.pressure = psat
        rho = value(model.prop_in.dens_mass_phase["Liq"])
        print("{}, {}, {}, {}, {}".format(c[0], psat, rho, c[1], c[2]))
        if c[0] < 644:
            assert(abs(rho-c[2]) < 0.4)
        else:
            assert(abs(rho-c[2]) < 5)

@pytest.mark.skipif(not iapws_available(), reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_vapor_density_sat():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("sat_prop.txt", col=14)
    for c in cond:
        model.prop_in.temperature.set_value(c[0])
        psat = value(model.prop_in.pressure_sat)
        model.prop_in.pressure = psat
        rho = value(model.prop_in.dens_mass_phase["Vap"])
        print("{}, {}, {}, {}, {}".format(c[0], psat, rho, c[1], c[2]))
        if c[0] < 645:
            assert(abs(rho-c[2])/c[2]*100 < 1.0)

@pytest.mark.skipif(not iapws_available(), reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_density():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("prop.txt", col=2)
    phase = read_data("prop.txt", col=13)
    for i, c in enumerate(cond):
        if phase[i][2] in ["liquid", "supercritical"]:
            p = "Liq"
        else:
            p = "Vap"
        model.prop_in.temperature.set_value(c[0])
        model.prop_in.pressure = c[1]
        rho = value(model.prop_in.dens_mass_phase[p])
        print("{}, {}, {}, {}, {}".format(c[0], c[1], rho, c[2], p))
        if c[0] - 646.86 < 0.01 and c[1] - 22000000 < 1:
            assert(abs(rho-c[2])/c[2]*100 < 10.0) # have to look into this
        else:
            assert(abs(rho-c[2])/c[2]*100 < 1.0)

@pytest.mark.skipif(not iapws_available(), reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_enthalpy():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("prop.txt", col=5)
    phase = read_data("prop.txt", col=13)
    for i, c in enumerate(cond):
        if phase[i][2] in ["liquid", "supercritical"]:
            p = "Liq"
        else:
            p = "Vap"
        model.prop_in.temperature.set_value(c[0])
        model.prop_in.pressure = c[1]
        h = value(model.prop_in.enth_mol_phase[p]/model.prop_in.mw/1000)
        print("{}, {}, {}, {}, {}".format(c[0], c[1], h, c[2], p))
        assert(abs(h-c[2])/c[2]*100 < 1)

@pytest.mark.skipif(not iapws_available(), reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_enthalpy_vapor_as_function_of_p_and_tau():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("prop.txt", col=5)
    phase = read_data("prop.txt", col=13)
    for i, c in enumerate(cond):
        if phase[i][2] in ["liquid", "supercritical"]:
            continue
        model.prop_in.temperature.set_value(c[0])
        h = value(model.prop_in.func_hvpt(c[1]/1000, 647.096/c[0]))
        print("{}, {}, {}, {}".format(c[0], c[1], h, c[2]))
        assert(abs(h-c[2])/c[2]*100 < 1)

@pytest.mark.skipif(not iapws_available(), reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_enthalpy_liquid_as_function_of_p_and_tau():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("prop.txt", col=5)
    rho_dat = read_data("prop.txt", col=2)
    phase = read_data("prop.txt", col=13)
    for i, c in enumerate(cond):
        p = phase[i][2]
        if p in ["vapor"]:
            continue
        model.prop_in.temperature.set_value(c[0])
        h = value(model.prop_in.func_hlpt(c[1]/1000, 647.096/c[0]))
        rho = value(model.prop_in.func_delta_liq(c[1]/1000, 647.096/c[0]))*322
        print("{}, {}, {}, {}, {}, {}, {}".format(c[0], c[1], h, c[2], rho_dat[i][2], rho, p))
        assert(abs(h-c[2])/c[2]*100 < 1)

@pytest.mark.skipif(not iapws_available(), reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_entropy():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("prop.txt", col=6)
    phase = read_data("prop.txt", col=13)
    for i, c in enumerate(cond):
        if phase[i][2] in ["liquid", "supercritical"]:
            p = "Liq"
        else:
            p = "Vap"
        model.prop_in.temperature.set_value(c[0])
        model.prop_in.pressure = c[1]
        rho = value(model.prop_in.entr_mol_phase[p]/model.prop_in.mw/1000)
        print("{}, {}, {}, {}, {}".format(c[0], c[1], rho, c[2], p))
        assert(abs(rho-c[2])/c[2]*100 < 1)

@pytest.mark.skipif(not iapws_available(), reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_speed_of_sound():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("prop.txt", col=9)
    phase = read_data("prop.txt", col=13)
    for i, c in enumerate(cond):
        if c[2] == "undefined":
            continue
        if (c[0] > 640 and c[0] < 659) and (c[1] > 2.1e7 and c[1] < 2.5e7):
            #near critical and non-analytic terms were omitted
            continue
        if phase[i][2] in ["liquid", "supercritical"]:
            p = "Liq"
        else:
            p = "Vap"
        model.prop_in.temperature.set_value(c[0])
        model.prop_in.pressure = c[1]
        rho = value(model.prop_in.speed_sound_phase[p])
        print("{}, {}, {}, {}, {}".format(c[0], c[1], rho, c[2], p))
        assert(abs(rho-c[2])/c[2]*100 < 1)

@pytest.mark.skipif(not iapws_available(), reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_cp():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("prop.txt", col=8)
    phase = read_data("prop.txt", col=13)
    for i, c in enumerate(cond):
        if c[2] == "undefined":
            continue
        if (c[0] > 640 and c[0] < 680) and (c[1] > 2.1e7 and c[1] < 3.1e7):
            #near critical and non-analytic terms were omitted
            continue
        if phase[i][2] in ["liquid", "supercritical"]:
            p = "Liq"
        else:
            p = "Vap"
        model.prop_in.temperature.set_value(c[0])
        model.prop_in.pressure = c[1]
        rho = value(model.prop_in.cp_mol_phase[p]/model.prop_in.mw/1000)
        print("{}, {}, {}, {}, {}".format(c[0], c[1], rho, c[2], p))
        assert(abs(rho-c[2])/c[2]*100 < 1)

@pytest.mark.skipif(not iapws_available(), reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_cv():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("prop.txt", col=7)
    phase = read_data("prop.txt", col=13)
    for i, c in enumerate(cond):
        if c[2] == "undefined":
            continue
        if (c[0] > 640 and c[0] < 680) and (c[1] > 2.1e7 and c[1] < 3.1e7):
            #near critical and non-analytic terms were omitted
            continue
        if phase[i][2] in ["liquid", "supercritical"]:
            p = "Liq"
        else:
            p = "Vap"
        model.prop_in.temperature.set_value(c[0])
        model.prop_in.pressure = c[1]
        rho = value(model.prop_in.cv_mol_phase[p]/model.prop_in.mw/1000)
        print("{}, {}, {}, {}, {}".format(c[0], c[1], rho, c[2], p))
        assert(abs(rho-c[2])/c[2]*100 < 1)
