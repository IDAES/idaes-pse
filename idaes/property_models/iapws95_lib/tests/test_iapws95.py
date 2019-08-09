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

import pytest
from pyomo.environ import *
from pyomo.opt import SolverFactory
from idaes.property_models import iapws95
import csv
import math
import os

# Set module level pyest marker
pytestmark = pytest.mark.iapws
prop_available = iapws95.iapws95_available()


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


@pytest.mark.slow
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
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


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_liquid_density_sat():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("sat_prop.txt", col=2)
    for c in cond:
        if c[0] > 645: # getting very close to ciritical point
            tol = 0.05
        else:
            tol = 0.001
        model.prop_in.temperature.set_value(c[0])
        model.prop_in.pressure = c[1]
        rho = value(model.prop_in.dens_mass_phase["Liq"])
        #print("{:.2f}, {:.2f}, {:.2f}, {:.2f}".format(c[0], c[1], rho, c[2]))
        assert(abs(rho-c[2])/c[2] < tol)


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_vapor_density_sat():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("sat_prop.txt", col=14)
    for c in cond:
        if c[0] > 645: # getting very close to ciritical point
            tol = 0.01
        else:
            tol = 0.001
        model.prop_in.temperature.set_value(c[0])
        model.prop_in.pressure = c[1]
        rho = value(model.prop_in.dens_mass_phase["Vap"])
        assert(abs(rho-c[2])/c[2] < tol)


@pytest.mark.slow
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_liquid_enthalpy_sat():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("sat_prop.txt", col=5)
    for c in cond:
        if c[0] > 645: # getting very close to ciritical point
            tol = 0.01
        else:
            tol = 0.001
        model.prop_in.pressure = c[1]
        enth = value(model.prop_in.enth_mol_sat_phase["Liq"]/model.prop_in.mw/1000.0)
        assert(abs((enth-c[2])/c[2]) < tol)


@pytest.mark.slow
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_vapor_enthalpy_sat():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("sat_prop.txt", col=17)
    for c in cond:
        if c[0] > 645: # getting very close to ciritical point
            tol = 0.01
        else:
            tol = 0.001
        model.prop_in.pressure = c[1]
        enth = value(model.prop_in.enth_mol_sat_phase["Vap"]/model.prop_in.mw/1000.0)
        assert(abs((enth-c[2])/c[2]) < tol)


@pytest.mark.slow
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_enthalpy_of_vaporization():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond_liq = read_data("sat_prop.txt", col=5)
    cond_vap = read_data("sat_prop.txt", col=17)
    for i, c in enumerate(cond_liq):
        if c[0] > 645: # getting very close to ciritical point
            tol = 0.05
        else:
            tol = 0.001
        model.prop_in.pressure.value = c[1]
        enth = value(model.prop_in.dh_vap_mol/model.prop_in.mw/1000.0)
        enth_dat = cond_vap[i][2] - c[2]
        if abs(enth_dat) > 1e-8:
            assert(abs((enth-enth_dat)/enth_dat) < tol)
        else:
            assert(abs(enth-enth_dat) < tol)
    #Over Critical Pressure
    model.prop_in.pressure = model.prop_in.config.parameters.pressure_crit*1.1
    enth = value(model.prop_in.dh_vap_mol/model.prop_in.mw/1000.0)
    assert(abs(enth) < 0.001)


@pytest.mark.slow
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
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
        if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 645:
            tol = 0.03 # steep part in sc region
        elif c[1] < 20:
            tol = 0.005 #very low pressure < 20 Pa
        else:
            tol = 0.001
        assert(abs(rho-c[2])/c[2] < tol)


@pytest.mark.slow
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
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
        rho = value(model.prop_in.dens_mass_phase[p])
        if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
            tol = 0.03 # steep part in sc region
        elif c[1] < 20:
            tol = 0.005 #very low pressure < 20 Pa
        else:
            tol = 0.0015
        assert(abs(h-c[2])/c[2] < tol)


@pytest.mark.slow
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
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
        rho = value(model.prop_in.dens_mass_phase["Vap"])
        if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
            tol = 0.03 # steep part in sc region
        else:
            tol = 0.003
        assert(abs(h-c[2])/c[2] < tol)


@pytest.mark.slow
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
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
        rho = value(model.prop_in.dens_mass_phase["Liq"])
        if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
            tol = 0.03 # steep part in sc region
        else:
            tol = 0.006
        assert(abs(h-c[2])/c[2] < tol)


@pytest.mark.slow
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
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
        s = value(model.prop_in.entr_mol_phase[p]/model.prop_in.mw/1000)
        rho = value(model.prop_in.dens_mass_phase[p])
        if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
            tol = 0.03 # steep part in sc region
        else:
            tol = 0.003
        assert(abs(s-c[2])/c[2] < tol)


@pytest.mark.slow
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.nocircleci()
def test_internal_energy():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters":model.prop_param})
    cond = read_data("prop.txt", col=4)
    phase = read_data("prop.txt", col=13)
    for i, c in enumerate(cond):
        if phase[i][2] in ["liquid", "supercritical"]:
            p = "Liq"
        else:
            p = "Vap"
        model.prop_in.temperature.set_value(c[0])
        model.prop_in.pressure = c[1]
        u = value(model.prop_in.energy_internal_mol_phase[p]/model.prop_in.mw/1000)
        rho = value(model.prop_in.dens_mass_phase[p])
        if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
            tol = 0.02 # steep part in sc region
        else:
            tol = 0.002
        assert(abs(u-c[2])/c[2] < tol)


@pytest.mark.slow
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
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
        w = value(model.prop_in.speed_sound_phase[p])
        rho = value(model.prop_in.dens_mass_phase[p])
        if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
            tol = 0.03 # steep part in sc region
        else:
            tol = 0.005
        assert(abs(w-c[2])/c[2] < tol)


@pytest.mark.slow
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
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
        cp = value(model.prop_in.cp_mol_phase[p]/model.prop_in.mw/1000)
        rho = value(model.prop_in.dens_mass_phase[p])
        if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
            tol = 0.03 # steep part in sc region
        else:
            tol = 0.005
        assert(abs(cp-c[2])/c[2] < tol)


@pytest.mark.slow
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
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
        cv = value(model.prop_in.cv_mol_phase[p]/model.prop_in.mw/1000)
        rho = value(model.prop_in.dens_mass_phase[p])
        if rho > 250 and rho < 420 and c[0] < 700 and c[0] > 640:
            tol = 0.03 # steep part in sc region
        else:
            tol = 0.003
        assert(abs(cv-c[2])/c[2] < tol)
