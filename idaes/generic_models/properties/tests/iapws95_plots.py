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

from pyomo.environ import *
from pyomo.core.base.external import AMPLExternalFunction
from pyomo.opt import SolverFactory
from idaes.generic_models.properties import iapws95
import numpy as np
import matplotlib.pyplot as plt

def make_model():
    model = ConcreteModel()
    model.prop_param = iapws95.Iapws95ParameterBlock()
    model.prop_in = iapws95.Iapws95StateBlock(default={"parameters": model.prop_param})
    return model

def plot_temperature_vapor_fraction(n=100, h=(0, 6000), p=10000):
    model = make_model()
    prop = model.prop_in
    hlist = [h[0] + i*(h[1] - h[0])/n for i in range(n)]
    vf = [0]*n
    T = [0]*n
    for i, h in enumerate(hlist):
        vf[i] = value(prop.func_vf(h, p))
        T[i] = 647.096/value(prop.func_tau(h, p))

    fig, ax1 = plt.subplots()
    ax1.set_ylim([-0.1, 1.1])
    ax1.set_xlabel('enthalpy (kJ/kg)')
    ax1.set_ylabel('vapor fraction', color=(0, 0.5, 0))
    ax1.plot(hlist, vf, color=(0, 0.5, 0))
    ax2 = ax1.twinx()
    ax2.set_ylabel('temperature (K)', color=(0, 0, 0.5))
    ax2.plot(hlist, T, color=(0, 0, 0.5))
    fig.tight_layout()
    plt.show()

def plot_psat(n=100, T=(200, 700)):
    model = make_model()
    prop = model.prop_in
    Tlist = [T[0] + i*(T[1] - T[0])/n for i in range(n)]
    p = [0]*n
    for i, T in enumerate(Tlist):
        p[i] = value(prop.func_p_sat(647.096/T))
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('temperature (K)')
    ax1.set_ylabel('saturation pressure (kPa)', color=(0, 0.5, 0))
    ax1.plot(Tlist, p, color=(0, 0.5, 0))
    fig.tight_layout()
    plt.show()

def plot_Tsat(n=100, p=(1, 25000)):
    model = make_model()
    prop = model.prop_in
    Plist = [p[0] + i*(p[1] - p[0])/n for i in range(n)]
    T = [0]*n
    for i, P in enumerate(Plist):
        T[i] = 647.096/value(prop.func_tau_sat(P))
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('pressure (kPa)')
    ax1.set_ylabel('saturation temperature (K)', color=(0, 0.5, 0))
    ax1.plot(Plist, T, color=(0, 0.5, 0))
    fig.tight_layout()
    plt.show()

def plot_hpt(n=100, T=(240, 1500), p=20000):
    model = make_model()
    prop = model.prop_in
    Tsat = 647.096/value(prop.func_tau_sat(p))
    print("Tsat = {}".format(Tsat))
    Tlist = [T[0] + i*(T[1] - T[0])/n for i in range(n)]
    h = [0]*n
    for i, T in enumerate(Tlist):
        if T > Tsat:
            f = prop.func_hvpt
        else:
            f = prop.func_hlpt
        h[i] = value(f(p, 647.096/T))
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('temperature (kPa)')
    ax1.set_ylabel('enthalpy (kJ/kg)', color=(0, 0.5, 0))
    ax1.plot(Tlist, h, color=(0, 0.5, 0))
    fig.tight_layout()
    plt.show()

if __name__ == '__main__':
    #p=22064
    p = 24233
    plot_temperature_vapor_fraction(p=p)
    plot_hpt(n=200, p=p)
    #plot_psat(n=200)
    #plot_psat(n=100, T=(647.09, 647.096))
    #plot_Tsat(n=200)
