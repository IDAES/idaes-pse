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
"""IDAES IAPWS-95 Steam properties

Dropped all critical enhancments and non-analytic terms ment to improve
accruacy near the critical point. These tend to cause singularities in the
equations, and it is assumend that we will try to avoid operating very close to
the critical point.

 References: (some of this is only used in the C++ part)
   International Association for the Properties of Water and Steam (2016).
       IAPWS R6-95 (2016), "Revised Release on the IAPWS Formulation 1995 for
       the Properties of Ordinary Water Substance for General Scientific Use,"
       URL: http://iapws.org/relguide/IAPWS95-2016.pdf
   Wagner, W.,  A. Pruss (2002). "The IAPWS Formulation 1995 for the
       Thermodynamic Properties of Ordinary Water Substance for General and
       Scientific Use." J. Phys. Chem. Ref. Data, 31, 387-535.
   Wagner, W. et al. (2000). "The IAPWS Industrial Formulation 1997 for the
       Thermodynamic Properties of Water and Steam," ASME J. Eng. Gas Turbines
       and Power, 122, 150-182.
   Akasaka, R. (2008). "A Reliable and Useful Method to Determine the
       Saturation State from Helmholtz Energy Equations of State." Journal of
       Thermal Science and Technology, 3(3), 442-451.
   International Association for the Properties of Water and Steam (2011).
       IAPWS R15-11, "Release on the IAPWS Formulation 2011 for the
       Thermal Conductivity of Ordinary Water Substance,"
       URL: http://iapws.org/relguide/ThCond.pdf
   International Association for the Properties of Water and Steam (2008).
       IAPWS R12-08, "Release on the IAPWS Formulation 2008 for the Viscosity
       of Ordinary Water Substance,"
       URL: http://iapws.org/relguide/visc.pdf
"""
__author__ = "John Eslick"

# Import Python libraries
import os

# Import Pyomo libraries
from pyomo.environ import (
    Expression,
    Param,
    RangeSet,
    Set,
    exp,
    sqrt,
    units as pyunits
)

# Import IDAES
from idaes.core import declare_process_block_class

import idaes
import idaes.logger as idaeslog
from idaes.generic_models.properties.helmholtz.helmholtz import (
    _available,
    _htpx,
    HelmholtzParameterBlockData,
    HelmholtzStateBlockData,
    HelmholtzThermoExpressions,
    PhaseType,
    StateVars,
    _StateBlock,
)

# Logger
_log = idaeslog.getLogger(__name__)
_so = os.path.join(idaes.bin_directory, "iapws95_external.so")


def iapws95_available():
    """Make sure the compiled IAPWS-95 functions are available. Yes, in Windows
    the .so extention is still used.
    """
    return _available(_so)


def htpx(T=None, P=None, x=None):
    """
    Convenience function to calculate steam enthalpy from temperature and
    either pressure or vapor fraction. This function can be used for inlet
    streams and initialization where temperature is known instead of enthalpy.
    User must provided values for two of T, P, or x.

    Args:
        T: Temperature with units (between 200 and 3000 K)
        P: Pressure with units (between 1 and 1e9 Pa), None if saturated vapor
        x: Vapor fraction [mol vapor/mol total] (between 0 and 1), None if
        superheated or subcooled

    Returns:
        Total molar enthalpy [J/mol].
    """
    prop = Iapws95StateBlock(default={"parameters": Iapws95ParameterBlock()})
    return _htpx(T=T, P=P, x=x, prop=prop,
                 Tmin=270*pyunits.K, Tmax=3e3*pyunits.K,
                 Pmin=1e-4*pyunits.kPa, Pmax=1e6*pyunits.kPa)


@declare_process_block_class("Iapws95ParameterBlock")
class Iapws95ParameterBlockData(HelmholtzParameterBlockData):
    CONFIG = HelmholtzParameterBlockData.CONFIG()

    def build(self):
        self._set_parameters(
            library=_so,
            eos_tag="iapws95",
            state_block_class=Iapws95StateBlock,
            component_list=Set(initialize=["H2O"]),
            phase_equilibrium_idx=Set(initialize=[1]),
            phase_equilibrium_list={1: ["H2O", ("Vap", "Liq")]},
            mw=Param(initialize=0.01801528,
                     doc="Molecular weight [kg/mol]",
                     units=pyunits.kg/pyunits.mol),
            temperature_crit=Param(
                initialize=647.096,
                doc="Critical temperature [K]",
                units=pyunits.K),
            pressure_crit=Param(initialize=2.2064e7,
                                doc="Critical pressure [Pa]",
                                units=pyunits.Pa),
            dens_mass_crit=Param(initialize=322,
                                 doc="Critical density [kg/m3]",
                                 units=pyunits.kg/pyunits.m**3),
            specific_gas_constant=Param(
                initialize=461.51805,
                doc="Water Specific Gas Constant [J/kg/K]",
                units=pyunits.J/pyunits.kg/pyunits.K),
            pressure_bounds=(0.1, 1e9),
            temperature_bounds=(250, 2500),
            enthalpy_bounds=(0, 1e5),
        )
        super().build()
        # Thermal conductivity parameters.
        # "Release on the IAPWS Formulation 2011 for the Thermal Conductivity
        # of Ordinary Water Substance"
        self.tc_L0 = Param(
            RangeSet(0, 5),
            initialize={
                0: 2.443221e0,
                1: 1.323095e1,
                2: 6.770357e0,
                3: -3.454586e0,
                4: 4.096266e-1,
            },
            doc="0th order thermal conductivity parameters",
            units=pyunits.K*pyunits.m/pyunits.W
        )

        self.tc_L1 = Param(
            RangeSet(0, 5),
            RangeSet(0, 6),
            initialize={
                (0, 0): 1.60397357,
                (1, 0): 2.33771842,
                (2, 0): 2.19650529,
                (3, 0): -1.21051378,
                (4, 0): -2.7203370,
                (0, 1): -0.646013523,
                (1, 1): -2.78843778,
                (2, 1): -4.54580785,
                (3, 1): 1.60812989,
                (4, 1): 4.57586331,
                (0, 2): 0.111443906,
                (1, 2): 1.53616167,
                (2, 2): 3.55777244,
                (3, 2): -0.621178141,
                (4, 2): -3.18369245,
                (0, 3): 0.102997357,
                (1, 3): -0.463045512,
                (2, 3): -1.40944978,
                (3, 3): 0.0716373224,
                (4, 3): 1.1168348,
                (0, 4): -0.0504123634,
                (1, 4): 0.0832827019,
                (2, 4): 0.275418278,
                (3, 4): 0.0,
                (4, 4): -0.19268305,
                (0, 5): 0.00609859258,
                (1, 5): -0.00719201245,
                (2, 5): -0.0205938816,
                (3, 5): 0.0,
                (4, 5): 0.012913842,
            },
            doc="1st order thermal conductivity parameters",
        )
        # Viscosity parameters
        # "Release on the IAPWS Formulation 2008 for the Viscosity of
        # Ordinary Water Substance "
        self.visc_H0 = Param(
            RangeSet(0, 4),
            initialize={0: 1.67752, 1: 2.20462, 2: 0.6366564, 3: -0.241605},
            doc="0th order viscosity parameters",
            units=1/pyunits.s/pyunits.Pa
        )

        self.visc_H1 = Param(
            RangeSet(0, 6),
            RangeSet(0, 7),
            initialize={
                (0, 0): 5.20094e-1,
                (1, 0): 8.50895e-2,
                (2, 0): -1.08374,
                (3, 0): -2.89555e-1,
                (4, 0): 0.0,
                (5, 0): 0.0,
                (0, 1): 2.22531e-1,
                (1, 1): 9.99115e-1,
                (2, 1): 1.88797,
                (3, 1): 1.26613,
                (4, 1): 0.0,
                (5, 1): 1.20573e-1,
                (0, 2): -2.81378e-1,
                (1, 2): -9.06851e-1,
                (2, 2): -7.72479e-1,
                (3, 2): -4.89837e-1,
                (4, 2): -2.57040e-1,
                (5, 2): 0.0,
                (0, 3): 1.61913e-1,
                (1, 3): 2.57399e-1,
                (2, 3): 0.0,
                (3, 3): 0.0,
                (4, 3): 0.0,
                (5, 3): 0.0,
                (0, 4): -3.25372e-2,
                (1, 4): 0.0,
                (2, 4): 0.0,
                (3, 4): 6.98452e-2,
                (4, 4): 0.0,
                (5, 4): 0.0,
                (0, 5): 0.0,
                (1, 5): 0.0,
                (2, 5): 0.0,
                (3, 5): 0.0,
                (4, 5): 8.72102e-3,
                (5, 5): 0.0,
                (0, 6): 0.0,
                (1, 6): 0.0,
                (2, 6): 0.0,
                (3, 6): -4.35673e-3,
                (4, 6): 0.0,
                (5, 6): -5.93264e-4,
            },
            doc="1st order viscosity parameters",
        )
        self.set_default_scaling("therm_cond_phase", 1e2, index="Liq")
        self.set_default_scaling("therm_cond_phase", 1e1, index="Vap")
        self.set_default_scaling("visc_d_phase", 1e5, index="Liq")
        self.set_default_scaling("visc_d_phase", 1e6, index="Vap")
        self.set_default_scaling("visc_k_phase", 1e5, index="Liq")
        self.set_default_scaling("visc_k_phase", 1e7, index="Vap")


@declare_process_block_class(
    "Iapws95StateBlock",
    block_class=_StateBlock,
    doc="""This is some placeholder doc.""",
)
class Iapws95StateBlockData(HelmholtzStateBlockData):
    """
    This is a property package for calculating thermophysical properties of
    water.
    """

    def build(self, *args):
        """
        Callable method for Block construction
        """
        super().build(*args)

        phlist = self.config.parameters.private_phase_list
        tau = self.tau
        delta = self.dens_phase_red  # this shorter name is from IAPWS

        # Phase Thermal conductiviy
        def rule_tc(b, p):
            L0 = self.config.parameters.tc_L0
            L1 = self.config.parameters.tc_L1
            return (
                sqrt(1.0 / tau) / sum(L0[i] * tau ** i for i in L0) *
                exp(delta[p] * sum((tau - 1)**i *
                    sum(L1[i, j] * (delta[p] - 1)**j for j in range(0, 6))
                for i in range(0, 5))))

        self.therm_cond_phase = Expression(
            phlist, rule=rule_tc, doc="Thermal conductivity [W/K/m]"
        )

        # Phase dynamic viscosity
        def rule_mu(b, p):
            H0 = self.config.parameters.visc_H0
            H1 = self.config.parameters.visc_H1
            # The units of this are really weird, so I am just going to append
            # units to the expression rather than give units to the parameters
            return (
                1e-4 * sqrt(1.0 / tau) / sum(H0[i] * tau ** i for i in H0) *
                exp(delta[p] * sum((tau - 1)**i *
                    sum(H1[i, j] * (delta[p] - 1)**j for j in range(0, 7))
                for i in range(0, 6))))

        self.visc_d_phase = Expression(
            phlist, rule=rule_mu, doc="Viscosity (dynamic)"
        )

        # Phase kinimatic viscosity
        def rule_nu(b, p):
            return self.visc_d_phase[p] / self.dens_mass_phase[p]

        self.visc_k_phase = Expression(
            phlist, rule=rule_nu, doc="Kinematic viscosity"
        )
