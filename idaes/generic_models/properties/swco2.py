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
"""IDAES Span-Wager CO2 properties

Dropped all critical enhancments and non-analytic terms ment to improve
accruacy near the critical point. These tend to cause singularities in the
equations, and it is assumend that we will try to avoid operating very close to
the critical point.

 References: (some of this is only used in the C++ part)
   Span, R., W. Wagner, W. (1996). "A New Equation of State for Carbon Dioxide
       Covering the Fluid Region from the Triple-Point Temperature to 1100 K at
       Pressures up to 800 MPa." J. Phys. Chem. Ref. Data, 25(6), 1509-1596.
   Akasaka, R. (2008). "A Reliable and Useful Method to Determine the
       Saturation State from Helmholtz Energy Equations of State." Journal of
       Thermal Science and Technology, 3(3), 442-451.
   Vesovic, V., W.A. Wakeham, G.A. Olchowy, J.V. Sengers, J.T.R. Watson, J.
       Millat, (1990). "The transport properties of carbon dioxide." J. Phys.
       Chem. Ref. Data, 19, 763-808.
   Fenghour, A., W.A. Wakeham, V. Vesovic, (1998). "The Viscosity of Carbon
       Dioxide." J. Phys. Chem. Ref. Data, 27, 31-44.
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
    log,
)

# Import IDAES
from idaes.core import (
    declare_process_block_class,
    StateBlock,
    StateBlockData,
)
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
_so = os.path.join(idaes.bin_directory, "swco2_external.so")


def swco2_available():
    """Make sure the compiled IAPWS-95 functions are available. Yes, in Windows
    the .so extention is still used.
    """
    return _available(_so)


def htpx(T, P=None, x=None):
    """
    Convenience function to calculate enthalpy from temperature and
    either pressure or vapor fraction. This function can be used for inlet
    streams and initialization where temperature is known instead of enthalpy.

    User must provided values for one (and only one) of arguments P and x.

    Args:
        T: Temperature [K] (between 200 and 3000)
        P: Pressure [Pa] (between 1 and 1e9), None if saturated vapor
        x: Vapor fraction [mol vapor/mol total] (between 0 and 1), None if
        superheated or subcooled

    Returns:
        Total molar enthalpy [J/mol].
    """
    prop = SWCO2StateBlock(default={"parameters": SWCO2ParameterBlock()})
    return _htpx(T=T, P=P, x=x, prop=prop, Tmin=200, Tmax=3e3, Pmin=0.1, Pmax=1e9)


@declare_process_block_class("SWCO2ParameterBlock")
class SWCO2ParameterBlockData(HelmholtzParameterBlockData):
    CONFIG = HelmholtzParameterBlockData.CONFIG()

    def build(self):
        self._set_parameters(
            library=_so,
            eos_tag="swco2",
            state_block_class=SWCO2StateBlock,
            component_list=Set(initialize=["CO2"]),
            phase_equilibrium_idx=Set(initialize=[1]),
            phase_equilibrium_list={1: ["CO2", ("Vap", "Liq")]},
            mw=Param(initialize=0.0440098, doc="Molecular weight [kg/mol]"),
            temperature_crit=Param(
                initialize=304.1282,
                doc="Critical temperature [K]",
            ),
            pressure_crit=Param(initialize=7.377e6, doc="Critical pressure [Pa]"),
            dens_mass_crit=Param(initialize=467.6, doc="Critical density [kg/m3]"),
            specific_gas_constant=Param(
                initialize=188.9241,
                doc="CO2 Specific Gas Constant [J/kg/K]",
            ),
            pressure_bounds=(0.1, 1e9),
            temperature_bounds=(210, 2500),
            enthalpy_bounds=(-2e4, 1e5),
        )

        super().build()

        # Thermal conductivity parameters.
        # Vesovic et al. (1990)
        self.tc_b = Param(
            RangeSet(0, 8),
            initialize={
                0: 0.4226159,
                1: 0.6280115,
                2: -0.5387661,
                3: 0.6735941,
                4: 0,
                5: 0,
                6: -0.4362677,
                7: 0.2255388,
            },
            doc="0 density limit thermal conductivity parameters",
        )

        self.tc_c = Param(
            RangeSet(1, 6),
            initialize={
                1: 2.387869e-2,
                2: 4.350794,
                3: -10.33404,
                4: 7.981590,
                5: -1.940558,
            },
            doc="0 density limit thermal conductivity parameters",
        )

        self.tc_d = Param(
            RangeSet(1, 5),
            initialize={
                1: 2.447164e-2,
                2: 8.705605e-5,
                3: -6.547950e-8,
                4: 6.594919e-11,
            },
            doc="residual thermal conductivity parameters",
        )

        # Viscosity parameters
        # "Fenghour et al. (1998) with critial enhancment Vesovic et al. (1990)
        self.visc_a = Param(
            RangeSet(0, 5),
            initialize={
                0: 0.235156,
                1: -0.491266,
                2: 5.211155e-2,
                3: 5.347906e-2,
                4: -1.537102e-2,
            },
            doc="0 density limit viscosity parameters",
        )

        # The indexing looks a little weird here, but it's from the source
        self.visc_d = Param(
            {(1, 1), (2, 1), (6, 4), (8, 1), (8, 2)},
            initialize={
                (1, 1): 0.4071119e-2,
                (2, 1): 0.7198037e-4,
                (6, 4): 0.2411697e-16,
                (8, 1): 0.2971072e-22,
                (8, 2): -0.1627888e-22,
            },
            doc="residual viscosity parameters",
        )


@declare_process_block_class(
    "SWCO2StateBlock",
    block_class=_StateBlock,
    doc="""This is some placeholder doc.""",
)
class SWCO2StateBlockData(HelmholtzStateBlockData):
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
        T = self.temperature
        rho = self.dens_mass_phase

        # Phase Thermal conductiviy
        def rule_tc(b, p):
            if p == "Liq":
                self.scaling_factor[self.therm_cond_phase[p]] = 1e1
            else:
                self.scaling_factor[self.therm_cond_phase[p]] = 1e2
            b = self.config.parameters.tc_b
            c = self.config.parameters.tc_c
            d = self.config.parameters.tc_d
            cint_over_k = (
                1.0 + exp(-183.5/T)*sum(c[i]*(T/100)**(2-i) for i in c)
            )
            Ts = T/251.196
            G = sum(b[i]/Ts**i for i in b)
            return 1e-3*(
                0.475598*sqrt(T)*(1 + 2.0/5.0*cint_over_k)/G +
                sum(d[i]*rho[p]**i for i in d)
            )

        self.therm_cond_phase = Expression(
            phlist,
            rule=rule_tc,
            doc="Thermal conductivity [W/K/m]",
        )

        # Phase dynamic viscosity
        def rule_mu(b, p):
            if p == "Liq":
                self.scaling_factor[self.visc_d_phase[p]] = 1e5
            else:
                self.scaling_factor[self.visc_d_phase[p]] = 1e6
            a = self.config.parameters.visc_a
            d = self.config.parameters.visc_d
            Ts = T/251.196
            return 1e-6*(
                1.00697*sqrt(T)/exp(sum(a[i]*log(Ts)**i for i in a)) +
                d[1,1]*rho[p] +
                d[2,1]*rho[p]**2 +
                d[6,4]*rho[p]**6/Ts**3 +
                d[8,1]*rho[p]**8 +
                d[8,2]*rho[p]**8/Ts
            )

        self.visc_d_phase = Expression(
            phlist, rule=rule_mu, doc="Viscosity (dynamic) [Pa*s]"
        )

        # Phase kinimatic viscosity
        def rule_nu(b, p):
            if p == "Liq":
                self.scaling_factor[self.visc_k_phase[p]] = 1e5
            else:
                self.scaling_factor[self.visc_k_phase[p]] = 1e7
            return self.visc_d_phase[p] / self.dens_mass_phase[p]

        self.visc_k_phase = Expression(
            phlist, rule=rule_nu, doc="Kinematic viscosity [m^2/s]"
        )
