##################################################################################
#                                                                                #
# CO2 EOS Expressions and Parameters:                                            #
#                                                                                #
# Span, R., and W. Wanger (1996). "A New Equation of State for Carbon Dioxide    #
#     Covering the Fluid Region from the Triple-Point Temperature to 1100 K as   #
#     Pressures up to 800 MPa." Journal of Physical and Chemical Reference Data, #
#     25, 1509.                                                                  #
#                                                                                #
# Vesovic, V., W.A. Wakeham, G.A. Olchowy, J.V. Sengers, J.T.R. Watson, J.       #
#     Millat, (1990). "The transport properties of carbon dioxide." J. Phys.     #
#     Chem. Ref. Data, 19, 763-808.                                              #
#                                                                                #
# Fenghour, A., W.A. Wakeham, V. Vesovic, (1998). "The Viscosity of Carbon       #
#     Dioxide." J. Phys. Chem. Ref. Data, 27, 31-44.                             #
#                                                                                #
# Mulero, A., I. Cachadina, Parra, M., "Recommended Correlations for the         #
#     Surface Tension of Common Fluids," J. Phys. Chem. Ref. Data 41, 043105     #
#     (2012).                                                                    #
#                                                                                #
##################################################################################

import pyomo.environ as pyo
from idaes.core.util.math import smooth_max
from idaes.models.properties.general_helmholtz.helmholtz_parameters import (
    WriteParameters,
)


def thermal_conductivity_rule(m):
    b = {
        0: 0.4226159,
        1: 0.6280115,
        2: -0.5387661,
        3: 0.6735941,
        4: 0,
        5: 0,
        6: -0.4362677,
        7: 0.2255388,
    }
    c = {
        1: 2.387869e-2,
        2: 4.350794,
        3: -10.33404,
        4: 7.981590,
        5: -1.940558,
    }
    d = {
        1: 24.47164,
        2: 8.705605e-2,
        3: -6.547950e-5,
        4: 6.594919e-8,
    }
    T = m.T_star / m.tau
    Ts = T / 251.196
    rho = m.rho_star * m.delta
    G = sum(b[i] / Ts**i for i in b)
    cint_over_k = 1.0 + pyo.exp(-183.5 / T) * sum(
        c[i] * (T / 100) ** (2 - i) for i in c
    )
    return (
        475.598 * pyo.sqrt(T) * (1 + 2.0 / 5.0 * cint_over_k) / G
        + d[1] * rho
        + d[2] * rho**2
        + d[3] * rho**3
        + d[4] * rho**4
    ) / 1e3


def viscosity_rule(m):
    a = {
        0: 0.235156,
        1: -0.491266,
        2: 5.211155e-2,
        3: 5.347906e-2,
        4: -1.537102e-2,
    }
    d = {
        1: 0.4071119e-2,
        2: 0.7198037e-4,
        3: 0.2411697e-16,
        4: 0.2971072e-22,
        5: -0.1627888e-22,
    }
    T = m.T_star / m.tau
    rho = m.delta * m.rho_star
    Ts = T / 251.196
    return (
        1.00697 * pyo.sqrt(T) / pyo.exp(sum(a[i] * pyo.log(Ts) ** i for i in a))
        + d[1] * rho
        + d[2] * rho**2
        + d[3] * rho**6 / Ts**3
        + d[4] * rho**8
        + d[5] * rho**8 / Ts
    )


def surface_tension_rule(m):
    x = smooth_max(1 - 1 / m.tau * m.T_star / m.Tc, 0, 1e-8)
    sigma = {0: 0.07863}
    n = {0: 2.471}
    return sum(sigma[i] * x ** n[i] for i in n) * 1000


def main():

    we = WriteParameters(parameters="co2.json")
    we.add(
        {
            "viscosity": viscosity_rule,
            "thermal_conductivity": thermal_conductivity_rule,
        }
    )
    we.write()


if __name__ == "__main__":
    main()
