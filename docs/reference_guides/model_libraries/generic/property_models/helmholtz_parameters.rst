Helmholtz EoS Parameter and Expression files
============================================

.. module:: idaes.models.properties.general_helmholtz.helmholtz_parameters

Helmholtz equations of state are difficult to implement in an equation oriented framework.
Switching from the native temperature and density state variables to a more useful set and solving
the phase equilibrium problem are difficult and may have multiple solutions, only one of which
is physically meaningful.  External functions allow for exact first and second derivatives to
be calculated and used by the solver, while still allowing procedural code to be used for the
solution. The external function also serve as decomposition making the models
more robust by splitting off difficult-to-solve sub-problems to be solved separately.

To accommodate addition of substances using various forms of Helmholtz free energy models,
Helmholtz free energy and transport properties are read into the external function library as
AMPL NL-files.  Basic parameters are read as a json file.

Generating the expression and parameter files needed is done via utilities provided by IDAES.
The parameter and expression files originate from a Python script and an optional json master
parameter file.

This section describes how to add or modify substances for the Helmholtz equation of state
library.

Class
-----

The WriteParameters class creates the expression and parameter files needed by the
external functions to define a substance.  Parameters and expression are defined in a
script and an optional main parameter file.  Either custom or predefined expressions
can be used. 

.. autoclass:: WriteParameters
    :members:

Master Parameter File
---------------------

While it is possible to the define equation of state and transport property expressions and
parameters entirely in a Python script we will focus here on the use of a main parameter
json file, which is read and converted into expressions and parameters that can be used
by the external functions.  The overall structure of the main json parameter file is
given below. The details of each section will be filled as we progress through this 
documentation::

    {
        "comp": "co2",
        "basic": {
            ...
        },
        "eos": {
            "reference": [
                "Span, R., and W. Wagner (1996). A New Equation of State for Carbon Dioxide",
                "    Covering the Fluid Region from the Triple-Point Temperature to 1100 K as",
                "    Pressures up to 800 MPa. Journal of Physical and Chemical Reference Data,",
                "    25, 1509."
            ],
            ...
        },
        "aux": {
            "reference": [
                "Span, R., and W. Wagner (1996). A New Equation of State for Carbon Dioxide",
                "    Covering the Fluid Region from the Triple-Point Temperature to 1100 K as",
                "    Pressures up to 800 MPa. Journal of Physical and Chemical Reference Data,",
                "    25, 1509."                                                              
            ],
            ...
        },
        "transport": {
            "thermal_conductivity": {
                "reference": [
                    "Vesovic, V., W.A. Wakeham, G.A. Olchowy, J.V. Sengers, J.T.R. Watson, J.",
                    "    Millat, (1990). The transport properties of carbon dioxide. J. Phys.",
                    "    Chem. Ref. Data, 19, 763-808."
                ]
                ...
            },
            "viscosity": {
                "reference": [
                    "Fenghour, A., W.A. Wakeham, V. Vesovic, (1998). The Viscosity of Carbon",
                    "     Dioxide. J. Phys. Chem. Ref. Data, 27, 31-44."
                ]
                ...
            },
            "surface_tension": {
                "reference": [
                    "Mulero, A., I. Cachadina, Parra, M., Recommended Correlations for the Surface ",
                    "    Tension of Common Fluids, J. Phys. Chem. Ref. Data 41, 043105, (2012)."
                ],
                ...
            }
        }
    }

By convention the main parameter file is named {comp}.json, where {comp} is the component name
(e.g. co2.json); although, the main parameter file name can be anything.

Script
------

An example of the most basic script for reading a main parameter file and generating
expression and parameter files is given below. Predefined transport expressions are 
not yet available, so the script will also commonly define transport property expressions, 
and that will be covered later::

    from idaes.models.properties.general_helmholtz.helmholtz_parameters import (
    WriteParameters,
    )


    def main():
        """Generate parameter and expression files."""
        we = WriteParameters(parameters="r227ea.json")
        we.write()


    if __name__ == "__main__":
        main()

Basic Parameters
----------------

The basic parameters required are listed in the table below.

=============== ======================================================================= ==========
Parameter       Description                                                             Units
=============== ======================================================================= ==========
``R``           specific gas constant (ideal gas constant/molecular weight)             kJ/kg/K
``MW``          molecular weight                                                        g/mol 
``T_star``      reducing temperature (usually critical temperature)                     K 
``rho_star``    reducing density (usually critical density)                             kg/m3
``Tc``          critical temperature                                                    K
``rhoc``        critical density                                                        kg/m3
``Pc``          critical pressure (recalculated to match Tc and rhoc)                   kPa
``Tt``          triple point temperature (can be approximate)                           K
``Pt``          triple point pressure (can be approximate)                              kPa
``rhot_l``      triple point liquid density (can be approximate)                        kg/m3
``rhot_v``      triple point vapor density (can be approximate)                         kg/m3
``P_min``       minimum pressure                                                        kPa
``P_max``       maximum pressure                                                        kPa
``rho_max``     maximum density                                                         kg/m3
``T_min``       minimum temperature                                                     K
``T_max``       maximum temperature                                                     K
=============== ======================================================================= ==========

A truncated example of the main json parameter file is given below to show the form for 
basic parameters::

    {
        "comp": "co2",
        "basic": {
            "R": 0.1889241,
            "MW": 44.0098,
            "T_star": 304.1282,
            "rho_star": 467.6,
            "Tc": 304.1282,
            "rhoc": 467.6,
            "Pc": 7377.3,
            "Tt": 216.592,
            "Pt": 517.95,
            "rhot_l": 1178.46,
            "rhot_v": 13.761,
            "P_min": 1e-9,
            "P_max": 5e5,
            "rho_max": 1600.0,
            "T_min": 216,
            "T_max": 1400
        },
        ...
    }

Reference State
---------------

Some properties such as enthalpy, entropy, and internal energy don't have a well defined
absolute value, and can only really be measured relative to a reference state. For the
purpose of the parameter files, equations of state are initially implemented using whatever
reference state is used in the original paper.  An offset can be set in the parameter files
to adjust the reference state if desired.

The ``calculate_reference_offset(delta, tau, s0, h0)`` method in the WriteParameters class 
can be used to calculate a new reference state.  The temperature, density, reference enthalpy,
and reference entropy are required.  The temperature and density can be calculated using the 
IDAES property functions possibly using the expression writer class. The reference enthalpy and
entropy are chosen arbitrarily. The reference state offset is dimensionless and is added to the
:math:`n^\circ_1` and :math:`n^\circ_2` parameters in the ideal part for the dimensionless Helmholtz free
energy correlation.  

The following example shows how to set the reference state in the main parameter file::

    {
        "comp": "co2",
        "basic": {
            ...
        },
        "eos": {
            ...
            "reference_state_offset": [
                -14.4979156224319,
                8.82013935801453 
            ],
            ...
        },
        ...
    }


Predefined Expressions
----------------------

Commonly used expressions are predefined and can be used by specifying expression types 
in the main parameter file.

Ideal Part of Dimensionless Helmholtz Free Energy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Type 1**

The predefined type 1 form of the ideal portion of the dimensionless Helmholtz free energy 
is shown below. Parameters are provided in the ``"eos"`` section of the main parameter file. 
This form of the expression requires a dictionary of the :math:`n^\circ_i$`` 
parameters as ``n0`` and the math `\gamma^\circ$` parameters as ``g0``.  The last term
index in the sum (:math:`h`) is provided as ``last_term_ideal``.  

.. math::

    \log_e \delta + n^\circ_1 + n^\circ_2 \tau + n^\circ_3 \log_e \tau + 
    \sum_{i = 4}^h n^\circ_i \log_e \left[ 1 - \exp(-\gamma^\circ_i  \tau)\right]

A truncated example of a main parameter file is provided below, which uses this the type 1
form of the dimensionless Helmholtz free energy::

    {
        "comp": "co2",
        "basic": {
            ...
        },
        "eos": {
            ...
            "n0": {
                "1": 8.37304456,
                "2": -3.70454304,
                "3": 2.5,
                "4": 1.99427042,
                "5": 0.62105248,
                "6": 0.41195293,
                "7": 1.04028922,
                "8": 0.08327678
            },
            "g0": {
                "4": 3.15163,
                "5": 6.11190,
                "6": 6.77708,
                "7": 11.32384,
                "8": 27.08792
            },
            "last_term_ideal": 8,
            "phi_ideal_type": 1,
            ...
        },
        ...
    }

**Type 2**

Parameters for the type 2 version are similar to type 1, but the ``last_term_ideal`` parameter is 
a list [h1, h2].

.. math::

    \log_e \delta + n^\circ_1 + n^\circ_2 \tau + n^\circ_3 \log_e \tau  + 
    \sum_{i = 4}^{h_1} n^\circ_i  \tau^{\gamma^\circ_i} + 
    \sum_{i = h_1 + 1}^{h_2} n^\circ_i \log_e \left[ 1 - \exp(-\gamma^\circ_i  \tau)\right]


**Type 3**

Parameters for the type 3 version are similar to type 1.

.. math::

    \log_e \delta + n^\circ_1 + n^\circ_2 \tau + n^\circ_3 \log_e \tau  + 
    \sum_{i = 4}^{h} n^\circ_i  \tau^{\gamma^\circ_i}

**Type 4**

Parameters for the type 4 version are similar to type 1

.. math::
     \log_e \delta + n^\circ_1 + n^\circ_2 \tau + (n^\circ_3 - 1) \log_e \tau + 
     \sum_{i = 4}^h n^\circ_i \log_e \left[ 1 - \exp(\frac{-\gamma^\circ_i  \tau}{T_c})\right]

     
Residual Part of Dimensionless Helmholtz Free Energy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Predefined residual term expressions are discussed here.  The ``last_term_residual`` 
parameter is a list [h1, h2, ...].  Equation parameters are given by dictionaries 
in the ``"eos"`` section.  Greek letter parameters are given by the parameters 
:math:`\alpha` = ``a``, :math:`\beta` = ``b``, :math:`\varepsilon` = ``e``, and 
:math:`\gamma` = ``g``.   

The example below shows how to use a predefined residual dimensionless Helmholtz 
free energy expression::

    {
        "comp": "co2",
        ...
        "eos": {
            ...
            "c": {
                "8": 1,
                ...
                "34": 6
            },
            "d": {
                "1": 1,
                ...
                "39": 3
            },
            "t": {
                "1": 0.00,
                ...
                "39": 3.00
            },
            "n": {
                "1": 3.88568232031610e-01,
                ...
                "42": 5.50686686128420e-02
            },
            "a": {
                "35": 25,
                ...
                "39": 20
            },
            "b": {
                "35": 325,
                ...
                "39": 275
            },
            "g": {
                "35": 1.16,
                ...
                "39": 1.22
            },
            "e": {
                "35": 1,
                "36": 1,
                "37": 1,
                "38": 1,
                "39": 1
            },
            "last_term_residual": [7, 34, 39],
            "phi_residual_type": 2
            ...
        }
        ...
    }

**Type 1**

.. math::

	\phi^r(\delta, \tau) = \sum^{h_1}_1 n_i \delta^{d_i} \tau^{t_i} + \sum_{h_1 + 1}^{h_2} n_i \delta^{d_i} \tau^{t_i} \exp(-\delta^{c_i})

**Type 2**

.. math::

	\phi^r(\delta, \tau) = \sum^{h_1}_{i=1} n_i \delta^{d_i} \tau^{t_1} + 
	\sum_{i =h_1 + 1}^{h_2} n_i \delta^{d_i} \tau^{t_1} \exp(-\delta^{c_i}) + 
	\sum_{i = h_2 + 1}^{h_3} n_i \delta^{d_i} \tau^{t_1} \exp\left[-\alpha_i(\delta - \varepsilon_i)^2 - \beta_i(\tau - \gamma_i)^2 \right]

**Type 3**

.. math::

	\phi^r(\delta, \tau) = 
	\sum^{h_1}_{i=1} n_i \delta^{d_i} \tau^{t_1} + 
	\sum_{i=h_1 + 1}^{h_2} n_i \delta^{d_i} \tau^{t_1} \exp(-\delta^{c_i})
	\sum_{i=h_2 + 1}^{h_3} n_i \delta^{d_i} \tau^{t_1} \exp(-\delta^{c_i}) \exp(-\tau^{b_i})

**Type 4**

.. math::

	\phi^r(\delta, \tau) = 
	\sum^{h_0}_{i=1} n_i \delta^{d_i} \tau^{t_1} + 
	\sum^m_{j=1} \left[ \exp(-\delta^j) \sum^{h_j}_{i = h_{j-1} + 1} n_i \delta^{d_i} \tau^{t_1} \right]


**Type 5**

.. math::
    \phi^r(\delta, \tau) = \sum^{h_1}_{i=1} n_i \delta^{d_i} \tau^{t_1} + 
	\sum_{i =h_1 + 1}^{h_2} n_i \delta^{d_i} \tau^{t_1} \exp(-\delta^{c_i}) + 
	\sum_{i = h_2 + 1}^{h_3} n_i \delta^{d_i} \tau^{t_1} \exp\left[-\alpha_i(\delta - \varepsilon_i)^2 - \beta_i(\tau - \gamma_i)^2 \right] +
    \sum_{i = h_3+1}^{h_4} n_i \delta^{d_i} \tau^{t_1} \exp\left[-\alpha_i(\delta - \varepsilon_i)^2 - \frac{1}{\beta_i(\tau - \gamma_i)^2 + b_i} \right]
    

Approximate Saturated Reduced Density
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To assist in the equilibrium calculation, curves for approximate saturated liquid 
and vapor density are required. These curves are typically provided in the papers 
where the equations of state are defined. There are two common forms for these 
curves, which we call type 1 and type 2 for convenience.

An example main parameter file entry for the approximate density curves is given 
below::

    {
        "comp": "co2",
        ...
        "aux": {
            "reference": [
                "Span, R., and W. Wagner (1996). A New Equation of State for Carbon Dioxide",
                "    Covering the Fluid Region from the Triple-Point Temperature to 1100 K as",
                "    Pressures up to 800 MPa. Journal of Physical and Chemical Reference Data,",
                "    25, 1509."                                                              
            ],
            "delta_l_sat_approx": {
                "c": 1,
                "n": {
                    "1": 1.9245108,
                    "2": -0.62385555,
                    "3": -0.32731127,
                    "4": 0.39245142
                },
                "t": {
                    "1": 0.34,
                    "2": 0.5,
                    "3": 1.6666666666666667,
                    "4": 1.8333333333333333
                },
                "type": 2
            },
            "delta_v_sat_approx": {
                "c": 1,
                "n": {
                    "1": -1.7074879,
                    "2": -0.82274670,
                    "3": -4.6008549,
                    "4": -10.111178,
                    "5": -29.742252
                },
                "t": {
                    "1": 0.340,
                    "2": 0.5,
                    "3": 1.0,
                    "4": 2.3333333333333333,
                    "5": 4.666666666666667
                },
                "type": 2
            }
        },
        ...
    }

*Type 1*

.. math:: 
 
    \delta = c + \sum_{i=1}^h n_i  \left(1 - \frac{T}{T_c} \right)^{t_i}
    
*Type 2*

.. math::

    \delta = c \exp \left[ \sum_{i=0}^h \left(1 - \frac{T}{T_c} \right)^{t_i} \right]

Surface Tension
~~~~~~~~~~~~~~~

The surface tension external functions provide surface tension in units of mN/m.  
There is only one form for surface tension and it is only a function of temperature.
The parameters required are ``Tc`` [K], ``s`` [mN/m] and ``n`` [dimensionless]. ``Tc`` 
is provided again in case the critical temperature differs slightly from the equation 
of state parameter. The number of terms is inferred from the ``s`` dictionary. This 
version of the surface tension expression maxes out at the critical surface tension, 
so as not to cause numerical issues at temperatures above critical. 

.. math::
    
    \sigma = \sum_{i = 0}^h s_i \max\left(1 - \frac{T}{T_c}, 0 \right)^{n_i}

An example surface tension section is given below::

    {
        "comp": "co2",
        "transport": {
            ...
            "surface_tension": {
                "reference": [
                    "Mulero, A., I. Cachadina, Parra, M., Recommended Correlations for the Surface ",
                    "    Tension of Common Fluids, J. Phys. Chem. Ref. Data 41, 043105, (2012)."
                ],
                "Tc": 304.1282,
                "s": {
                    "0": 78.63
                },
                "n": {
                    "0": 2.471
                },
                "type": 1
            }
        }
    }

Viscosity
~~~~~~~~~

Currently, there are no predefined viscosity expressions.

Thermal Conductivity
~~~~~~~~~~~~~~~~~~~~

Currently, there are no predefined viscosity expressions.

Custom Expressions
------------------

Although any custom expressions can be defined rather than using the predefined ones,
we will focus here on viscosity and thermal conductivity since predefined expressions
are not available from them yet.  The water example demonstrates all the concepts.

Variables
~~~~~~~~~

Custom expressions are written on Pyomo models. With the ``delta`` (reduced  density)
and ``tau`` (1/reduced temperature).  The easiest way to add an expression is to write a
rule for it then use the ``add()`` method to add the expression to the WriteParameters
class.

Units
~~~~~

Most expression used are dimensionless.  Only viscosity (microPa*s), thermal conductivity 
(mW/m/K) and surface tension (mN/m) have units.

Property Functions
~~~~~~~~~~~~~~~~~~

If expressions require calculation of other properties (e.g. pressure or heat capacity)
they can call other property functions, as demonstrated in the example up next.

Example
~~~~~~~

The example for water below shows how custom expressions can be defined. Note that heat
capacity, viscosity, and isothermal compressibility are used in the thermal conductivity 
expression by calling property functions::

    import math
    import pyomo.environ as pyo
    from pyomo.common.fileutils import find_library
    from idaes.core.util.math import smooth_max
    from idaes.models.properties.general_helmholtz.helmholtz_parameters import (
        WriteParameters,
    )

    def thermal_conductivity_rule(m):
        """Thermal conductivity rule

        International Association for the Properties of Water and Steam (2011).
            IAPWS R15-11, "Release on the IAPWS Formulation 2011 for the
            Thermal Conductivity of Ordinary Water Substance,"
            URL: http://iapws.org/relguide/ThCond.pdf
        """

        L0 = {
            0: 2.443221e-3,
            1: 1.323095e-2,
            2: 6.770357e-3,
            3: -3.454586e-3,
            4: 4.096266e-4,
        }

        L1 = {
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
        }
        delta = m.delta
        tau = m.tau
        big_lam = 177.8514
        qdb = 1.0 / 0.40
        nu = 0.630
        gamma = 1.239
        xi0 = 0.13
        big_gam0 = 0.06
        Tbr = 1.5
        Tb = 1 / tau
        lib_path = find_library("general_helmholtz_external")
        m.cp = pyo.ExternalFunction(library=lib_path, function="cp")
        m.cv = pyo.ExternalFunction(library=lib_path, function="cv")
        m.mu = pyo.ExternalFunction(library=lib_path, function="mu")
        m.itc = pyo.ExternalFunction(
            library=lib_path, function="itc"
        )  # isothermal compressibility [1/MPa]
        deltchi = smooth_max(
            delta**2
            * (m.itc("h2o", delta, tau) - m.itc("h2o", delta, 1 / Tbr) * Tbr / Tb)
            * m.Pc
            / 1000,
            0,
            eps=1e-8,
        )
        xi = xi0 * (deltchi / big_gam0) ** (nu / gamma)
        y = qdb * xi
        kappa = m.cp("h2o", delta, tau) / m.cv("h2o", delta, tau)
        Z = (
            2.0
            / math.pi
            / y
            * (
                ((1 - 1.0 / kappa) * pyo.atan(y) + 1.0 / kappa * y)
                - (1 - pyo.exp(-1 / (1.0 / y + y**2 / 3 / delta**2)))
            )
        )
        cpb = m.cp("h2o", delta, tau) / m.R
        mub = m.mu("h2o", delta, tau)

        lam2 = big_lam * delta * Tb * cpb / mub * Z
        return lam2 + pyo.sqrt(1.0 / m.tau) / sum(
            L0val * m.tau**i for i, L0val in L0.items()
        ) * pyo.exp(
            m.delta
            * sum(
                (m.tau - 1) ** i * sum(L1[i, j] * (m.delta - 1) ** j for j in range(0, 6))
                for i in range(0, 5)
            )
        )


    def viscosity_rule(mvisc):
        """Viscosity calculation for water

        International Association for the Properties of Water and Steam (2008).
            IAPWS R12-08, "Release on the IAPWS Formulation 2008 for the Viscosity
            of Ordinary Water Substance, URL: http://iapws.org/relguide/visc.pdf
        """

        H0 = {0: 1.67752, 1: 2.20462, 2: 0.6366564, 3: -0.241605}
        H1 = {
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
        }

        return (
            1e2
            * pyo.sqrt(1.0 / mvisc.tau)
            / sum(H0val * mvisc.tau**i for i, H0val in H0.items())
            * pyo.exp(
                mvisc.delta
                * sum(
                    (mvisc.tau - 1) ** i
                    * sum(H1[i, j] * (mvisc.delta - 1) ** j for j in range(0, 7))
                    for i in range(0, 6)
                )
            )
        )


    def main():
        """Generate property and expression files."""

        we = WriteParameters(parameters="h2o.json")
        we.add(
            {
                "viscosity": viscosity_rule,
                "thermal_conductivity": thermal_conductivity_rule,
            }
        )
        we.write()


    if __name__ == "__main__":
        main()


