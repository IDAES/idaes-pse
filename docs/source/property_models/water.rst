Water/Steam
===========

Two property modules are available for pure water and steam properties.  The
property modules use the same calculations and yield consistent results, but one
uses pressure and molar enthalpy as state variables and the other uses
temperature, pressure, and vapor fraction as the state variables.

.. toctree::
    :maxdepth: 1

    iapws_ph
    iapws_tpx

Units
-----

The water property modules are in SI units (m, kg, s, mol).

Methods
-------

These methods use the IAPWS-95 formulation for scientific use for thermodynamic
properties (:ref:`Wagner and Pruss, 2002 <wagner-2002>`; :ref:`IAPWS, 2016
<iapws-2016>`). To solve the phase equilibrium the method of :ref:`Akasaka
(2008) <akasaka-2008>` was used. For solving these equations some relations from
the IAPWS-97 formulation for industrial use are used as initial values
(:ref:`Wagner et al., 2002 <wagner-2002>`). The industrial formulation is
slightly discontinuous between different regions, so it may not be suitable for
optimization. In addition to thermodynamic quantities, viscosity and thermal
conductivity are calculated (:ref:`IAPWS, 2008 <iapws-2008>`;
:ref:`IAPWS, 2011 <iapws-2011>`).


External Functions
------------------

The IAPWS-95 formulation uses density and temperature as state variables. For
most applications those sate variables are not the most convenient choices. Using
other sate variables requires solving an equation or two to get density and
temperature from the chosen state variables. This equation can have numerous
roots. Rather than solve these equations as part of the full process simulation,
external functions were developed that can solve the equations required to change
state variables and guarantee the correct roots.

The external property functions are written in C++ and complied such that they
can be called by AMPL solvers.  See the :ref:`installation instructions
<install:Installation Instructions>` for information about compiling these
functions.  The external functions provide both first and second derivatives for
all property function calls, however at phase transitions some of these functions
may be non-smooth.

IDAES Framework Wrappers
------------------------

A wrapper for these functions is provided with provides compatibility for with
the IDAES framework. Some methods for dealing with non-smoothness are also
included in the IDAES wrappers.

Pressure-Enthalpy
~~~~~~~~~~~~~~~~~

Temperature-Pressure-Vapor Fraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Expressions Common All State Variable Sets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. _iapws-2016:

International Association for the Properties of Water and Steam (2016).
IAPWS R6-95 (2016), "Revised Release on the IAPWS Formulation 1995 for
the Properties of Ordinary Water Substance for General Scientific Use,"
URL: http://iapws.org/relguide/IAPWS95-2016.pdf

.. _wagner-2002:

Wagner, W.,  A. Pruss (2002). "The IAPWS Formulation 1995 for the
Thermodynamic Properties of Ordinary Water Substance for General and
Scientific Use." J. Phys. Chem. Ref. Data, 31, 387-535.

.. _wagner-2000:

Wagner, W. et al. (2000). "The IAPWS Industrial Formulation 1997 for the
Thermodynamic Properties of Water and Steam," ASME J. Eng. Gas Turbines
and Power, 122, 150-182.

.. _akasaka-2008:

Akasaka, R. (2008). "A Reliable and Useful Method to Determine the Saturation
State from Helmholtz Energy Equations of State." Journal of Thermal
Science and Technology, 3(3), 442-451.

.. _iapws-2011:

International Association for the Properties of Water and Steam (2011).
IAPWS R15-11, "Release on the IAPWS Formulation 2011 for the
Thermal Conductivity of Ordinary Water Substance,"
URL: http://iapws.org/relguide/ThCond.pdf.

.. _iapws-2008:

International Association for the Properties of Water and Steam (2008).
IAPWS R12-08, "Release on the IAPWS Formulation 2008 for the Viscosity of
Ordinary Water Substance,"
URL: http://iapws.org/relguide/visc.pdf.
