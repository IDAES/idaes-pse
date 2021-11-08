Solid Oxide Fuel Cell/Solid Oxide Electrolysis Cell (SOFC/SOEC)
===============================================================

.. Warning::

  This is a beta testing model and is still in development. Currently reactions
  are missing that allow for hydrocarbon reforming in the cell.  The primary
  current use of this model is modeling of SOECs.


This is a 2D steady-state or dynamic model of a solid oxide fuel or electrolysis
cell.

Modes
-----

The model can be operated in SOFC or SOEC mode in either countercurrent or
co-current flow. The model operate in isothermal model with an overall energy
balance or in non-isothermal mode with microscopic energy balances. Isothermal
mode may be appropriate for steady-state SOEC applications.

Model Structure
---------------

The model is made up of several sub-models described in the table below.
Equations that involve connections between the sub-models are defined at the
unit model level.  A generic channel and generic electrode model where used
to simplify the fc, oc, fe, and oe blocks.

+---------------------+---------------------+
| Block               | Description         |
+=====================+=====================+
| fc                  | fuel channel        |
+---------------------+---------------------+
| oc                  | oxygen channel      |
+---------------------+---------------------+
| fe                  | fuel electrode      |
+---------------------+---------------------+
| oe                  | oxygen electrode    |
+---------------------+---------------------+
| el                  | electrolyte         |
+---------------------+---------------------+
| ic                  | interconnect        |
+---------------------+---------------------+


Variables, Expressions, and Parameters
--------------------------------------



These quantities show up in the channel and electrode sub-models.  Some quantities
also appear in the electrolyte and interconnect.



+---------------------+------------------+------------------------+------------+-----------------------------------+
| Symbol              | Model Name       | Blocks                 | Units      | Description                       |
+=====================+==================+========================+============+===================================+
| :math:`P`           | pressure         | fc, oc, fe, oe         | Pa         | Pressure                          |
+---------------------+------------------+------------------------+------------+-----------------------------------+
| :math:`T`           | temperature      | fc, oc, fe, eo, el, ic | K          | Temperature                       |
+---------------------+------------------+------------------------+------------+-----------------------------------+



Model Equations
---------------

Mass Balances
~~~~~~~~~~~~~

General Equations
"""""""""""""""""

Concentration

..math ::

  c_i = \frac{x_i P}{T R}

Mole Fraction

.. math::

  1 = \sum_{i \in C} x_i


Channels
""""""""

In the channel, z is the dimensionless length in the direction of flow, with
:math:`z=0` at the fluid inlet and :math:`z=1` at the fluid outlet.

Channel area for flow.

.. math::

  A_c = l_x l_w

Flow

.. math ::

  f = A_c u / v

.. math ::

  \frac{\partial c_i}{\partial t} = -\frac{1}{l_c}\frac{\partial c_i u}{\partial z} - \frac{J_{xi}}{l_x}

Where :math:`J_{xi}` is the flux of component :math:`i` to the electrode in mol/m:math:`^2`/s.

Electrodes
""""""""""

Effective Diffusivity

..math ::

  D_{ie} = \frac{\epsilon}{\tau} D_{im}

.. math::

  


Energy Balances
~~~~~~~~~~~~~~~

Electrochemical Equations
~~~~~~~~~~~~~~~~~~~~~~~~~


Properties
~~~~~~~~~~

Ideal gas properties are used in this unit model.  The property package is
currently built into the model, and should be compatible with port containing
temperature, pressure, mole fraction, and molar flowrate. The inclusion of
property calculations in the unit model is an attempt to maximize efficiency
and tractability. This model is spatially discretized in two dimensions leading
to a large system of equations.  Future versions of this model will support
general property packages.

Pressure, Temperature, and Volume
"""""""""""""""""""""""""""""""""

The ideal gas law is used.

.. math::

  Pv = RT

Where :math:`v` is the molar volume.

Enthalpy and Entropy
""""""""""""""""""""

The ideal gas enthalpy and entropy for pure components are estimated using
relations and parameters from "NIST Chemistry Webbook" https://webbook.nist.gov/.

Pure component enthalpy is given by:

.. math ::

  H_i = A_i \frac{T}{1000} +
  \frac{B_i}{2} \left( \frac{T}{1000} \right)^2 +
  \frac{C_i}{3} \left( \frac{T}{1000} \right)^3 +
  \frac{D_i}{4} \left( \frac{T}{1000} \right)^4 -
  E_i \left( \frac{T}{1000} \right)^{-1} +
  F_i

Pure component entropy is given by:

.. math ::

  S_i = A_i \log_e \left( \frac{T}{1000} \right) +
  B_i \left( \frac{T}{1000} \right) +
  \frac{C_i}{2} \left( \frac{T}{1000} \right)^2 +
  \frac{D_i}{3} \left( \frac{T}{1000} \right)^3 -
  \frac{E_i}{2} \left( \frac{T}{1000} \right)^{-2} +
  G_i

Pure component internal energy is given by:

.. math ::

  U_i = H_i - Pv

The mixture enthalpy is:

.. math ::

  H_{m} = \sum_{i \in C} x_i H_i

The mixture internal energy is:

.. math ::

  U_{m} = \sum_{i \in C} x_i U_i

The mixture entropy is:

.. math ::

  S_{m} = \sum_{i \in C} \left( x_i S_{i} + R x_i \log_e x_i \right)


Diffusion Coefficients
""""""""""""""""""""""

Binary diffusion coefficients are calculated using Equation 11-3.2 and
characteristic length and Lenard-Jones energy parameters from "The Properties of
Gases and Liquids" 5th Ed. by Poling, Prausnitz, and O'Connell 2001.

.. math::

  D_{AB} = \frac{2.66T^\frac{3}{2}}{100 PM_{AB}\sigma_{AB}^2\Omega_D}

Where :math:`D_{AB}` is the binary diffusion coefficient for components A and B,
:math:`P` is pressure in Pa, and :math:`T` is the temperature in K.  Calculations
of other quantities are shown below.

.. math::

  M_{AB} = \left(\frac{1}{M_A} + \frac{1}{M_B}\right)^{-1}

Where :math:`M_A` is the molecular weight of component A in g/mol.

.. math::

  \varepsilon_{AB} = \left(\varepsilon_A \varepsilon_B \right)^\frac{1}{2}

Where :math:`\varepsilon_A` is the characteristic Lenard-Jones energy of component A.

.. math::

  \sigma_{AB} = \frac{\sigma_A + \sigma_B}{2}

Where :math:`\sigma_A` is the characteristic length of component A.

.. math::

  \Omega_D = \frac{1.06036}{\left(kT/\varepsilon_{AB}\right)^{0.15610}} +
    \frac{0.19300}{\exp\left(0.47635kT/\varepsilon_{AB}\right)} +
    \frac{1.035587}{\exp\left(1.52996kT/\varepsilon_{AB}\right)}
    \frac{1.76474}{\exp\left(3.89411kT/\varepsilon_{AB}\right)}

Where :math:`k` is Boltzmann's constant (:math:`k` is included in the tabulated
parameters for math:`\varepsilon`).

The diffusion coefficient of a component in a mixture is approximated by the
equation below taken from "Transport Phenomena," by Bird, Stewart and Lightfoot
1960. This assumes components other than :math:`i` move with the same velocity.

.. math::

  D_{im} = \left(1 - x_i \right)
    \left( \sum_{j \in C, j \ne i} \frac{x_j}{D_{ij}} \right)^{-1}

Where C is the set of components.
