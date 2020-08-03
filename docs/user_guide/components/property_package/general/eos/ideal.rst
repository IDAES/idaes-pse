Ideal Gases and Liquids (``Ideal``)
===================================

.. contents:: Contents 
    :depth: 2

Introduction
------------

Ideal behavior represents the simplest possible equation of state that ensures thermodynamic consistency between different properties.

Mass Density by Phase
---------------------

The following equation is used for both liquid and vapor phases, where :math:`p` indicates a given phase:

.. math:: \rho_{mass, p} = \rho_{mol, p} \times MW_p

where :math:`MW_p` is the mixture molecular weight of phase :math:`p`.

Molar Density by Phase
----------------------

For the vapor phase, the Ideal Gas Equation is used to calculate the molar density;

.. math:: \rho_{mol, Vap} = \frac{P}{RT}

whilst for the liquid phase the molar density is the weighted sum of the pure component liquid densities:

.. math:: \rho_{mol, Liq} = \sum_j{x_{Liq, j} \times \rho_{Liq, j}}

where :math:`x_{Liq, j}` is the mole fraction of component :math:`j` in the liquid phase.

Molar Enthalpy by Phase
-----------------------

For both liquid and vapor phases, the molar enthalpy is calculated as the weighted sum of the component molar enthalpies for the given phase:

.. math:: h_{mol, p} = \sum_j{x_{p, j} \times h_{mol, p, j}}

where :math:`x_{p, j}` is the mole fraction of component :math:`j` in the phase :math:`p`.

Component Molar Enthalpy by Phase
---------------------------------

Component molar enthalpies by phase are calculated using the pure component method provided by the users in the property package configuration arguments.

Molar Entropy by Phase
-----------------------

For both liquid and vapor phases, the molar entropy is calculated as the weighted sum of the component molar entropies for the given phase:

.. math:: s_{mol, p} = \sum_j{x_{p, j} \times s_{mol, p, j}}

where :math:`x_{p, j}` is the mole fraction of component :math:`j` in the phase :math:`p`.

Component Molar Entropy by Phase
---------------------------------

Component molar entropies by phase are calculated using the pure component method provided by the users in the property package configuration arguments.

Component Fugacity by Phase
---------------------------

For the vapor phase, ideal behavior is assumed:

.. math:: f_{Vap, j} = P

For the liquid phase, Raoult's Law is used:

.. math:: f_{Liq, j} = P_{sat, j}

Component Fugacity Coefficient by Phase
---------------------------------------

Ideal behavior is assumed, so all :math:`\phi_{p, j} = 1` for all components and phases.

Molar Gibbs Energy by Phase
---------------------------

For both liquid and vapor phases, the molar Gibbs energy is calculated as the weighted sum of the component molar Gibbs energies for the given phase:

.. math:: g_{mol, p} = \sum_j{x_{p, j} \times g_{mol, p, j}}

where :math:`x_{p, j}` is the mole fraction of component :math:`j` in the phase :math:`p`.

Component Gibbs Energy by Phase
-------------------------------

Component molar Gibbs energies are calculated using the definition of Gibbs energy:

.. math:: g_{mol, p, j} = h_{mol, p, j} - s_{mol, p, j} \times T
