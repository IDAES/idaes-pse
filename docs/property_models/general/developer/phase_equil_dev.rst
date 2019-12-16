Developing Phase Equilibrium Methods
====================================

.. contents:: Contents 
    :depth: 3

Handling phase equilibrium and phase transitions within an equation oriented framework can be challenging as it is necessary to ensure that all constraints and variables has feasible solution at all states. When dealing with disappearing phases and correlations that can become ill-defined or singular outside of the two phase envelope, it is necessary to either bound the problem to the two-phase region or reformulate the problem.

The IDAES Generic Property Package Framework provides support for reformulating the problem by defining an "equilibrium temperature" (`self._teq`) at which all phase equilibrium calculations are performed. Issues surrounding phase transitions can be avoided by  providing a definition for the equilibrium temperature that satisfies the following constraints:

.. math:: T_{\text{bubble}} <= T_{eq} <= T{\text{dew}}

The Phase Equilibrium module allows users to provide a definition for the equilibrium temperature, along with any necessary instructions on how to initialize the components associated with this definition.

A Phase Equilibrium module consists of two methods , which are described below.

`phase_equil(self)`
-------------------

The `phase_equil` method is responsible for defining the variables and constraints necessary for calculating the equilibrium temperature, and at a minimum must contain one constraint relating the equilibrium temperature (`self._teq`) to the system temperature (`self.temperature`).

`phase_equil_initialization(self)`
----------------------------------

This method is called by the Generic Property Package Framework initialization routine and should initialize the constraints associated with the phase equilibrium definition.

Note that the Generic Property Package Framework beings by deactivating all constraints in the `StateBlock` so the first step in the `phase_equil_initialization` method should be to activate any constraints defined in `phase_equil`. Additionally, this method may calculate initial values for any supporting variables defined in `phase_equil` based on variables that have already been initialized (primarily `temperature` and bubble and dew points if used). Developers should be careful however to fully understand the initialization sequence of the Generic Property Package Framework to understand which variables may have been initialized at this point.
