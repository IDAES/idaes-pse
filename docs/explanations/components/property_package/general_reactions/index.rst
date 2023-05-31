Generic Reaction Package Framework
==================================

Contents
--------

.. toctree::
    :maxdepth: 1

    generic_reactions
    rate_rxns
    equil_rxns
    method_libraries

Introduction
------------

.. note::
    The generic reaction package framework is still under development. Whilst the current framework is functional, features are still being developed and added.

    The generic reaction package framework builds upon the existing framework for implementing reaction packages within IDAES, and will not prevent the use of custom written reaction packages in the future. Whilst it is hoped that the generic framework will be able to handle most common applications, users with more unusual systems or those solving computationally intensive problems may need to write custom reaction packages for their cases.

The Generic Reaction Package Framework breaks down reaction packages into a number of components which can be assembled in a modular fashion. Users need only provide those components which they require for their system of interest, and components can be drawn from libraries of existing components or provided by the user as custom code. Details on how to set up the definition of a reaction package using the generic framework are given :ref:`here<explanations/components/property_package/general_reactions/generic_reactions:Defining Reaction Packages>`.

The components which make up a generic reaction package are as follows:

1. Choose a base set of :ref:`units of measurement<explanations/components/property_package/general_reactions/generic_reactions:Units of Measurement>` for the property package.
2. :ref:`Associate<explanations/components/property_package/general_reactions/generic_reactions:Linking to a Thermophysical Property Package>` the reaction package with an appropriate thermodynamic property package. The thermodynamic property package must use the same set of base units of measurement,
3. Define the :ref:`basis of the reaction terms<explanations/components/property_package/general_reactions/generic_reactions:Setting Reaction Basis>` for the reaction package.
4. Define the :ref:`rate-based reactions<explanations/components/property_package/general_reactions/rate_rxns:Defining Rate-Based Reactions>` of interest in the system.
5. Define the :ref:`equilibrium-based reactions<explanations/components/property_package/general_reactions/equil_rxns:Defining Equilibrium Reactions>` of interest in the system. Note that phase equilibrium is generally handled in the thermodynamic property package.

The following sections will describe how to define a reaction package using the Generic Reaction Package Framework along with the libraries of sub-models currently available. Finally, the :ref:`developers<explanations/components/property_package/general/developers:Developing New Property Libraries>` section describes how to go about defining your own custom components to use when creating custom property packages.

.. note::
   Within most IDAES models "parameters" are in fact defined as Pyomo 'Vars' (i.e. variables) which are fixed at their defined values. Whilst `Params` would seem to be the logical choice for these, parameter estimation problems require the parameters being estimated to be defined as `Vars` so that the solver is free to vary them.

