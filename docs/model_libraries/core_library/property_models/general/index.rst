Generic Property Package Framework
==================================

Contents
--------

.. toctree::
    :maxdepth: 1

    generic_definition
    component_def
    phase_def
    state_definition
    phase_equilibrium
    developers

Introduction
------------

.. note::
    The generic property package framework is still under development. Whilst the current framework is functional, features are still being developed and added.

    The generic property package framework builds upon the existing framework for implementing property packages within IDAES, and will not prevent the user of custom written property packages in the future. Due to the complex nature of thermophysical property calculations, the generic property framework cannot support all possible materials and applications. Whilst it is hoped that the generic framework will be able to handle most common applications, users with more unusual systems or those solving computationally intensive problems may need to write custom property packages for their cases.

Property packages represent the core of any process model, and having a suitable property package is key to successfully modeling any process system. However, developing property packages is a significant challenge even for experienced modelers as they involve large numbers of tightly coupled constraints and parameters. The goal of the IDAES Generic Property Package Framework is to provide a flexible platform on which users can build property packages for common types of systems by calling upon libraries of modular sub-models to build up complex property calculations with the least effort possible.

The Generic Property Package Framework breaks down property packages into a number of components which can be assembled in a modular fashion. Users need only provide those components which they require for their system of interest, and components can be drawn from libraries of existing components or provided by the user as custom code. Details on how to set up the definition of a property package using the generic framework are given :ref:`here<model_libraries/core_library/property_models/general/generic_definition:Defining Property Packages>`.

The components which make up a generic property package are as follows:

1. Define the :ref:`components<model_libraries/core_library/property_models/general/component_def:Defining Components>` which make up the material of interest, including methods for calculating the pure component properties of interest in the system.
2. Define the :ref:`phases of interest<model_libraries/core_library/property_models/general/phase_def:Defining Phases>` for the application, including equations of state and other phase specific decisions.
3. Choose the set of :ref:`state variables<model_libraries/core_library/property_models/general/state_definition:State Definition>` you wish to use and a reference state for the system.
4. (Optional) Define any :ref:`phase equilibria<model_libraries/core_library/property_models/general/phase_equilibrium:Defining Phase Equilibria>` which occurs in the system and methods associated with calculating this.

The following sections will describe how to define a property package using the Generic Property Package Framework along with the libraries of sub-models currently available. Finally, the :ref:`developers<model_libraries/core_library/property_models/general/developers:Developing New Property Libraries>` section describes how to go about defining your own custom components to use when creating custom property packages.
