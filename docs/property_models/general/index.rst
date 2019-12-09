Generic Property Package Framework
==================================

Contents
--------

.. toctree::
    :maxdepth: 1

    generic_definition
    generic_options
    state_definition
    eos
    pure
    bubble_dew
    phase_equilibrium
    developers

Introduction
------------

.. note::
    The generic property package framework is still under development. Whilst the current framework is functional, features are still being developed and added to increase functionality.

    The generic property package framework builds upon the existing framework for implementing property packages within IDAES, and will not prevent the user of custom written property packages in the future. However, it is envisioned that the generic property package framework will provide a more streamlined interface for developing property packages in most circumstances, and it is hoped that most property packages will migrate to using the generic property framework in the future.

Property packages represent the core of any process model, and having a suitable property package is key to successfully modeling any process system. However, developing property packages is a significant challenge even for the most experienced modelers, as they involve complex, non-linear equations. The goal of the IDAES Generic Property Package Framework is to provide a flexible platform on which users can build custom property packages by calling upon libraries of modular sub-models to build up complex property calculations with the least effort possible.

The Generic Property Package Framework breaks down property packages into a number of components which can be assembled in a modular fashion. Users need only provide those components which they require for their system of interest, and components can be drawn from libraries of existing components or provided by the user as custom code. The components which make up a generic property package are as follows:

1. :ref:`Definition<property_models/general/generic_definition:Defining Property Packages>` of the component list and phases of interest, along with any phase equilibrium the user wishes to include.
2. A definition of the :ref:`variables<property_models/general/state_definition:Defining State Variables>` the user wishes to use to define the state of their material (state variables), along with any bounds on these.
3. An :ref:`equation of state<property_models/general/eos:Equations of State>` to describe each phase within the users property package.
4. :ref:`Correlations<property_models/general/pure:Defining Pure Component Properties>` for the pure component properties of each component in the users system. Correlations are only required for those properties the user will user within their model.
5. :ref:`Methods<property_models/general/bubble_dew:Bubble and Dew Point Methods>` to calculate the bubble and dew points of mixtures.
6. A :ref:`formulation<property_models/general/phase_equilibrium:Phase Equilibrium Formulations>` to use for defining any phase equilibrium within the users system.

The following section will describe how to define a property package using the Generic Property Package Framework along with the libraries of sub-models currently available. Finally, the :ref:`developers<property_models/general/developers:Developing New Property Libraries>` section describes how to go about defining your own custom components to use when creating custom property packages.

