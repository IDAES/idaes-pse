Configuration Options
=====================

.. contents:: Contents 
    :depth: 2

The following configuration options are available in the IDAES Generic Property Package Framework.

Mandatory Configuration Options
-------------------------------

Users must provide a selection for the following options in all property packages.

Component List
^^^^^^^^^^^^^^
**Argument:** `config.component_list`

The list of chemical species of interest in the material.

Phase List
^^^^^^^^^^
**Argument:** `config.phase_list`

The list of thermodynamic phases that should be included in the model. Phases may or may not interact via phase equilibrium.

State Definition
^^^^^^^^^^^^^^^^
**Argument:** `config.state_definition`

An IDAES state definition module which creates the desired set of state variables along with any necessary auxiliary variables.

Equation of State
^^^^^^^^^^^^^^^^^
**Argument:** `config.equation_of_state`

A `dict` indicating an equation of state to use for each phase in the property package. The expected form is:

.. code:: python

    config.equation_of_state = {'phase_1': eos_1, 'phase_2': eos_2, ...}

Each phase in `config.phase_list` must be assigned an equation of state, which should take the form of an IDAES equation of state module which defines methods for calculating all thermophysical and transport properties.

Additional Configuration Options
--------------------------------

The following configuration options are not necessary, but are useful or required in some circumstances.

Phase Component Dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^
**Argument:** `config.phase_component_dict`

The option allows users to specify different component lists for each phase in their system. This is useful in circumstances where certain species will only ever appear in a given phase (e.g. a non-condensible gas). The expected form of this argument is:

.. code:: python

    config.phase_component_dict = {'phase_1': [list of components in phase 1], 'phase_2': [list of components in phase 2], ...}

Component lists for each phase must be a subset of `config.component_list`, and all components in `config.component_list` should appear in at least one phase.

State Bounds
^^^^^^^^^^^^
**Argument:** `config.state_bounds`

The option allows users to specify custom bounds on the state variables in their property package during construction. This is important for bounding the resulting problem and ensuring solutions do not stray outside the regions over which property parameters were fitted. The expected form of this argument is:

.. code:: python

    config.state_bounds = {'state_var_1': (lower, upper), 'state_var_2': (lower, upper), ...}

Users should consult the documentation for the state definition they are using to determine the state variables which can be bounded.

Phase Equilibrium Formulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Argument:** `config.phase_equilibrium_formulation`

The option allows users to specify the formulation to use for expressing phase equilibrium in their property package. This argument should be an IDAES phase equilibrium module which creates constraint describing the equilibria between phases. If the user wishes to include phase equilibria in their property package, both this argument and the `phase_equilibrium_dict` argument must be provided.

Phase Equilibrium Dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Argument:** `config.phase_equilibrium_dict`

The option allows users to specify which components in their system are in equilibrium between different phases. The expected form of this argument is:

.. code:: python

    config.phase_equilibrium_dict = {id1: [component, (phase_1, phase_2)], id2: [component, (phase_1, phase_2)], ...}

Here the `id` is used to identify each phase equilibrium reaction, `component` identifies the component in equilibrium and `phase1` and `phase2` identify the two phases over which this component should be in equilibrium. For cases where a given component is in equilibrium across more than 2 phases, multiple entries for the component are required identifying each pair of phases which should be in equilibrium (this is the reason for the `id` to identify reactions rather than just component name).

If the user wishes to include phase equilibria in their property package, both this argument and the `phase_equilibrium_formulation` argument must be provided.

Bubble Temperature
^^^^^^^^^^^^^^^^^^
**Argument:** `config.temperature_bubble`

This argument allows users to specify a method for calculating the bubble temperature of the mixture in their property package.

Dew Temperature
^^^^^^^^^^^^^^^
**Argument:** `config.temperature_dew`

This argument allows users to specify a method for calculating the dew temperature of the mixture in their property package.

Bubble Pressure
^^^^^^^^^^^^^^^
**Argument:** `config.pressure_bubble`

This argument allows users to specify a method for calculating the bubble pressure of the mixture in their property package.

Dew Pressure
^^^^^^^^^^^^
**Argument:** `config.pressure_dew`

This argument allows users to specify a method for calculating the dew pressure of the mixture in their property package.

Pure Component Property Options
-------------------------------

The remaining options allow users to select methods to use for calculating each pure component property, and users must provide a selection for every method that will be used within their process flowsheet. A full list of supported pure component properties can be found :ref:`here<model_libraries/core_library/property_models/general/pure:Supported Properties>`.


