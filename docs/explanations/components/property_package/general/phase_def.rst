Defining Phases
===============

The second step in defining a property package using the Generic Property Package Framework is to define the phases of interest in the system. Due to the equation-oriented nature of the IDAES modeling framework, it is necessary to define any phases the user believes may be important *a priori* as it is not possible to determine what phases should be included on-the-fly. Phases are defined using `IDAES Phase objects<explanations/components/property_package/phase:Phase Object>`, and are automatically constructed using the `phases` configuration argument from the `GenericParameterBlock`.

The `phases` Argument
---------------------

Each `GenericParameterBlock` has a configuration argument named `phases` which is used to construct the `Phase` objects and populate them with instructions on how to calculate thermophysical properties for that phase. The `phases` configuration argument is expected to be a `dict` where the keys are the names for the phases of interest and the values are a configuration arguments for the named phase (which are passed to the `Phase` object as it is instantiated).

.. code-block:: python

    "phases": {
        "phase_1": {
            "type": Phase,
            "equation_of_state": EoS,
            "equation_of_state_options": {},
            "parameter_data": {}},
        "phase_2": {
            "type": Phase,
            "equation_of_state": EoS,
            "equation_of_state_options": {},
            "parameter_data": {}}}

Type Argument
^^^^^^^^^^^^^

Each phase in the `phases` argument must be assigned a valid phase type from those supported by the IDAES Framework (e.g. LiquidPhase, SolidPhase, VaporPhase). This should be provided using the `type` argument.

Equations of State
^^^^^^^^^^^^^^^^^^

Equations of state (or equivalent methods) describe the relationship between different thermophysical properties within a mixture and ensure that the behavior of these are thermodynamically consistent. Each phase must be assigned an Equation of State (or equivalent method) in the form of a Python module which will assemble the necessary variables, constraints and expressions associated with the desired approach.

A wide range of equations of states are available in literature for different applications and levels of rigor, and the IDAES Generic Property Package Framework provides a number of prebuilt modules for users, which are listed below.

Equation of state packages may allow for user options (e.g., choosing a specific type of cubic equation of state). The options are set using the `equation_of_state_options` argument, and the options available are described in the documentation of each equation of state module.

Equation of State Libraries
"""""""""""""""""""""""""""

.. toctree::
    :maxdepth: 1

    eos/ideal
    eos/cubic
    
Critical Properties of Mixtures
"""""""""""""""""""""""""""""""

Calculation of the critical properties of mixtures depends on the equation of state being used. As the Modular Property Package framework allows users to specify different equations of state for each phase, the following logic is used to determine which equation of state to use for calculating critical properties.

1. If a vapor-liquid equilibrium pair is defined in the `"phases_in_equilibrium"` configuration argument, then the liquid phase from this pair is used (IDAES generally assumed supercritical fluids are liquid like).
2. If no vapor-liquid equilibrium pair is defined, then the first liquid phase define is used.
3. If no liquid phases are defined then the vapor phase is used (it is assumed there will only be one vapor phase).
4. If no vapor phase is defined, then a `PropertyPackageError` is returned as there is no suitable phase for calculating critical properties.

Note that not all Equations of State are suitable for calculating critical properties as many cannot represent the critical conditions. The Ideal equation of State is one common example of an equation of state that DOES NOT support critical properties.

During initialization, the critical properties of the mixture are approximated using the mole fraction weighted sum of the component critical properties (note that users must provide values for all component critical properties (i.e., `compress_fact_crit`, `dens_mol_crit`, `pressure_crit` and `temperature_crit`) for initialization if critical properties are to be calculated.

Phase-Specific Parameter
^^^^^^^^^^^^^^^^^^^^^^^^

In some cases, a property package may include parameters which are specific to a given phase. In these cases, these parameters are stored as part of the associated `Phase` object and the values of these set using the `parameter_data` argument when declaring the phase. This is done in the same fashion as for :ref:`component specific parameters<explanations/components/property_package/general/component_def:Parameter Data>`.

Phases with Partial Component Lists
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In many applications a mixture will contain species that only appear in a single phase (either by nature or assumption). Common examples include crystalline solids and non-condensable gases. The IDAES Generic Property Package Framework provides support for these behaviors and allows users to specify phase-specific component lists (i.e., a list of components which appear in a given phase).

This is done by providing a phase with a `component_list` argument, which provides a `list` of component names which appear in the phase. The framework automatically validates the `component_list` argument to ensure that it is a sub-set of the master component list for the property package, and will inform the user if an unrecognized component is included. If a phase is not provided with a `component_list` argument it is assumed that all components defined in the master component list may be present in the phase.
