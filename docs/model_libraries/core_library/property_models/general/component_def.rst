Defining Components
===================

The first step in defining a generic property package is to describe each of the chemical species of interest within the system, including methods for calculating the necessary thermophysical properties of the pure component. Components are defined using :ref:`IDAES Component objects<core/comp:Component Objects>`, and are automatically constructed using the `components` configuration argument from the `GenericParameterBlock`.

The `components` Argument
-------------------------

Each `GenericParameterBlock` has a configuration argument named `components` which is used to construct the `Component` objects and populate them with instructions on how to calculate thermophysical properties for that component. The `components` configuration argument is expected to be a dict-of-dicts, where the keys are the names for the chemical species of interest, and the values are a dict of configuration arguments for the named component (which are passed to the `Component` object as it is instantiated).

.. code-block:: python

    "components": {
        "species_1": {options},
        "species_2": {options}}

Configuration Arguments
-----------------------

The configuration arguments for each chemical species are used to define methods for calculating pure component properties and defining the parameters associated with these. A full list of the supported configuration arguments for `Component` objects can be found :ref:`here<core/comp:Component Class>`.

Valid Phases
^^^^^^^^^^^^

In many cases, a given chemical species can only exist in certain phases; the most common example being ionic solids which dissociate upon dissolution (thus forming new ionic species in an aqueous phase). For each component, the user can set a list of the valid phase types for the component (liquid, vapor and/or solid) using the `valid_phase_types" configuration argument. This configuration argument should be a list containing `PhaseType` `Enums` (imported from `idaes.core.phases`) indicating the types of phases in which this component can exist.

This information is used by the Generic Property Framework to automatically determine the valid phase-component pairs for the user defined system. Users can override this automatic definition by providing a component list for a given phase in the definition of each `Phase` as discussed later (note however that user-defined phase-component lists are validated against the valid phases, and an exception will be raised if a component is assigned in a phase for which it is not valid).

Pure Component Property Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Most methods for calculating the thermophysical properties of materials start from estimating the properties of each component in its pure form, before applying mixing rules to determine the properties of the mixture. Pure component properties generally take the form of empirical correlations as a function of material state (generally temperature) derived from experimental data. Data and correlations for many components are readily available in literature. However due to the empirical nature of these correlations and the wide range of data available, different sources use different forms for their correlations.

Within the IDAES Generic Property Package Framework, pure component property correlations can be provided as either Python functions or classes;

* functions are used for self-contained correlations with hard-coded parameters,
* classes are used for more generic correlations which require associated parameters.

When providing a method via the `components` configuration argument, users can either provide a pointer to the desired class/method directly, or to a Python module containing a class or method with the same name as the property to be calculated. More details on the uses of these and how to construct your own can be found in the :ref:`developer documentation<model_libraries/core_library/property_models/general/developers:Developing New Property Libraries>`.

Pure Component Libraries
""""""""""""""""""""""""

As a starting point for users, the IDAES Generic Property Package Framework contains a library of some common methods for calculating properties of interest. These libraries are organized by source, and are listed below.

.. note::
    Users should be careful about mixing-and-matching methods from different libraries, especially for the same component. Thermodynamic properties are intrinsically coupled, thus many correlations are also linked and often share parameters. Mixing-and-matching correlations may result in two correlations using parameters with the same name but with different expectations.

    Additionally, sources often use different approaches for defining the thermodynamic reference state of the material, thus users need to ensure that a consistent reference state is being used when combining methods from different sources.

.. toctree::
    :maxdepth: 1

    pure/NIST
    pure/Perrys
    pure/RPP

Phase Equilibrium Formulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For those applications involving phase equilibria, there are number of different approaches that can be taken to specify the equilibrium condition. For example, equilibrium may be described in terms of an empirical partitioning coefficient or in terms of fugacities in each phase. To allow users to specify the approach they wish to use, each `Component` object contains a `phase_equilibrium_form` configuration argument.

As a given system may incorporate multiple phase equilibria, the `phase_equilibrium_form` argument should be a `dict` with keys beings a tuple of interacting phases and values being a Python method describing how the equilibrium condition should be defined. A simple example for a VLE system is shown below:

.. code-block:: python

    "phase_equilibrium_form": {("Vap", "Liq"): fugacity}

The IDAES Generic Property Package Framework contains a library of common forms for the equilibrium condition, which is described :ref:`here<model_libraries/core_library/property_models/general/pe/pe_forms:Library of Common Equilibrium Forms>`.

Parameter Data
^^^^^^^^^^^^^^

Most pure component property correlations depend upon empirical parameters which need to be specified by the user. All the in-built property libraries built these parameters automatically expect the user to provide values these parameters via the `parameter_data` configuration argument. The `parameter_data` configuration argument should be a `dict` with keys being the name of the required parameters and the values being a value or dict of values to use when initializing the parameter (i.e. the dict must have keys which match the indexing set of the parameter).

.. note::
    A `dict` is used for specifying parameter values to allow users greater flexibility in defining their own methods with custom parameters.

Additionally, the following quantities are properties of the component (i.e. not a function of state) and are included in the component parameters.

* Molecular weight: "mw"
* Critical Pressure: "pressure_crit"
* Critical Temperature: "temperature_crit"
