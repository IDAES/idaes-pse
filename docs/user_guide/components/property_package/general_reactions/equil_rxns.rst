Defining Equilibrium Reactions
==============================


The `equilibrium_reactions` Argument
------------------------------------

Each `GenericReactionParameterBlock` has a configuration argument named `equilibrium_reactions` which is used to define equilibrium reactions and specify how to calculate properties associated with these. The `equilibrium_reactions` configuration argument is expected to be a dict-of-dicts, where the keys are the names for the equilibrium reactions, and the values are a dict of configuration arguments for that reaction. Note that reaction names must be unique across both rate-based and equilibrium reactions, as all reactions are indexed by name.

.. code-block:: python

    "equilibrium_reactions": {
        "reaction_1": {options},
        "reaction_2": {options}}

Configuration Arguments
-----------------------

The configuration arguments for each equilibrium reaction are used to define methods for calculating reaction properties and defining the parameters associated with these. A full list of the supported configuration arguments is given below:

* stoichiometry (required)
* equilibrium_form (required)
* concentration_form
* heat_of_reaction
* equilibrium_constant

Stoichiometry
^^^^^^^^^^^^^

The `stoichiometry` configuration argument is used to define which components take part in a reaction, and is a required argument. The `stoichiometry` argument should be a dict where the keys are phase-component pairs and the values are the stoichiometric coefficient for that pair. Users need only provide values for those components that take part in the reaction - all undeclared phase-component pairs will be assumed to have a value of 0. An example of defining the reaction stoichiometry is given below, where in phase_1 component_1 is converted to component_2 in a 1:1 ratio:

.. code-block:: python

    "stoichiometry": {
        ("phase_1", "component_1"): -1,
        ("phase_1", "component_2"): 1}

Concentration Form
^^^^^^^^^^^^^^^^^^

See :ref:`rate reaction<user_guide/components/property_package/general_reactions/rate_rxns:Concentration Form>` documentation.

Other Reaction Properties
^^^^^^^^^^^^^^^^^^^^^^^^^

The remaining configuration arguments are used to define how different properties should be calculated for each reaction. The `equilibrium_form` argument is required, however all other properties need only be defined if needed for the user's application. These arguments should be provided as either Python functions or classes;

* functions are used for self-contained correlations with hard-coded parameters,
* classes are used for more generic correlations which require associated parameters.

A list of the libraries of methods available in the IDAES Framework can be found :ref:`here<user_guide/components/property_package/general_reactions/method_libraries:Reaction Module Libraries>`.

