Defining Rate-Based Reactions
=============================


The `rate_reactions` Argument
-----------------------------

Each `GenericReactionParameterBlock` has a configuration argument named `rate_reactions` which is used to define rate-based reactions and specify how to calculate properties associated with these. The `rate_reactions` configuration argument is expected to be a dict-of-dicts, where the keys are the names for the rate-based reactions, and the values are a dict of configuration arguments for that reaction. Note that reaction names must be unique across both rate-based and equilibrium reactions, as all reactions are indexed by name.

.. code-block:: python

    "rate_reactions": {
        "reaction_1": {options},
        "reaction_2": {options}}

Configuration Arguments
-----------------------

The configuration arguments for each rate-based reaction are used to define methods for calculating reaction properties and defining the parameters associated with these. A full list of the supported configuration arguments is given below:

* stoichiometry (required)
* rate_form (required)
* concentration_form
* heat_of_reaction
* rate_constant


Stoichiometry
^^^^^^^^^^^^^

The `stoichiometry` configuration argument is used to define which components take part in a reaction, and is a required argument. The `stoichiometry` argument should be a dict where the keys are phase-component pairs and the values are the stoichiometric coefficient for that pair. Users need only provide values for those components that take part in the reaction - all undeclared phase-component pairs will be assumed to have a value of 0. An example of defining the reaction stoichiometry is given below, where in phase_1 component_1 is converted to component_2 in a 1:1 ratio:

.. code-block:: python

    "stoichiometry": {
        ("phase_1", "component_1"): -1,
        ("phase_1", "component_2"): 1}

Concentration Form
^^^^^^^^^^^^^^^^^^

Many common rate forms can be written using a number of different bases, such as molarity, molality or partial pressure. The `concentration_form` configuration argument is used in these cases to determine what basis to use for the concentration terms in the rate form and automatically write the correct expression (and determine units for the associated parameters. The `concentration_form` configuration argument must be an instance of a `ConcentrationForm` `Enum` (imported from idaes.models.properties.modular_properties.base.utility), and the following forms are currently available:

* molarity: ConcentrationForm.molarity
* activity: ConcentrationForm.activity
* molality: ConcentrationForm.molality
* mole fractions: ConcentrationForm.moleFraction
* mass fractions: ConcentrationForm.massFraction
* partial pressure: ConcentrationForm.partialPressure

Other Reaction Properties
^^^^^^^^^^^^^^^^^^^^^^^^^

The remaining configuration arguments are used to define how different properties should be calculated for each reaction. The `rate_form` argument is required, however all other properties need only be defined if needed for the user's application. These arguments should be provided as either Python functions or classes;

* functions are used for self-contained correlations with hard-coded parameters,
* classes are used for more generic correlations which require associated parameters.

A list of the libraries of methods available in the IDAES Framework can be found :ref:`here<explanations/components/property_package/general_reactions/method_libraries:Reaction Module Libraries>`.
