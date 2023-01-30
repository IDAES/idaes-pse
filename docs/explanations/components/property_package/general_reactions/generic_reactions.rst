Defining Reaction Packages
==========================

.. contents:: Contents
    :depth: 2

Introduction
------------

In order to create and use a property package using the IDAES Generic Reaction Package Framework, users must provide a definition for the system they wish to model. The framework supports two approaches for defining the property package, which are described below, both of which are equivalent in practice.

Units of Measurement
--------------------

As with generic thermophysical property packages, when defining a reaction package using the generic framework users must define the base units for the reaction package (see :ref:`link<reference_guides/core/uom:Unit Sets>`). The approach for setting the base units and units for all parameters is the same as for thermophysical property packages and depends on the approach used to define the reaction package.

Config Dictionary
-----------------

The most common way to use the Generic Reaction Package Framework is to create an instance of the `GenericReactionParameterBlock` component and provide it with a dictionary of configuration arguments, as shown below:

.. code-block:: python

    m = ConcreteModel()

    m.fs = FlowsheetBlock()

    m.fs.thermo_properties = PhysicalParameterBlock()

    m.fs.reaction_properties = GenericReactionParameterBlock(property_package=m.fs.thermo_properties, **config_dict)

In the above example, the PhysicalParameterBlock object can be from any thermophysical property package suitable for the user's application.

Users need to populate `config_dict` with the desired options for their system as described in the other parts of this documentation. An example of a configuration dictionary can be found later on this page. For details on each configuration option, please see the relevant documentation.

Using this approach, units of measurement are defined using the `base_units` option in the configuration dictionary. Users must provide units for the 5 core quantities, and may also provide units for the other 2 SI base quantities (if required). For details on other configuration options, please see the relevant documentation.

Linking to a Thermophysical Property Package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As state information is defined by thermophysical property packages in IDAES, each reaction package must be linked to an appropriate thermophysical property package. This linkage is used by the reaction package to find the state information required to calculate the reaction properties, and thus the thermophysical property package must support all the properties required by the reaction package.

Setting Reaction Basis
^^^^^^^^^^^^^^^^^^^^^^

Many reaction properties (e.g. reaction rates) can be defined on different bases, such as a mass or molar basis. All properties within a package must use the same basis, which can be set using the "reaction_basis" configuration argument (see below). This must be done using the MaterialFlowBasis `Enum`, which can be imported from `idaes.core`.

Configuration Example
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from pyomo.environ import units as pyunits

    from idaes.core import MaterialFlowBasis

    config_dict = {
        "base_units": {"time": pyunits.s,
                       "length": pyunits.m,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "temperature": pyunits.K},
        "rate_reactions": {
            "R1": {"stoichiometry": {("Liq", "A"): -1,
                                     ("Liq", "B"): -1,
                                     ("Liq", "C"): 2},
                   "heat_of_reaction": constant_dh_rxn,
                   "rate_constant": arrhenius,
                   "rate_form": power_law_rate,
                   "concentration_form": ConcentrationForm.moleFraction,
                   "parameter_data": {
                       "dh_rxn_ref": (-10000, pyunits.J/pyunits.mol),
                       "arrhenius_const": (1, pyunits.mol/pyunits.m**3/pyunits.s),
                       "energy_activation": (1000, pyunits.J/pyunits.mol)}}},
        "equilibrium_reactions": {
            "R2": {"stoichiometry": {("Liq", "B"): -1,
                                     ("Liq", "C"): -1,
                                     ("Liq", "D"): 1},
                   "heat_of_reaction": constant_dh_rxn,
                   "equilibrium_constant": van_t_hoff,
                   "equilibrium_form": power_law_equil,
                   "concentration_form": ConcentrationForm.moleFraction,
                   "parameter_data": {
                       "dh_rxn_ref": (-20000, pyunits.J/pyunits.mol),
                       "k_eq_ref": (100, None),
                       "T_eq_ref": (350, pyunits.K)}}}}

Class Definition
----------------

Alternatively, the IDAES Generic Reaction Package Framework supports defining classes derived from the IDAES `GenericReactionParameterData` class with methods for defining configuration options and parameters.

Users can define two methods that are called automatically when an instance of the property package is created:

1. `configure`, which defines the users selection of sub-models, and
2. `parameters`, which defines the parameters necessary for the selected property methods.

A basic outline of a user defined Reaction Parameter Block is shown below.

.. code-block:: python

    @declare_process_block_class("UserReactionParameterBlock")
    class UserReactionParameterData(GenericReactionParameterData):
        def configure(self):
            # Set configuration options
            self.config.option_1 = value

        def parameters(self):
            # Define parameters
            self.param_1 = Var(index_set, initialize=value)

Users should populate the `configure` and `parameters` methods as discussed below.

Configure
^^^^^^^^^

The 'configure` method is used to assign values to the configuration arguments, using the format `self.config.option_name = value`. Users will also need to set the units of measurement in the property package metadata.

Parameters
^^^^^^^^^^

The `parameters` method is used to construct all the parameters associated with the property calculations and to specify values for these. The list of necessary parameters is based on the configuration options and the selected methods. Each method lists their necessary parameters in their documentation. Users need only define those parameters required by the options they have chosen.
