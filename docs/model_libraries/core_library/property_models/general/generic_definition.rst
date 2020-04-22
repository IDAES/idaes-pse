Defining Property Packages
==========================

.. contents:: Contents 
    :depth: 2

Introduction
------------

In order to create and use a property package using the IDAES Generic Property Package Framework, users must provide a definition for the material they wish to model. The framework supports two approaches for defining the property package, which are described below, both of which are equivalent in practice.

Config Dictionary
-----------------

The most common way to use the Generic Property Package Framework is to create an instance of the `GenericParameterBlock` component and provide it with a dictionary of configuration arguments, as shown below:

.. code-block:: python

    m = ConcreteModel()

    m.fs = FlowsheetBlock()

    m.fs.properties = GenericParameterBlock(default=config_dict)

Users need to populate `config_dict` with the desired options for their system as described in the other parts of this documentation. An example of a configuration dictionary for a benzene-toluene VLE system is shown below. For details on each configuration option, please see the relevant documentation.

.. code-block:: python

    config_dict = {
        "components": {
            'benzene': {
                "dens_mol_liq_comp": Perrys,
                "enth_mol_liq_comp": Perrys,
                "enth_mol_ig_comp": RPP,
                "pressure_sat_comp": RPP,
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                "parameter_data": {
                    "mw": 78.1136E-3,  # [1]
                    "pressure_crit": 48.9e5,  # [1]
                    "temperature_crit": 562.2,  # [1]
                    "dens_mol_liq_comp_coeff": {'1': 1.0162*1e3,  # [2] pg. 2-98
                                                '2': 0.2655,
                                                '3': 562.16,
                                                '4': 0.28212},
                    "cp_mol_ig_comp_coeff": {'A': -3.392E1,  # [1]
                                             'B': 4.739E-1,
                                             'C': -3.017E-4,
                                             'D': 7.130E-8},
                    "cp_mol_liq_comp_coeff": {'1': 1.29E2,  # [2]
                                              '2': -1.7E-1,
                                              '3': 6.48E-4,
                                              '4': 0,
                                              '5': 0},
                    "enth_mol_form_liq_comp_ref": 49.0e3,  # [3]
                    "enth_mol_form_vap_comp_ref": 82.9e3,  # [3]
                    "pressure_sat_comp_coeff": {'A': -6.98273,  # [1]
                                                'B': 1.33213,
                                                'C': -2.62863,
                                                'D': -3.33399}}},
            'toluene': {
                "dens_mol_liq_comp": Perrys,
                "enth_mol_liq_comp": Perrys,
                "enth_mol_ig_comp": RPP,
                "pressure_sat_comp": RPP,
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                "parameter_data": {
                    "mw": 92.1405E-3,  # [1]
                    "pressure_crit": 41e5,  # [1]
                    "temperature_crit": 591.8,  # [1]
                    "dens_mol_liq_comp_coeff": {'1': 0.8488*1e3,  # [2] pg. 2-98
                                                '2': 0.26655,
                                                '3': 591.8,
                                                '4': 0.2878},
                    "cp_mol_ig_comp_coeff": {'A': -2.435E1,
                                             'B': 5.125E-1,
                                             'C': -2.765E-4,
                                             'D': 4.911E-8},
                    "cp_mol_liq_comp_coeff": {'1': 1.40E2,  # [2]
                                              '2': -1.52E-1,
                                              '3': 6.95E-4,
                                              '4': 0,
                                              '5': 0},
                    "enth_mol_form_liq_comp_ref": 12.0e3,  # [3]
                    "enth_mol_form_vap_comp_ref": 50.1e3,  # [3]
                    "pressure_sat_comp_coeff": {'A': -7.28607,  # [1]
                                                'B': 1.38091,
                                                'C': -2.83433,
                                                'D': -2.79168}}}},
        "phases":  {'Liq': {"type": LiquidPhase,
                            "equation_of_state": ideal},
                    'Vap': {"type": VaporPhase,
                            "equation_of_state": ideal}},
        "state_definition": FcPh,
        "state_bounds": {"flow_mol_comp": (0, 1000),
                         "temperature": (273.15, 450),
                         "pressure": (5e4, 1e6),
                         "enth_mol": (1e4, 2e5)},
        "pressure_ref": 1e5,
        "temperature_ref": 300,
        "phases_in_equilibrium": [("Vap", "Liq")],
        "phase_equilibrium_formulation": {("Vap", "Liq"): smooth_VLE},
        "temperature_bubble": bubble_temp_ideal,
        "temperature_dew": dew_temp_ideal,
        "pressure_bubble": bubble_press_ideal,
        "pressure_dew": dew_press_ideal}

Data Sources:

1. The Properties of Gases and Liquids (1987), 4th edition, Chemical Engineering Series - Robert C. Reid
2. Perry's Chemical Engineers' Handbook 7th Ed. (converted to J/mol.K, mol/m^3)
3. Engineering Toolbox, https://www.engineeringtoolbox.com, Retrieved 1st December, 2019

Class Definition
----------------

Alternatively, the IDAES Generic Property Package Framework supports defining classes derived from the IDAES `GenericParameterData` with methods for defining configuration options and parameters.

Users can define two methods which are called automatically when an instance of the property package is created:

1. `configure`, which defines the users selection of sub-models, and
2. `parameters`, which defines the parameters necessary for the selected property methods.

A basic outline of a user defined Property Parameter Block is shown below.

.. code-block:: python

    @declare_process_block_class("UserParameterBlock")
    class UserParameterData(GenericParameterData):
        def configure(self):
            # Set configuration options
            self.config.option_1 = value

        def parameters(self):
            # Define parameters
            self.param_1 = Var(index_set, initialize=value)

Users should populate the `configure` and `parameters` methods as discussed below.

Configure
^^^^^^^^^

The 'configure` method is used to assign values to the configuration arguments, using the format `self.config.option_name = value`.

Parameters
^^^^^^^^^^

The `parameters` method is used to construct all the parameters associated with the property calculations and to specify values for these. The list of necessary parameters is based on the configuration options and the selected methods. Each method lists their necessary parameters in their documentation. Users need only define those parameters required by the options they have chosen.

Property parameters can be defined as either Pyomo `Params` or `Vars` depending upon the users needs and application. Whilst `Params` would seem to be the logical choice, be aware that for parameter estimation problems the parameters being estimated need to be defined as `Vars` (so that the solver is free to vary them). 

.. note::

   If using `Params`, users should consider whether these should be `mutable` or not - `Params` that are not mutable have their value defined upon creation and this cannot be changed later.

   If using `Vars`, remember that you will need to fix the value unless you are trying to estimate the value of that parameter.

Examples
--------

Examples of using the IDAES Generic Property Package Framework can be found in the `idaes/property_models/core/examples` folder.
