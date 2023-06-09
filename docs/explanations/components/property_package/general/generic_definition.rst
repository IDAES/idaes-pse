Defining Property Packages
==========================

.. contents:: Contents 
    :depth: 2

Introduction
------------

In order to create and use a property package using the IDAES Generic Property Package Framework, users must provide a definition for the material they wish to model. The framework supports two approaches for defining the property package that are described below, both of which are equivalent in practice.

Units of Measurement
--------------------

When defining a property package using the generic framework, users must define the base units for the property package (see :ref:`link<reference_guides/core/uom:Unit Sets>`). The approach for setting the base units depends on the approach used to define the property package, and is discussed in more detail in each section.

The Generic Property Package Framework includes the necessary code to convert between different units of measurement as required, allowing users to combine property methods with different sets of units into a single property package. In these cases, each property method is written in its natural units (including parameters), and the final result is automatically converted to the base units.

For example, the Antoine equation is generally written with pressure in bars and temperature in either Kelvin or Celsius (depending on source). Using the generic property framework, the users provide the Antoine coefficients in their original units (i.e. bar and Kelvin/Celsius) and the property calculation is written in these units. However, the final result (saturation pressure) is then converted to the base units specified in the property package definition.

Property Parameters
-------------------

Thermophysical property models all depend upon a set of parameters to describe the fundamental behavior of the system. For the purposes of the Generic Property Framework, these parameters are grouped into three types:

1. Component-specific parameters - these are parameters that are specific to a given chemical species, and are defined in the `parameter_data` argument for each component and stored in the associated `Component` block. Examples of these parameters include those used to calculate the ideal, pure component properties.
2. Phase-specific parameters - these are parameters that are specific to a given phase, and are defined in the `parameter_data` argument for each phase and stored in the associated `Phase` block. These types of parameters are relatively uncommon.
3. Package-wide parameters - these are parameters that are not necessarily confined to a single phase or species, and are defined in the `parameter_data` argument of the overall property package and stored in the main `Physical Parameter` block. Examples of these types of parameters include binary interaction parameters, which involve multiple species and can be used in multiple phases.

Config Dictionary
-----------------

The most common way to use the Generic Property Package Framework is to create an instance of the `GenericParameterBlock` component and provide it with a dictionary of configuration arguments, as shown below:

.. code-block:: python

    m = ConcreteModel()

    m.fs = FlowsheetBlock()

    m.fs.properties = GenericParameterBlock(default=config_dict)

Users need to populate `config_dict` with the desired options for their system as described in the other parts of this documentation. An example of a configuration dictionary for a benzene-toluene VLE system is shown below.

Using this approach, units of measurement are defined using the `base_units` option in the configuration dictionary. Users must provide units for the 5 core quantities, and may also provide units for the other 2 SI base quantities (if required). For details on other configuration options, please see the relevant documentation.

.. code-block:: python

    from pyomo.environ import units as pyunits

    config_dict = {
        "base_units": {"time": pyunits.s,
                       "length": pyunits.m,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "temperature": pyunits.K},
        "components": {
        'benzene': {"type": Component,
                    "elemental_composition": {"C": 6, "H": 6},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "pressure_sat_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                    "parameter_data": {
                        "mw": (78.1136E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (48.9e5, pyunits.Pa),  # [1]
                        "temperature_crit": (562.2, pyunits.K),  # [1]
                        "dens_mol_liq_comp_coeff": {
                            '1': (1.0162, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.2655, None),
                            '3': (562.16, pyunits.K),
                            '4': (0.28212, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (-3.392E1, pyunits.J/pyunits.mol/pyunits.K),  # [1]
                            'B': (4.739E-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-3.017E-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (7.130E-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (1.29E2, pyunits.J/pyunits.kmol/pyunits.K),  # [2]
                            '2': (-1.7E-1, pyunits.J/pyunits.kmol/pyunits.K**2),
                            '3': (6.48E-4, pyunits.J/pyunits.kmol/pyunits.K**3),
                            '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                            '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                        "enth_mol_form_liq_comp_ref": (
                            49.0e3, pyunits.J/pyunits.mol),  # [3]
                        "enth_mol_form_vap_comp_ref": (
                            82.9e3, pyunits.J/pyunits.mol),  # [3]
                        "pressure_sat_comp_coeff": {'A': (-6.98273, None),  # [1]
                                                    'B': (1.33213, None),
                                                    'C': (-2.62863, None),
                                                    'D': (-3.33399, None)}}},
        'toluene': {"type": Component,
                    "elemental_composition": {"C": 7, "H": 8},
                    "dens_mol_liq_comp": Perrys,
                    "enth_mol_liq_comp": Perrys,
                    "enth_mol_ig_comp": RPP,
                    "pressure_sat_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                    "parameter_data": {
                        "mw": (92.1405E-3, pyunits.kg/pyunits.mol),  # [1]
                        "pressure_crit": (41e5, pyunits.Pa),  # [1]
                        "temperature_crit": (591.8, pyunits.K),  # [1]
                        "dens_mol_liq_comp_coeff": {
                            '1': (0.8488, pyunits.kmol*pyunits.m**-3),  # [2] pg. 2-98
                            '2': (0.26655, None),
                            '3': (591.8, pyunits.K),
                            '4': (0.2878, None)},
                        "cp_mol_ig_comp_coeff": {
                            'A': (-2.435E1, pyunits.J/pyunits.mol/pyunits.K),  # [1]
                            'B': (5.125E-1, pyunits.J/pyunits.mol/pyunits.K**2),
                            'C': (-2.765E-4, pyunits.J/pyunits.mol/pyunits.K**3),
                            'D': (4.911E-8, pyunits.J/pyunits.mol/pyunits.K**4)},
                        "cp_mol_liq_comp_coeff": {
                            '1': (1.40E2, pyunits.J/pyunits.kmol/pyunits.K),  # [2]
                            '2': (-1.52E-1, pyunits.J/pyunits.kmol/pyunits.K**2),
                            '3': (6.95E-4, pyunits.J/pyunits.kmol/pyunits.K**3),
                            '4': (0, pyunits.J/pyunits.kmol/pyunits.K**4),
                            '5': (0, pyunits.J/pyunits.kmol/pyunits.K**5)},
                        "enth_mol_form_liq_comp_ref": (
                            12.0e3, pyunits.J/pyunits.mol),  # [3]
                        "enth_mol_form_vap_comp_ref": (
                            50.1e3, pyunits.J/pyunits.mol),  # [3]
                        "pressure_sat_comp_coeff": {'A': (-7.28607, None),  # [1]
                                                    'B': (1.38091, None),
                                                    'C': (-2.83433, None),
                                                    'D': (-2.79168, None)}}}},
        "phases":  {'Liq': {"type": LiquidPhase,
                            "equation_of_state": ideal},
                    'Vap': {"type": VaporPhase,
                            "equation_of_state": ideal}},
        "state_definition": FcPh,
        "state_bounds": {
            # Note format is (lower, nominal, upper, units)
            "flow_mol": (0, 100, 1000, pyunits.mol/pyunits.s),
            "temperature": (273.15, 300, 450, pyunits.K),
            "pressure": (5e4, 1e5, 1e6, pyunits.Pa)},
        "pressure_ref": (1e5, pyunits.Pa),
        "temperature_ref": (300, pyunits.K),
        "phases_in_equilibrium": [("Vap", "Liq")],
        "phase_equilibrium_state": {("Vap", "Liq"): smooth_VLE},
        "bubble_dew_method": IdealBubbleDew}

Data Sources:

1. The Properties of Gases and Liquids (1987), 4th edition, Chemical Engineering Series - Robert C. Reid
2. Perry's Chemical Engineers' Handbook 7th Ed.
3. Engineering Toolbox, https://www.engineeringtoolbox.com, Retrieved 1st December, 2019

Class Definition
----------------

Alternatively, the IDAES Generic Property Package Framework supports defining classes derived from the IDAES `GenericParameterData` with methods for defining configuration options and parameters.

Users can define two methods that are called automatically when an instance of the property package is created:

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

The 'configure` method is used to assign values to the configuration arguments, using the format `self.config.option_name = value`. Users will also need to set the units of measurement in the property package metadata.

Parameters
^^^^^^^^^^

The `parameters` method is used to construct all the parameters associated with the property calculations and to specify values for these. The list of necessary parameters is based on the configuration options and the selected methods. Each method lists their necessary parameters in their documentation. Users need only define those parameters required by the options they have chosen.

Examples
--------

Examples of using the IDAES Generic Property Package Framework can be found in the `idaes/property_models/core/examples` folder.
