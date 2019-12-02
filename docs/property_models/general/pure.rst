Defining Pure Component Properties
==================================

Most methods for calculating the thermophysical properties of materials start from estimating the properties of each component in its pure form, before applying mixing rules to determine the properties of the mixture. Pure component properties generally take the form of empirical correlations as a function of material state (generally temperature) derived from experimental data. Data and correlations for many component components are readily available in literature. However due to the empirical nature of these correlations and the wide range of data available, different sources use different forms for their correlations.

Within the IDAES Generic Property Package Framework, pure component property correlations are provided in the form of Python methods which return a Pyomo expression relating the pure component property to the material state (using the :ref:`standard naming conventions<standards:Standard Variable Names>`. IDAES provides a number of libraries of containing common forms for these correlations, and a list of the libraries currently supported by IDAES is given below.

A list of all the pure component properties currently supported by the IDAES Generic Property Package Framework can be found after the list of pure component libraries.

Pure Component Libraries
------------------------

.. toctree::
    :maxdepth: 1

    pure/NIST
    pure/Perrys
    pure/RPP

Supported Properties
--------------------

The following pure component properties are supported by IDAES Generic Property Package Framework.

.. csv-table::
   :header: "Property", "Method", "Arguments"

   "Ideal Gas Molar Heat Capacity", "`cp_mol_ig`", "component, temperature"
   "Ideal Gas Molar Enthalpy", "`enth_mol_ig`", "component, temperature"
   "Ideal Gas Molar Entropy", "`entr_mol_ig`", "component, temperature"
   "Ideal Liquid Molar Heat Capacity", "`cp_mol_liq`", "component, temperature"
   "Ideal Liquid Molar Enthalpy", "`enth_mol_liq`", "component, temperature"
   "Ideal Liquid Molar Entropy", "`entr_mol_liq`", "component, temperature"
   "Liquid Molar Density", "`dens_mol_liq`", "component, temperature"
   "Saturation Pressure", "`pressure_sat`", "component, temperature"
