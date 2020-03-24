Developing Pure Component Methods
=================================

.. contents:: Contents 
    :depth: 3

The most common task developers of new property packages will need to do is writing methods for new pure component property calculations. Most equation of state type approaches rely on a set of calculations for pure components under ideal conditions which are then modified to account for mixing and deviations from ideality. These pure component property calculations tend to be empirical correlations based on experimental data (generally as functions of temperature) and due to their empirical nature a wide range of forms have been used in literature.

In order to support different forms for these calculations, the IDAES Generic Property Package Framework uses Python methods to define the form of pure component property calculations. This allows developers and users to easily enter the form they wish to use for their application with a minimum amount of code.

Naming Methods
--------------

The IDAES Generic Property Package Framework supports two ways of providing pure component property methods:

1. Providing the method directly - users may directly provide their method of choice as a config argument (`config.property_name`) in the `PropertyParameterBlock`, in which case the method can use any name the user desires.
2. Providing a library module - alternatively, users can provide a module containing a library of methods as the config argument (`config.property_name`), in which case the framework searches the module for a method with the same name as the property (and the config argument). E.g., for the property `enth_mol_phase_comp` the method name would be `enth_mol_phase_comp` (as would the associated config argument).

Method Arguments
----------------

.. note::

    Currently, the IDAES Generic Property Package Framework assumes pure component property calculations will be a function of only temperature. If additional functionality is required, please contact the IDAES Developers.

Currently, all pure component property methods in the IDAES Generic Property Package Framework take three arguments:

1. A reference to the `StateBlock` where the method will be used (generally `self`),
2. An element of a component list,
3. A pointer to the `temperature` variable to be used in the calculation. By using a pointer rather than an absolute reference (i.e. `self.temperature`), this allows the method to be applied at different temperatures as necessary (e.g. the reference temperature).

Method Parameters
-----------------

Pure component property methods all depend on a number of parameters, often derived from empirical data. In order to avoid duplication of parameters and facilitate parameter estimation studies, all property parameters are stored in the `PropertyParameterBlock` and each `StateBlock` contains a reference to its associated parameter block (`self.params`).

For pure component property methods, parameter names are define in the associated methods thus developers can choose any name they desire. However, the IDAES standard is to use the name of the property appended with `_coeff` and developers are encouraged to follow this convention.

Method Body
-----------

The body of the pure component property method should assemble an expression describing the specified quantity for the component given in the method arguments. This expression should involve Pyomo components from the `StateBlock` (i.e. `self`), the associated `PropertyParameterBlock` (`self.params`) and be returned in the final step of the method.

Example
-------

Below is an example of a pure component property method for the molar heat capacity of a component in the (ideal) gas phase with the form :math:`c_{\text{p, ig}, j} = A + B \times T`.

.. code:: python

    def cp_mol_ig_comp(self, component, temperature):
        # Method named using standard naming convention
        # Arguments are self, a component and temperature

        # Return an expression involving temperature and parameters
        return (self.params.cp_mol_ig_comp_coeff[component, "A"] +
                self.params.cp_mol_ig_comp_coeff[component, "B"]*temperature)

Note that the method only returns an expression representing the R.H.S. of the correlation.
