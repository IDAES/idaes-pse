Developing Equation of State Modules
====================================

.. contents:: Contents 
    :depth: 3

The central part of any property package are the equations of state or equivalent models which describe how the mixture behaves under the conditions of interest. For systems with multiple phases and phase equilibrium, each phase must have its own equation of state (or equivalent), which must provide information on phase equilibrium which is compatible with the other phases in the system.

Equations of State and Multiple Phases
--------------------------------------

The IDAES Generic Property Package Framework requires users to assign an equation of state module for each phase in their system, thus equations of state can be written for specific phases (e.g. an ideal gas equation of state). In some cases, developers may wish to write equations of state for multiple phases, and the generic framework supports this by indexing all properties by phase.

Developers are encouraged to add checks to their methods to ensure their equations of state are only applied to phases where they are appropriate (e.g. an ideal gas equation of state should raise an exception if the phase argument is not "Vap").

General Structure
-----------------

Equation of State Modules in the IDAES Generic Property Package Framework are files (modules) containing a number of methods which describe the behavior of the material. These method define how each of the properties associated with a given phase should be calculated, and the list of properties supported for a given phase is limited by the methods provided by the developer of the equation of state.

Phase Equilibrium
-----------------

When calculating phase equilibrium, the IDAES Generic Property Package Framework uses the general form :math:`\Phi^e_{\text{phase 1, j}} = \Phi^e_{\text{phase 2, j}}` where :math:`\Phi^e_{p, j}` is the fugacity of component :math:`j` in phase :math:`p` calculated at the equilibrium temperatures (:math:`T_{eq}`, variable name `self._teq`). The equilibrium temperature is calculated using the users' choice of phase equilibrium formulation and determines how the property package will handle phase transitions.

All equation of state methods should contain a method for calculating fugacity if they are to support phase equilibrium calculations.

Accessing Pure Component Property Methods
-----------------------------------------

In most cases, property calculations in the equation of state methods will require calculations of the pure component properties for the system. These can be accessed using `get_method` (imported from from `idaes.property_models.core.generic.generic_property`) using the form `get_method(self, "property_name")`. This will return the **method** defined by the user in the `PropertyParameterBlock` for the named property, which can then be used in the equation of state methods (note that users will need to call the method and provide it with the required arguments - generally `self`, component and a pointer to temperature).

Common Methods
--------------

For equations of state that support multiple phases, there may be certain calculations and/or variables that are common to all phases. To support this (and avoid duplication of these), equation of state methods should contain a method named `common` which implements any component which are common to multiple phases. This method should also contain checks to ensure that these components have not already been created for another phase in the system (to avoid duplication). In cases where there are no common components, this method can `pass`.

Mixture Property Methods
------------------------

The main part of an equation of state method are a set of methods which describe properties of the mixture for a given phase. Any mixture property that the property package needs to support must be defined as a method in the equation of state module, which returns an expression for the given property (construction of the actual Pyomo component will be handled by the core framework code).

Mixture properties can be defined in any way the developer desires, and can cross-link and reference other mixture properties as required. Developers should recall that the State Definition method should have defined the following properties which can be used in mixture property correlations:

* pressure
* temperature
* mole_frac_phase_comp
* phase_frac

Other state variables **may** have been defined by the user's choice of State Definition, however this cannot be guaranteed. Developers may chose to assume that certain state variables will be present, but this will limit the application of their equation of state module to certain state definitions which should be clearly documented.

Example
-------

Below is an example method for a method in an equation of state module for calculating molar density that supports both liquid and vapor phases.

.. code:: python

    def dens_mol_phase(b, phase):
        if phase == "Vap":
            return b.pressure/(b.params.gas_const*b.temperature)
        elif phase == "Liq":
            return sum(b.mole_frac_phase_comp[phase, j] *
                       get_method(b, "dens_mol_liq_comp")(b, j, b.temperature)
                       for j in b.params.component_list)
        else:
            raise PropertyNotSupportedError("Phase not supported")
