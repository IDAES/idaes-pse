Developing State Definitions
============================

.. contents:: Contents 
    :depth: 3

The primary purpose of the State Definition method is to define the state variables which will be used to describe the state of the mixture in the property package. However, a number of other key aspects of the property package definition are tied to the choice of state variables and must be declared here as well.

State definitions are defined as Python modules with two methods and one list, which are describe below.

`define_state(self)`
--------------------

The first method in a State Definition module is the `define_state` method. This method is used to define the state variables and associated components and methods. The `define_state` method must define the following things:

State Variables
^^^^^^^^^^^^^^^

The most important part of a State Definition module is the definition of the state variables that should be used in the resulting property package. The choice of state variables is up to the module developer, however the set of variables selected must contain sufficient information to fully define the extensive and intensive state of the material. That is, if all the state variables are fixed, the resulting set of variables and constraints should form a square problem (i.e. 0 degrees of freedom). Beyond this requirement however, developers may choose any combination of state variables they wish.

State variables should be defined as Pyomo `Vars` with names drawn from the IDAES naming standard, and should include initial values and bounds. The Generic Property Package Framework includes an optional user input of bounds for the state variables (`config.state_bounds`) which developers are encouraged to make use of when setting bounds and initializing variables.

`define_state_vars`
^^^^^^^^^^^^^^^^^^^

In order to inform the IDAES Process Modeling Framework of which variables should be considered state variable, developers are required to define a method named `define_state_vars`. This method should return a `dict` where the keys are a string identifier for each state variable and the values being pointers to the associated `Var` component. For example:

.. code:: python

    def define_state_vars_state_definition():
        return {"flow_mol": self.flow_mol,
                "mole_frac_comp": self.mole_frac_comp,
                "pressure": self.pressure,
                "temperature": self.temperature,}
    self.define_state_vars = define_state_vars_state_definition

Auxiliary Variables
^^^^^^^^^^^^^^^^^^^

Whilst the developer is free to choose any set of state variable they wish to define their system, there are certain properties/quantities associated with material state that are frequently used in process models. For example, most property calculation methods drawn upon empirical correlations for pure component properties which are most commonly expressed as functions of temperature (and sometimes pressure). Additionally, multiphase systems often require knowledge of the volume fractions of each phase present.

To ensure that these properties/quantities are available when required, it is required that State Definition modules define the following quantities if they are not already one of the state variables chosen:

* `temperature` - the temperature of the mixture,
* `pressure` - the pressure of the mixture,
* `mole_frac_phase_comp` - mole fraction of the mixture by phase and component (even if only one phase is present),
* `phase_frac` - volume fractions of each phase (even if only one phase is present).

These quantities can be defined as either Pyomo `Vars` with associated `Constraints`, or as Pyomo `Expressions` as the developer desires. Developers may choose to include additional auxiliary variables as required by their needs (e.g. different forms of flow rates).

Supporting Constraints
^^^^^^^^^^^^^^^^^^^^^^

Depending upon the choice of state and auxiliary variables, developers may need to include a number of supporting constraints in their State Definitions. Common examples include constraints for the sum of mole fractions in the system, and relationships between different types of flow rates. Any number of constraints can be included by the developer to suit their needs, subject to the limitations of degrees of freedom.

However, developers need to be aware of the difference between inlet and outlet states and how this affects which constraints can be written. In the case of inlet states, all state variables are defined by the upstream process and thus no constraint can be written that involves only state variables (e.g. sum of mole fractions). For outlet (and intermediate) states however, it is often necessary to include these types of constraints to fully define the system. The IDAES Process Modeling Framework uses the `config.defined_state` configuration argument to indicate situations where the state variables should be considered fully defined (e.g. inlets) which can be used in `if` statements to determine whether a constraint should be included.

`always_flash`
^^^^^^^^^^^^^^

Whilst the set of state variables chosen must be sufficient for fully defining the state of the material, depending on the set of state variables chosen information of the phase separation (if applicable) may or may not be explicitly included. For example, using total flow rate and composition along with pressure and specific enthalpy is sufficient to define the state of the material, however it does not explicitly describe the phase fractions of the system. In these cases, it is necessary to perform a flash calculation at every state in the system to determine the phase fractions. However, If the state is defined in terms of flow rates by phase and component along with pressure and specific enthalpy, information on the phase separation is already included in the state definition and flash calculations are not required where the state is fully defined (i.e. `config.state_defined` is True).

To inform the Generic Property Package Framework of whether phase equilibrium calculations should be included when `config.state_defined` is True, all State Definitions are required to include a component named `always_flash` which is a boolean indicating whether equilibrium calculations should always be included (True) or only included when the state is not fully defined (False).

`get_material_flow_terms(phase, comp)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to automate the construction of the material balance equations, the IDAES Process Modeling Framework expects property packages to provide expressions for the flow terms in these equations. This is done via the `get_material_flow_terms` method which should return an expression involving variables in the StateBlock which should be used as the flow term in the material balances.

There are many forms this expression can take depending upon the state variables chosen and how the developer wishes to formulate the material balance equations, and the framework endeavors to support as many of these as possible. Material flow terms are defined on a phase-component basis (i.e. a separate expression for each component in each phase). An example of a `get_material_flow_term` using flow rate and mole fractions by phase is shown below.

.. code:: python

    def get_material_flow_terms_definition(phase, component):
        return self.flow_mol_phase[phase] * self.mole_frac_phase_comp[phase, component]
    self.get_material_flow_terms = get_material_flow_terms_definition

`get_enthalpy_flow_terms(phase)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the same way that `get_material_flow_terms` is used to automate construction of the material balance equations, automating the construction of the energy balance equations requires a `get_enthalpy_flow_terms` method. This method should return an expression for the enthalpy flow terms involving variables in the StateBlock.

There are many forms for the enthalpy flow terms as well, and developers may choose whichever best suits their needs. Enthalpy flow terms are defined on a phase basis, and an example is shown below using flow rate and specific enthalpy by phase.

.. code:: python

    def get_enthalpy_flow_terms_definition(phase):
        return self.flow_mol_phase[phase] * self.enth_mol_phase[phase]
    self.get_enthalpy_flow_terms = get_enthalpy_flow_terms_definition

`get_material_density_terms(phase, component)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For dynamic system, calculation of the material holdups also requires a material density term which is defined using the `get_material_density_terms` method. This method is defined in a similar fashion to the `get_material_flow_terms` method and is also defined on a phase-component basis.

`get_energy_density_terms(phase)`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For dynamic system, calculation of the energy holdups also requires an energy density term which is defined using the `get_energy_density_terms` method. This method is defined in a similar fashion to the `get_enthalpy_flow_terms` method and is also defined on a phase basis. Note however that the energy density term should only include internal energy contributions, and not the full enthalpy density (i.e. excluding the PV term).

`get_material_flow_basis()`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To automate generation of some terms in the balance equations, the IDAES Process Modeling Framework needs to know the basis (mass, mole or other) of the flow terms.  This is defined in the State Definition by providing a `get_material_flow_basis` method which returns a `MaterialFlowBasis` `Enum` (importable from `idaes.core`). E.g.:

.. code:: python

    def get_material_flow_basis_definition():
        return MaterialFlowBasis.molar
    self.get_material_flow_basis = get_material_flow_basis_definition

`default_material_balance_type()`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The IDAES Process Modeling Framework allows property packages to specify a default form for the material balance equations to be used if the modeler does not specify a form. Whilst not strictly required, developers are strongly encouraged to define a default form for the material balance equations.

To set the default material balance type, the State Definition must implement a method which returns a `MaterialBalanceType` `Enum` (importable from `idaes.core`. E.g.:

.. code:: python

    def default_material_balance_type_definition():
        return MaterialBalanceType.componentTotal
    self.default_material_balance_type = default_material_balance_type_definition

`default_energy_balance_type()`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The IDAES Process Modeling Framework allows property packages to specify a default form for the energy balance equations to be used if the modeler does not specify a form. Whilst not strictly required, developers are strongly encouraged to define a default form for the energy balance equations.

To set the default energy balance type, the State Definition must implement a method which returns an `EnergyBalanceType` `Enum` (importable from `idaes.core`. For an example, see `default_material_balance_type` above.

`define_port_members()`
^^^^^^^^^^^^^^^^^^^^^^^

In some situations, it is desirable to pass additional information between unit operations in a model beyond just the state variables. In these circumstance, the developer may define a `define_port_members` method which describes the information to be passed in `Ports` connecting units. This method should return a `dict` with a form similar to that of `define_state_vars`. Note that developers must also ensure that any additional information passed in `Ports` does not result in an over-specified problem, generally by excluding certain constraints in StateBlocks where `config.defined_state` is True.

If this method is not defined, `Ports` will default to using the variables described in `define_state_vars` instead.

`define_display_vars()`
^^^^^^^^^^^^^^^^^^^^^^^

Developers may also define a `define_display_vars` method which is used by the IDAES `report` methods to determine what information should be displayed for each state. The `define_display_vars` method should return a `dict` containing the information to display with the keys being the display name for the information and value being the quantity to display (similar to the `define_state_Vars` method). If this method is not defined then the `define_state_vars` method is used by the `report` methods instead.

`state_initialization(self)`
----------------------------

The `state_initialization` method is called as part of the Generic Property Package Framework `initialize` method and is expected to set initial guesses for any auxiliary variables defined by the State Definition based on the current values of the state variables. Note that the state variables will have been provided with initial guesses for the current state of the material from the process models, and thus will likely not be at their pre-defined initial conditions.

`self.do_not_initialize`
------------------------

The `do_not_initialize` component is a list containing a list of `Constraint` names which should remain deactivated during initialization of the StateBlock and only reactivated during the final step on initialization. Common examples of these are those constraints that are only written for outlet Blocks (i.e. those when `config.defined_state` is False), such as overall sum of mole fraction constraints.

