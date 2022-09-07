1D Control Volume Class
=======================

.. contents:: Contents
    :depth: 2

The ControlVolume1DBlock block is used for systems with one spatial dimension where material flows parallel to the spatial domain. Examples of these types of unit operations include plug flow reactors and pipes. ControlVolume1DBlock blocks are discretized along the length domain and contain one StateBlock and one ReactionBlock (if applicable) at each point in the domain (including the inlet and outlet).

.. module:: idaes.core.base.control_volume1d

.. autoclass:: ControlVolume1DBlock
    :members:

.. autoclass:: ControlVolume1DBlockData
    :members:

ControlVolume1DBlock Equations
-------------------------------

This section documents the variables and constraints created by each of the methods provided by the ControlVolume0DBlock class.

* :math:`t` indicates time index
* :math:`x` indicates spatial (length) index
* :math:`p` indicates phase index
* :math:`j` indicates component index
* :math:`e` indicates element index
* :math:`r` indicates reaction name index

Most terms within the balance equations written by ControlVolume1DBlock are on a basis of per unit length (e.g. :math:`mol/m \cdot s`).

add_geometry
^^^^^^^^^^^^

The add_geometry method creates the normalized length domain for the control volume (or a reference to an external domain). All constraints in ControlVolume1DBlock assume a normalized length domain, with values between 0 and 1.

This method also adds variables and constraints to describe the geometry of the control volume. ControlVolume1DBlock does not support varying dimensions of the control volume with time at this stage.

**Variables**

.. csv-table::
   :header: "Variable Name", "Symbol", "Indices", "Conditions"
   :widths: 25, 10, 10, 30

   "length_domain", ":math:`x`", "None", "None"
   "volume", ":math:`V`", "None", "None"
   "area", ":math:`A`", "None", "None"
   "length", ":math:`L`", "None", "If `length_var` argument is provided, a reference to the provided component is made in place of creating a new variable"

**Constraints**

`geometry_constraint`:

.. math:: V = A \times L

add_phase_component_balances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Material balances are written for each component in each phase (e.g. separate balances for liquid water and steam). Physical property packages may include information to indicate that certain species do not appear in all phases, and material balances will not be written in these cases (if `has_holdup` is True holdup terms will still appear for these species, however these will be set to 0).

**Variables**

.. csv-table::
   :header: "Variable Name", "Symbol", "Indices", "Conditions"
   :widths: 25, 15, 10, 30

   "material_holdup", ":math:`M_{t,x,p,j}`", "t, x, p, j", "has_holdup = True"
   "phase_fraction", ":math:`\phi_{t,x,p}`", "t, x, p", "has_holdup = True"
   "material_accumulation", ":math:`\frac{\partial M_{t,x,p,j}}{\partial t}`", "t, x, p, j", "dynamic = True"
   "_flow_terms", ":math:`F_{t, x, p, j}`", "t, x, p, j", None
   "material_flow_dx", ":math:`\frac{\partial F_{t,x,p,j}}{\partial x}`", "t, x, p, j", None
   "rate_reaction_generation", ":math:`N_{kinetic,t,x,p,j}`", "t, x, p ,j", "has_rate_reactions = True"
   "rate_reaction_extent", ":math:`X_{kinetic,t,x,r}`", "t, x, r", "has_rate_reactions = True"
   "equilibrium_reaction_generation", ":math:`N_{equilibrium,t,x,p,j}`", "t, x, p ,j", "has_equilibrium_reactions = True"
   "equilibrium_reaction_extent", ":math:`X_{equilibrium,t,x,r}`", "t, x, r", "has_equilibrium_reactions = True"
   "phase_equilibrium_generation", ":math:`N_{pe,t,x,p,j}`", "t, x, p ,j", "has_phase_equilibrium = True"
   "mass_transfer_term", ":math:`N_{transfer,t,x,p,j}`", "t, x, p ,j", "has_mass_transfer = True"


**Constraints**

`material_balances(t, x, p, j)`:

.. math:: L \times \frac{\partial M_{t, x, p, j}}{\partial t} = fd \times \frac{\partial F_{t, x, p, j}}{\partial x} + L \times N_{kinetic, t, x, p, j} + L \times N_{equilibrium, t, x, p, j} + L \times N_{pe, t, x, p, j} + L \times N_{transfer, t, x, p, j} + L \times N_{custom, t, x, p, j}

:math:`fd` is a flow direction term, which allows for material flow to be defined in either direction. If material flow is defined as `forward`, :math:`fd = -1`, otherwise :math:`fd = 1`.

The :math:`N_{custom, t, x, p, j}` term allows the user to provide custom terms (variables or expressions) in both mass and molar basis which will be added into the material balances, which will be converted as necessary to the same basis as the material balance (by multiplying or dividing by the component molecular weight). The basis of the material balance is determined by the physical property package, and if undefined (or not mass or mole basis), an Exception will be returned.

`material_flow_linking_constraints(t, x, p, j)`:

This constraint is an internal constraint used to link the extensive material flow terms in the StateBlocks into a single indexed variable. This is required as Pyomo.DAE requires a single indexed variable to create the associated DerivativeVars and their numerical expansions.

If `has_holdup` is True, `material_holdup_calculation(t, x, p, j)`:

.. math:: M_{t, x, p, j} = \rho_{t, x, p, j} \times A \times \phi_{t, x, p}

where :math:`\rho_{t, x, p ,j}` is the density of component :math:`j` in phase :math:`p` at time :math:`t` and location :math:`x`.

If `dynamic` is True:

Numerical discretization of the derivative terms, :math:`\frac{\partial M_{t,x,p,j}}{\partial t}`, will be performed by Pyomo.DAE.

If `has_rate_reactions` is True, `rate_reaction_stoichiometry_constraint(t, x, p, j)`:

.. math:: N_{kinetic, t, x, p, j} = \alpha_{r, p, j} \times X_{kinetic, t, x, r}

where :math:`\alpha_{r, p. j}` is the stoichiometric coefficient of component :math:`j` in phase :math:`p` for reaction :math:`r` (as defined in the PhysicalParameterBlock).

If `has_equilibrium_reactions` argument is True, `equilibrium_reaction_stoichiometry_constraint(t, x, p, j)`:

.. math:: N_{equilibrium, t, x, p, j} = \alpha_{r, p, j} \times X_{equilibrium, t, x, r}

where :math:`\alpha_{r, p. j}` is the stoichiometric coefficient of component :math:`j` in phase :math:`p` for reaction :math:`r` (as defined in the PhysicalParameterBlock).

add_total_component_balances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Material balances are written for each component across all phases (e.g. one balance for both liquid water and steam). Physical property packages may include information to indicate that certain species do not appear in all phases, and material balances will not be written in these cases (if `has_holdup` is True holdup terms will still appear for these species, however these will be set to 0).

**Variables**

.. csv-table::
   :header: "Variable Name", "Symbol", "Indices", "Conditions"
   :widths: 25, 15, 10, 30

   "material_holdup", ":math:`M_{t,x,p,j}`", "t, x, p, j", "has_holdup = True"
   "phase_fraction", ":math:`\phi_{t,x,p}`", "t, x, p", "has_holdup = True"
   "material_accumulation", ":math:`\frac{\partial M_{t,x,p,j}}{\partial t}`", "t, x, p, j", "dynamic = True"
   "_flow_terms", ":math:`F_{t, x, p, j}`", "t, x, p, j", None
   "material_flow_dx", ":math:`\frac{\partial F_{t,x,p,j}}{\partial x}`", "t, x, p, j", None
   "rate_reaction_generation", ":math:`N_{kinetic,t,x,p,j}`", "t, x, p ,j", "has_rate_reactions = True"
   "rate_reaction_extent", ":math:`X_{kinetic,t,x,r}`", "t, x, r", "has_rate_reactions = True"
   "equilibrium_reaction_generation", ":math:`N_{equilibrium,t,x,p,j}`", "t, x, p ,j", "has_equilibrium_reactions = True"
   "equilibrium_reaction_extent", ":math:`X_{equilibrium,t,x,r}`", "t, x, r", "has_equilibrium_reactions = True"
   "mass_transfer_term", ":math:`N_{transfer,t,x,p,j}`", "t, x, p ,j", "has_mass_transfer = True"


**Constraints**

`material_balances(t, x, p, j)`:

.. math:: L \times \sum_p{\frac{\partial M_{t, x, p, j}}{\partial t}} = fd \times \sum{\frac{\partial F_{t, x, p, j}}{\partial x}} + L \times \sum_p{N_{kinetic, t, x, p, j}} + L \times \sum_p{N_{equilibrium, t, x, p, j}} + L \times \sum_p{N_{transfer, t, x, p, j}} + L \times N_{custom, t, x, j}

:math:`fd` is a flow direction term, which allows for material flow to be defined in either direction. If material flow is defined as `forward`, :math:`fd = -1`, otherwise :math:`fd = 1`.

The :math:`N_{custom, t, x, j}` term allows the user to provide custom terms (variables or expressions) in both mass and molar basis which will be added into the material balances, which will be converted as necessary to the same basis as the material balance (by multiplying or dividing by the component molecular weight). The basis of the material balance is determined by the physical property package, and if undefined (or not mass or mole basis), an Exception will be returned.

`material_flow_linking_constraints(t, x, p, j)`:

This constraint is an internal constraint used to link the extensive material flow terms in the StateBlocks into a single indexed variable. This is required as Pyomo.DAE requires a single indexed variable to create the associated DerivativeVars and their numerical expansions.

If `has_holdup` is True, `material_holdup_calculation(t, x, p, j)`:

.. math:: M_{t, x, p, j} = \rho_{t, x, p, j} \times A \times \phi_{t, x, p}

where :math:`\rho_{t, x, p ,j}` is the density of component :math:`j` in phase :math:`p` at time :math:`t` and location :math:`x`.

If `dynamic` is True:

Numerical discretization of the derivative terms, :math:`\frac{\partial M_{t,x,p,j}}{\partial t}`, will be performed by Pyomo.DAE.

If `has_rate_reactions` is True, `rate_reaction_stoichiometry_constraint(t, x, p, j)`:

.. math:: N_{kinetic, t, x, p, j} = \alpha_{r, p, j} \times X_{kinetic, t, x, r}

where :math:`\alpha_{r, p. j}` is the stoichiometric coefficient of component :math:`j` in phase :math:`p` for reaction :math:`r` (as defined in the PhysicalParameterBlock).

If `has_equilibrium_reactions` argument is True, `equilibrium_reaction_stoichiometry_constraint(t, x, p, j)`:

.. math:: N_{equilibrium, t, x, p, j} = \alpha_{r, p, j} \times X_{equilibrium, t, x, r}

where :math:`\alpha_{r, p. j}` is the stoichiometric coefficient of component :math:`j` in phase :math:`p` for reaction :math:`r` (as defined in the PhysicalParameterBlock).

add_total_element_balances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Material balances are written for each element in the mixture.

**Variables**

.. csv-table::
   :header: "Variable Name", "Symbol", "Indices", "Conditions"
   :widths: 25, 15, 10, 30

   "element_holdup", ":math:`M_{t,x,e}`", "t, x, e", "has_holdup = True"
   "phase_fraction", ":math:`\phi_{t,x,p}`", "t, x, p", "has_holdup = True"
   "element_accumulation", ":math:`\frac{\partial M_{t,x,e}}{\partial t}`", "t, x, e", "dynamic = True"
   "elemental_mass_transfer_term", ":math:`N_{transfer,t,x,e}`", "t, x, e", "has_mass_transfer = True"
   "elemental_flow_term", ":math:`F_{t,x,e}`", "t, x, e", None


**Constraints**

`elemental_flow_constraint(t, x, e)`:

.. math:: F_{t,x,e} = \sum_p{\sum_j{F_{t,x,p,j} \times n_{j, e}}}

where :math:`n_{j, e}` is the number of moles of element :math:`e` in component :math:`j`.

`element_balances(t, x, e)`:

.. math:: L \times \frac{\partial M_{t, x, e}}{\partial t} = fd \times \frac{\partial F_{t, x, e}}{\partial x} + L \times N_{transfer, t, p, j} + L \times N_{custom, t, e}

:math:`fd` is a flow direction term, which allows for material flow to be defined in either direction. If material flow is defined as `forward`, :math:`fd = -1`, otherwise :math:`fd = 1`.

The :math:`N_{custom, t, x, e}` term allows the user to provide custom terms (variables or expressions) which will be added into the material balances.

If `has_holdup` is True, `elemental_holdup_calculation(t, x, e)`:

.. math:: M_{t, x, e} = \rho_{t, x, p, j} \times A \times \phi_{t, x, p}

where :math:`\rho_{t, x, p ,j}` is the density of component :math:`j` in phase :math:`p` at time :math:`t` and location :math:`x`.

If `dynamic` is True:

Numerical discretization of the derivative terms, :math:`\frac{\partial M_{t,x,p,j}}{\partial t}`, will be performed by Pyomo.DAE.

add_total_enthalpy_balances
^^^^^^^^^^^^^^^^^^^^^^^^^^^

A single enthalpy balance is written for the entire mixture at each point in the spatial domain.

**Variables**

.. csv-table::
   :header: "Variable Name", "Symbol", "Indices", "Conditions"
   :widths: 25, 15, 10, 30

   "energy_holdup", ":math:`E_{t,x,p}`", "t, x, p", "has_holdup = True"
   "phase_fraction", ":math:`\phi_{t,x,p}`", "t, x, p", "has_holdup = True"
   "energy_accumulation", ":math:`\frac{\partial E_{t,x,p}}{\partial t}`", "t, x, p", "dynamic = True"
   "_enthalpy_flow", ":math:`H_{t,x,p}`", "t, x, p", "None"
   "enthalpy_flow_dx", ":math:`\frac{\partial H_{t,x,p}}{\partial x}`", "t, x, p", "None"
   "heat", ":math:`Q_{t,x}`", "t, x", "has_heat_transfer = True"
   "work", ":math:`W_{t,x}`", "t, x", "has_work_transfer = True"
   "enthalpy_transfer", ":math:`H_{transfer,t,x}`", "t, x", "has_enthalpy_transfer = True"

**Expressions**

`heat_of_reaction(t, x)`:

.. math:: Q_{rxn, t, x} = sum_r{X_{kinetic, t, x, r} \times \Delta H_{rxn, r}} + sum_r{X_{equilibrium, t, x, r} \times \Delta H_{rxn, r}}

where :math:`Q_{rxn, t, x}` is the total enthalpy released by both kinetic and equilibrium reactions, and :math:`\Delta H_{rxn, r}` is the specific heat of reaction for reaction :math:`r`.

**Parameters**

.. csv-table::
   :header: "Parameter Name", "Symbol", "Default Value"
   :widths: 25, 15, 15

   "scaling_factor_energy", ":math:`s_{energy}`", "1E-6"

**Constraints**

`enthalpy_balance(t)`:

.. math:: s_{energy} \times L \times \sum_p{\frac{\partial E_{t, x, p}}{\partial t}} = s_{energy} \times fd \ times \sum_p{\frac{\partial H_{t, x, p}}{\partial x}} + s_{energy} \times L \times Q_{t,x} + s_{energy} \times L \times W_{t,x} + s_{energy} \times L \times H_{transfer,t,x} + s_{energy} \times L \times Q_{rxn, t, x} + s_{energy} \times L \times E_{custom, t, x}

:math:`fd` is a flow direction term, which allows for material flow to be defined in either direction. If material flow is defined as `forward`, :math:`fd = -1`, otherwise :math:`fd = 1`.

The :math:`E_{custom, t, x}` term allows the user to provide custom terms  which will be added into the energy balance.

`enthalpy_flow_linking_constraints(t, x, p)`:

This constraint is an internal constraint used to link the extensive enthalpy flow terms in the StateBlocks into a single indexed variable. This is required as Pyomo.DAE requires a single indexed variable to create the associated DerivativeVars and their numerical expansions.

If `has_holdup` is True, `enthalpy_holdup_calculation(t, x, p)`:

.. math:: E_{t, x, p} = u_{t, x, p} \times A \times \phi_{t, x, p}

where :math:`u_{t, x, p}` is the internal density (specific internal energy) of phase :math:`p` at time :math:`t` and location :math:`x`.

If `dynamic` is True:

Numerical discretization of the derivative terms, :math:`\frac{\partial E_{t,x,p}}{\partial t}`, will be performed by Pyomo.DAE.

add_total_pressure_balances
^^^^^^^^^^^^^^^^^^^^^^^^^^^

A single pressure balance is written for the entire mixture at all points in the spatial domain.

**Variables**

.. csv-table::
   :header: "Variable Name", "Symbol", "Indices", "Conditions"
   :widths: 25, 15, 10, 30

   "pressure", ":math:`P_{t,x}`", "t, x", None
   "pressure_dx", ":math:`\frac{\partial P_{t,x}}{\partial x}`", "t, x", "None"
   "deltaP", ":math:`\Delta P_{t,x}`", "t, x", "has_pressure_change = True"

**Parameters**

.. csv-table::
   :header: "Parameter Name", "Symbol", "Default Value"
   :widths: 25, 15, 15

   "scaling_factor_pressure", ":math:`s_{pressure}`", "1E-4"

**Constraints**

`pressure_balance(t, x)`:

.. math:: 0 = s_{pressure} \times fd \times \frac{\partial P_{t,x}}{\partial x} + s_{pressure} \times L \times \Delta P_{t,x} + s_{pressure} \times L \times \Delta P_{custom, t, x}

:math:`fd` is a flow direction term, which allows for material flow to be defined in either direction. If material flow is defined as `forward`, :math:`fd = -1`, otherwise :math:`fd = 1`.

The :math:`\Delta P_{custom, t, x}` term allows the user to provide custom terms  which will be added into the pressure balance.

`pressure_linking_constraint(t, x)`:

This constraint is an internal constraint used to link the pressure terms in the StateBlocks into a single indexed variable. This is required as Pyomo.DAE requires a single indexed variable to create the associated DerivativeVars and their numerical expansions.
