0D Control Volume Class
=======================

.. contents:: Contents
    :depth: 2

The ControlVolume0DBlock block is the most commonly used Control Volume class, and is used for
systems where there is a well-mixed volume of fluid, or where variations in spatial domains are
considered to be negligible. ControlVolume0DBlock blocks generally contain two
:ref:`StateBlocks <reference_guides/core/physical_property_class:State Blocks>` - one for the incoming material and one for
the material within and leaving the volume - and one
:ref:`StateBlocks <reference_guides/core/reaction_property_class:Reaction Blocks>`.

.. module:: idaes.core.base.control_volume0d

.. autoclass:: ControlVolume0DBlock
    :members:

.. autoclass:: ControlVolume0DBlockData
    :members:

ControlVolume0DBlock Equations
-------------------------------

This section documents the variables and constraints created by each of the methods provided by
the ControlVolume0DBlock class.

* :math:`t` indicates time index
* :math:`p` indicates phase index
* :math:`j` indicates component index
* :math:`e` indicates element index
* :math:`r` indicates reaction name index

add_geometry
^^^^^^^^^^^^

The add_geometry method creates a single variable within the control volume named `volume` indexed
by time (allowing for varying volume over time). A number of other methods depend on this variable
being present, thus this method should generally be called first.

**Variables**

.. csv-table::
   :header: "Variable Name", "Symbol", "Indices", "Conditions"
   :widths: 25, 10, 10, 30

   "volume", :math:`V_t`, "t", "None"

**Constraints**

No additional constraints

add_phase_component_balances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Material balances are written for each component in each phase (e.g. separate balances for liquid
water and steam). Physical property packages may include information to indicate that certain species
do not appear in all phases, and material balances will not be written in these cases (if
`has_holdup` is True holdup terms will still appear for these species, however these will be set to 0).

**Variables**

.. csv-table::
   :header: "Variable Name", "Symbol", "Indices", "Conditions"
   :widths: 25, 15, 10, 30

   "material_holdup", ":math:`M_{t,p,j}`", "t, p, j", "has_holdup = True"
   "phase_fraction", ":math:`\phi_{t,p}`", "t, p", "has_holdup = True"
   "material_accumulation", ":math:`\frac{\partial M_{t,p,j}}{\partial t}`", "t, p, j", "dynamic = True"
   "rate_reaction_generation", ":math:`N_{kinetic,t,p,j}`", "t, p ,j", "has_rate_reactions = True"
   "rate_reaction_extent", ":math:`X_{kinetic,t,r}`", "t, r", "has_rate_reactions = True"
   "equilibrium_reaction_generation", ":math:`N_{equilibrium,t,p,j}`", "t, p ,j", "has_equilibrium_reactions = True"
   "equilibrium_reaction_extent", ":math:`X_{equilibrium,t,r}`", "t, r", "has_equilibrium_reactions = True"
   "phase_equilibrium_generation", ":math:`N_{pe,t,p,j}`", "t, p ,j", "has_phase_equilibrium = True"
   "mass_transfer_term", ":math:`N_{transfer,t,p,j}`", "t, p ,j", "has_mass_transfer = True"


**Constraints**

`material_balances(t, p, j)`:

.. math:: \frac{\partial M_{t, p, j}}{\partial t} = F_{in, t, p, j} - F_{out, t, p, j} + N_{kinetic, t, p, j} + N_{equilibrium, t, p, j} + N_{pe, t, p, j} + N_{transfer, t, p, j} + N_{custom, t, p, j}

The :math:`N_{custom, t, p, j}` term allows the user to provide custom terms (variables or expressions) in both mass and molar basis which will be added into the material balances, which will be converted as necessary to the same basis as the material balance (by multiplying or dividing by the component molecular weight). The basis of the material balance is determined by the physical property package, and if undefined (or not mass or mole basis), an Exception will be returned.

If `has_holdup` is True, `material_holdup_calculation(t, p, j)`:

.. math:: M_{t, p, j} = \rho_{t, p, j} \times V_{t} \times \phi_{t, p}

where :math:`\rho_{t, p ,j}` is the density of component :math:`j` in phase :math:`p` at time :math:`t`

If `dynamic` is True:

Numerical discretization of the derivative terms, :math:`\frac{\partial M_{t,p,j}}{\partial t}`, will be performed by Pyomo.DAE.

If `has_rate_reactions` is True, `rate_reaction_stoichiometry_constraint(t, p, j)`:

.. math:: N_{kinetic, t, p, j} = \alpha_{r, p, j} \times X_{kinetic, t, r}

where :math:`\alpha_{r, p. j}` is the stoichiometric coefficient of component :math:`j` in phase :math:`p` for reaction :math:`r` (as defined in the PhysicalParameterBlock).

If `has_equilibrium_reactions` argument is True, `equilibrium_reaction_stoichiometry_constraint(t, p, j)`:

.. math:: N_{equilibrium, t, p, j} = \alpha_{r, p, j} \times X_{equilibrium, t, r}

where :math:`\alpha_{r, p. j}` is the stoichiometric coefficient of component :math:`j` in phase :math:`p` for reaction :math:`r` (as defined in the PhysicalParameterBlock).

add_total_component_balances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Material balances are written for each component across all phases (e.g. one balance for both liquid water and steam). Most terms in the balance equations are still indexed by both phase and component however. Physical property packages may include information to indicate that certain species do not appear in all phases, and material balances will not be written in these cases (if `has_holdup` is True holdup terms will still appear for these species, however these will be set to 0).

**Variables**

.. csv-table::
   :header: "Variable Name", "Symbol", "Indices", "Conditions"
   :widths: 25, 15, 10, 30

   "material_holdup", ":math:`M_{t,p,j}`", "t, p, j", "has_holdup = True"
   "phase_fraction", ":math:`\phi_{t,p}`", "t, p", "has_holdup = True"
   "material_accumulation", ":math:`\frac{\partial M_{t,p,j}}{\partial t}`", "t, p, j", "dynamic = True"
   "rate_reaction_generation", ":math:`N_{kinetic,t,p,j}`", "t, p ,j", "has_rate_reactions = True"
   "rate_reaction_extent", ":math:`X_{kinetic,t,r}`", "t, r", "has_rate_reactions = True"
   "equilibrium_reaction_generation", ":math:`N_{equilibrium,t,p,j}`", "t, p ,j", "has_equilibrium_reactions = True"
   "equilibrium_reaction_extent", ":math:`X_{equilibrium,t,r}`", "t, r", "has_equilibrium_reactions = True"
   "mass_transfer_term", ":math:`N_{transfer,t,p,j}`", "t, p ,j", "has_mass_transfer = True"


**Constraints**

`material_balances(t, j)`:

.. math:: \sum_p{\frac{\partial M_{t, p, j}}{\partial t}} = \sum_p{F_{in, t, p, j}} - \sum_p{F_{out, t, p, j}} + \sum_p{N_{kinetic, t, p, j}} + \sum_p{N_{equilibrium, t, p, j}} + \sum_p{N_{pe, t, p, j}} + \sum_p{N_{transfer, t, p, j}} + N_{custom, t, j}

The :math:`N_{custom, t, j}` term allows the user to provide custom terms (variables or expressions) in both mass and molar basis which will be added into the material balances, which will be converted as necessary to the same basis as the material balance (by multiplying or dividing by the component molecular weight). The basis of the material balance is determined by the physical property package, and if undefined (or not mass or mole basis), an Exception will be returned.

If `has_holdup` is True, `material_holdup_calculation(t, p, j)`:

.. math:: M_{t, p, j} = \rho_{t, p, j} \times V_{t} \times \phi_{t, p}

where :math:`\rho_{t, p ,j}` is the density of component :math:`j` in phase :math:`p` at time :math:`t`

If `dynamic` is True:

Numerical discretization of the derivative terms, :math:`\frac{\partial M_{t,p,j}}{\partial t}`, will be performed by Pyomo.DAE.

If `has_rate_reactions` is True,, `rate_reaction_stoichiometry_constraint(t, p, j)`:

.. math:: N_{kinetic, t, p, j} = \alpha_{r, p, j} \times X_{kinetic, t, r}

where :math:`\alpha_{r, p. j}` is the stoichiometric coefficient of component :math:`j` in phase :math:`p` for reaction :math:`r` (as defined in the PhysicalParameterBlock).

If `has_equilibrium_reactions` argument is True, `equilibrium_reaction_stoichiometry_constraint(t, p, j)`:

.. math:: N_{equilibrium, t, p, j} = \alpha_{r, p, j} \times X_{equilibrium, t, r}

where :math:`\alpha_{r, p. j}` is the stoichiometric coefficient of component :math:`j` in phase :math:`p` for reaction :math:`r` (as defined in the PhysicalParameterBlock).

add_total_element_balances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Material balances are written for each element in the mixture.

**Variables**

.. csv-table::
   :header: "Variable Name", "Symbol", "Indices", "Conditions"
   :widths: 25, 15, 10, 30

   "element_holdup", ":math:`M_{t,e}`", "t, e", "has_holdup = True"
   "phase_fraction", ":math:`\phi_{t,p}`", "t, p", "has_holdup = True"
   "element_accumulation", ":math:`\frac{\partial M_{t,e}}{\partial t}`", "t, e", "dynamic = True"
   "elemental_mass_transfer_term", ":math:`N_{transfer,t,e}`", "t, e", "has_mass_transfer = True"

**Expressions**

`elemental_flow_in(t, p, e)`:

.. math:: F_{in,t,p,e} = \sum_j{F_{in, t, p, j}} \times n_{j, e}

`elemental_flow_out(t, p, e)`:

.. math:: F_{out,t,p,e} = \sum_j{F_{out, t, p, j}} \times n_{j, e}

where :math:`n_{j, e}` is the number of moles of element :math:`e` in component :math:`j`.

**Constraints**

`element_balances(t, e)`:

.. math:: \frac{\partial M_{t, e}}{\partial t} = \sum_p{F_{in, t, p, e}} - \sum_p{F_{out, t, p, e}} + \sum_p{N_{transfer, t, e}} + N_{custom, t, e}

The :math:`N_{custom, t, e}` term allows the user to provide custom terms (variables or expressions) which will be added into the material balances.

If `has_holdup` is True, `elemental_holdup_calculation(t, e)`:

.. math:: M_{t, e} = V_{t} \times \sum_{p, j}{\phi_{t, p} \times \rho_{t, p, j} \times n_{j, e}}

where :math:`\rho_{t, p ,j}` is the density of component :math:`j` in phase :math:`p` at time :math:`t`

If `dynamic` is True:

Numerical discretization of the derivative terms, :math:`\frac{\partial M_{t,e}}{\partial t}`, will be performed by Pyomo.DAE.

add_total_enthalpy_balances
^^^^^^^^^^^^^^^^^^^^^^^^^^^

A single enthalpy balance is written for the entire mixture.

**Variables**

.. csv-table::
   :header: "Variable Name", "Symbol", "Indices", "Conditions"
   :widths: 25, 15, 10, 30

   "energy_holdup", ":math:`E_{t,p}`", "t, p", "has_holdup = True"
   "phase_fraction", ":math:`\phi_{t,p}`", "t, p", "has_holdup = True"
   "energy_accumulation", ":math:`\frac{\partial E_{t,p}}{\partial t}`", "t, p", "dynamic = True"
   "heat", ":math:`Q_{t}`", "t", "has_heat_transfer = True"
   "work", ":math:`W_{t}`", "t", "has_work_transfer = True"
   "enthalpy_transfer", ":math:`H_{transfer,t}`", "t", "has_enthalpy_transfer = True"

**Expressions**

`heat_of_reaction(t)`:

.. math:: Q_{rxn, t} = sum_r{X_{kinetic, t, r} \times \Delta H_{rxn, r}} + sum_r{X_{equilibrium, t, r} \times \Delta H_{rxn, r}}

where :math:`Q_{rxn, t}` is the total enthalpy released by both kinetic and equilibrium reactions, and :math:`\Delta H_{rxn, r}` is the specific heat of reaction for reaction :math:`r`.

**Parameters**

.. csv-table::
   :header: "Parameter Name", "Symbol", "Default Value"
   :widths: 25, 15, 15

   "scaling_factor_energy", ":math:`s_{energy}`", "1E-6"

**Constraints**

`enthalpy_balance(t)`:

.. math:: s_{energy} \times \sum_p{\frac{\partial E_{t, p}}{\partial t}} = s_{energy} \times \sum_p{H_{in, t, p}} - s_{energy} \times \sum_p{H_{out, t, p}} + s_{energy} \times Q_t + s_{energy} \times W_t + + s_{energy} \times H_{transfer,t} +s_{energy} \times Q_{rxn, t} + s_{energy} \times E_{custom, t}

The :math:`E_{custom, t}` term allows the user to provide custom terms  which will be added into the energy balance.

If `has_holdup` is True, `enthalpy_holdup_calculation(t, p)`:

.. math:: E_{t, p} = u_{t, p} \times V_{t} \times \phi_{t, p}

where :math:`u_{t, p}` is the internal energy density (specific internal energy) of phase :math:`p` at time :math:`t`

If `dynamic` is True:

Numerical discretization of the derivative terms, :math:`\frac{\partial E_{t,p}}{\partial t}`, will be performed by Pyomo.DAE.

add_total_pressure_balances
^^^^^^^^^^^^^^^^^^^^^^^^^^^

A single pressure balance is written for the entire mixture.

**Variables**

.. csv-table::
   :header: "Variable Name", "Symbol", "Indices", "Conditions"
   :widths: 25, 15, 10, 30

   "deltaP", ":math:`\Delta P_{t}`", "t", "has_pressure_change = True"

**Parameters**

.. csv-table::
   :header: "Parameter Name", "Symbol", "Default Value"
   :widths: 25, 15, 15

   "scaling_factor_pressure", ":math:`s_{pressure}`", "1E-4"

**Constraints**

`pressure_balance(t)`:

.. math:: 0 = s_{pressure} \times P_{in, t} - s_{pressure} \times P_{out, t} + s_{pressure} \times \Delta P_t + s_{pressure} \times \Delta P_{custom, t}

The :math:`\Delta P_{custom, t}` term allows the user to provide custom terms  which will be added into the pressure balance.
