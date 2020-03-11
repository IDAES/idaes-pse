Moving Bed Reactor
====================================

The IDAES Moving Bed Reactor (MBR) model represents a unit operation where two material streams – one phase is solid while the other is gas or liquid – pass through a linear reactor vessel while undergoing some chemical reaction(s). The two streams have opposite flow direction (counter-flow). This model requires modeling the system in one spatial dimension.

Degrees of Freedom
------------------

MBRs generally have at least 2 (or more) degrees of freedom. Typically fixed variables are reactor length, and diameter. 


Model Structure
---------------

The core MBR unit model consists of two ControlVolume1DBlock Blocks (named gas_phase and solid_phase), each with one Inlet Port (named gas_inlet and solid_inlet) and one Outlet Port (named gas_outlet and solid_outlet).

Construction Arguments
----------------------

MBR unit models have construction arguments passed to the Holdup Block specific to the unit model, the gas phase and the solid phase. Arguments specific to the unit model are: 

* dynamic - indicates whether this model will be dynamic or not (default = useDefault). Valid values are:

		============== ========================================
		Argument       Description
		============== ========================================
		useDefault     get flag from parent (default = False)
		True           set as a dynamic model
		False          set as a steady-state model 
		============== ========================================

* has_holdup -  argument to indicate whether holdup terms should be constructed or not (default = False). Must be True if dynamic = True.
* finite_elements - sets the number of finite elements to use when discretizing the spatial domains (default = 10). This is used for both gas and solid phase domains.
* length_domain_set - (optional) list of point to use to initialize a new ContinuousSet if length_domain is not provided (default = [0.0, 1.0]).
* transformation_method - argument to specify the DAE transformation method for the domain (default = “dae.finite_difference”); should be compatible with the Pyomo DAE TransformationFactory.
* transformation_scheme - argument to specify the scheme to use for the selected DAE transformation method (default = “BACKWARD”); should be compatible with the Pyomo DAE TransformationFactory. See Pyomo documentation for supported schemes.
* collocation_points - sets the number of collocation points to use when discretizing the spatial domains (default = 3, collocation methods only).
* flow_type - flow configuration of Moving Bed (default="counter_current", gas side flows from 0 to 1, solid side flows from 1 to 0).
* material_balance_type - indicates what type of mass balance should be constructed (default = MaterialBalanceType.componentTotal). All valid values are:

		==================================== ==============================
		rgument       						 Description
		==================================== ==============================
		MaterialBalanceType.none     		 exclude material balances
		MaterialBalanceType.componentPhase   use phase component balances
		MaterialBalanceType.componentTotal   use total component balances
		MaterialBalanceType.elementTotal     use total element balances
		MaterialBalanceType.total            use total material balance 
		==================================== ==============================

* energy_balance_type - indicates what type of energy balance should be constructed (default = EnergyBalanceType.enthalpyTotal). All valid values are:

		==================================== =========================================
		Argument       						 Description
		==================================== =========================================
		EnergyBalanceType.none				 exclude energy balances
		EnergyBalanceType.enthalpyTotal		 single enthalpy balance for material
		EnergyBalanceType.enthalpyPhase		 enthalpy balances for each phase
		EnergyBalanceType.energyTotal		 single energy balance for material
		EnergyBalanceType.energyPhase		 energy balances for each phase
		==================================== =========================================

* area_definition - argument defining whether area variable should be vary across spatial domian or not (default = DistributedVars.variant).

Arguments specific to the gas and solid phases are:

* momentum_balance_type - indicates what type of momentum balance should be constructed (default = MomentumBalanceType.none). Valid values are:

		==================================== =========================================
		Argument       						 Description
		==================================== =========================================
		MomentumBalanceType.none			 exclude momentum balances
		MomentumBalanceType.pressureTotal	 single pressure balance for material
		MomentumBalanceType.pressurePhase	 pressure balances for each phase
		MomentumBalanceType.momentumTotal	 single momentum balance for material
		MomentumBalanceType.momentumPhase	 momentum balances for each phase
		==================================== =========================================

* has_pressure_change - indicates whether terms for pressure change should be constructed (default = False).
* pressure_drop_type - indicates what type of pressure drop correlation should be used (default = "simple_correlation"). Valid values are:

		=========================== =========================================
		Argument       				Description
		=========================== =========================================
		"simple_correlation"		use a simplified pressure drop correlation
		"ergun_correlation"			use the Ergun equation 
		=========================== =========================================
		
  For the gas phase, the momentum balance is set to: CONFIG.gas_phase_config.momentum_balance_type = MomentumBalanceType.pressureTotal.

* has_phase_equilibrium - argument to enable phase equilibrium on the gas side (default = False).
* has_equilibrium_reactions - indicates whether terms for equilibrium controlled reactions should be constructed (default = False).
* property_package - property package to use when constructing Property Blocks (default = ‘use_parent_value’). This is provided as a Property Parameter Block by the Flowsheet when creating the model. If a value is not provided, the Holdup Block will try to use the default property package if one is defined.
* property_package_args - set of arguments to be passed to the Property Blocks when they are created.
* reaction_package - reaction package to use when constructing Property Blocks (default = ‘use_parent_value’). This is provided as a Property Parameter Block by the Flowsheet when creating the model. If a value is not provided, the Holdup Block will try to use the default property package if one is defined.
* reaction_package_args - set of arguments to be passed to the Property Blocks when they are created.

Variables
---------

PFR units add the following additional Variables:

====================== ======= ===============================================================
Variable               Name    Notes
====================== ======= ===============================================================
:math:`L`              length  Reference to control_volume.length
:math:`A`              area    Reference to control_volume.area
:math:`V`              volume  Reference to control_volume.volume
:math:`Q_{t,x}`        heat    Only if has_heat_transfer = True, reference to holdup.heat
:math:`\Delta P_{t,x}` deltaP  Only if has_pressure_change = True, reference to holdup.deltaP
====================== ======= ===============================================================

Constraints
-----------

MBR units write the following additional Constraints at all points in the spatial domain:

.. math:: X_{t,x,r} = A \times r_{t,x,r}

where :math:`X_{t,x,r}` is the extent of reaction of reaction :math:`r` at point :math:`x` and time :math:`t`, :math:`A` is the cross-sectional area of the reactor and :math:`r_{t,r}` is the volumetric rate of reaction of reaction :math:`r` at point :math:`x` and time :math:`t` (from the outlet StateBlock).

Addtional Constraints
---------------------

The MBR model writes the following additional Constraints to describe the correlation used to compute the pressure drop in the reactor. 

If `pressure_drop_type` is `simple_correlation`:

.. math:: - \frac{\delta P - \delta z}{} = \left( \rho_{p} - \rho_{g} \right( \times a_{E} \times u_{g}

where :math:`P` is the system pressure, :math:`z` is the spatial domain, :math:`\rho_{p}` is the density of the particles, :math:`\rho_{g}` is the density of the gas, :math:`u_{g}` is the superficial velocity of the gas, and :math:`a_{E}` is a parameter having the value of 0.2. 

If `pressure_drop_type` is `ergun_correlation`:

.. math:: - \frac{\delta P - \delta z}{} = \frac{150 \times \mu_{g} \times {\left 1 - \epsilon \right}^{2}}{\epsilon^{3} \times {d_{p}}^2} + \frac{1.75 \times \left 1 - \epsilon \right \rho_{g} {\left u_{g} + u_{s} \right}^2}{\epsilon^{3} \times d_{p}}

where :math:`\mu_{g}` is the gas viscosity, :math:`\epsilon` is the reactor bed voidage, :math:`u_{s}` is the velocity of the solids, and :math:`d_{p}` is the diameter of the solid particles.

MBR Class
---------

.. module:: tasks.task-3.task-3.2.chemical_looping.unit_models.moving_bed

.. autoclass:: MBR
    :members:

MBRData Class
-------------

.. autoclass:: MBRData
    :members:
