Cubic Equations of State
========================

This property package implements a general form of a cubic equation of state which can be used for most cubibc-type equations of state. This package supports phase equilibrium calucations with a smooth phase transition formulation that makes it amenable for equation oriented optimization. The following equations of state are currently supported:

* Peng-Robinson
* Soave-Redlich-Kwong

**Flow basis**: Molar

**Units**: SI units

**State Variables**: 

The state block uses the following state variables:

.. raw:: html
	<head>
   	</head>
	
   	<body>
      	<ul>
 	<li>Total molar flow rate (mol/s) - <code> <font color="red"> flow_mol </font> </code>
 	<li>Temperature (K) - <code> <font color="red"> temperature </font> </code>
	<li>Presure (Pa) - <code> <font color="red"> pressure</font></code>
	<li>Mole fraction of the mixture - <code> <font color="red"> mole_frac_comp</font></code> 
	</ul>
	</body>

Inputs
------

When instantiating the parameter block that uses this particular state block, 1 optional argument can be passed:

.. raw:: html
	<head>

   	</head>
   	<body>
      	<ul>
 	<li><code> <font color="red"> valid_phase </font> </code> - <font color="green">"Liq"</font>  or <font color="green">"Vap"</font> or <font color="green">("Liq", "Vap") </font> or <font color="green">("Vap", "Liq") </font> </li></li>
	</ul>
	</body>

The ``valid_phase`` argument denotes the valid phases for a given set of inlet conditions. For example, if the user knows a priori that the it will only be a single phase (for example liquid only), then it is best not to include the complex flash equilibrium constraints in the model. If the user does not specify any option, then the package defaults to a 2 phase assumption meaning that the constraints to compute the phase equilibrium will be computed.

Degrees of Freedom
------------------
In general, the general cubic equation of state has a number of degrees of freedom equal to 2 + the number of components in the system (total flow rate, temperature, pressure and N-1 mole fractions). In soem cases (primarily inlets to units), this is increased by 1 due to the removal of a constraint on the sum of mole fractions.

VLE Model with Smooth Phase Transition
--------------------------------------

The flash equations consists of the following equations depending on the choice of state variables selected by the user. 

If the state variables are total flow, mole fraction, temperature, and pressure, then the following constraints are implemented:

.. math:: F^{in} = F^{liq} + F^{vap}
.. math:: z_{i}^{in}F^{in} = x_{i}^{liq}F^{liq} + y_{i}^{vap}F^{vap}

If the state variables are component flow rates, temperature, and pressure, then the following constraints are implemented:

.. math:: F^{in}_{i} = F^{liq}_{i} + F^{vap}_{i}

The equilibrium condition, the fugacity of the vapor and liquid phase are defined as follows:

.. math:: f_{i}^{vap} = f_{i}^{liq}
.. math:: f_{i}^{vap} = y_{i}\phi_{i}P
.. math:: f_{i}^{liq} = x_{i}p^{sat}_{i}\nu_{i}

The equilibrium constraint is written as a generic constraint such that it can be extended easily for non-ideal gases and liquids. As this property package only supports an ideal gas, the fugacity coefficient (:math:`\phi_{i}`) for the vapor phase is 1 and hence the expression reduces to :math:`y_{i}P`.

Typically, the flash calculations are computed at a given temperature, :math:`T`. However, the flash calculations become trivial if the given conditions do not fall in the two phase region. For simulation only studies, the user may know a priori the condition of the stream but when the same set of equations are used for optimization, there is a high probablity that the specifications can transcend the phase envelope and hence the flash equations included may be trivial in the single phase region (i.e. liquid or vapor only). To circumvent this problem, property packages in IDAES that support VLE will compute the flash calculations at an "equilibrium" temperature :math:`T_{eq}`. The equilibrium temperature is computed as follows:

.. math:: T_{1} = max(T_{bubble}, T) 
.. math:: T_{eq} = min(T_{1}, T_{dew})

where :math:`T_{eq}` is the equilibrium temperature at which flash calculations are computed, :math:`T` is the stream temperature, :math:`T_{1}` is the intermediate temperature variable, :math:`T_{bubble}` is the bubble point temperature of mixture, and :math:`T_{dew}` is the dew point temperature of the mixture. Note that, in the above equations, approximations are used for the max and min functions as follows:

.. math:: T_{1} = 0.5{[T + T_{bubble} + \sqrt{(T-T_{bubble})^2 + \epsilon_{1}^2}]}
.. math:: T_{eq} = 0.5{[T_{1} + T_{dew} - \sqrt{(T-T_{dew})^2 + \epsilon_{2}^2}]}

where :math:`\epsilon_1` and :math:`\epsilon_2` are smoothing parameters (mutable). The default values are 0.01 and 0.0005 respectively. It is recommended that :math:`\epsilon_1` > :math:`\epsilon_2`. Please refer to reference 4 for more details. Therefore, it can be seen that if the stream temperature is less than that of the bubble point temperature, the VLE calucalations will be computed at the bubble point. Similarly, if the stream temperature is greater than the dew point temperature, then the VLE calculations are computed at the dew point temperature. For all other conditions, the equilibrium calcualtions will be computed at the actual temperature. 

Additional constraints are included in the model to compute the thermodynamic properties based on the cubci equation of state, such as compressibility factors, fugacities, enthalpies and entropies. All of the equations come from "THe Properties of Gases and Liquids, 4th Edition" by Reid, Prausnitz and Poling. Please note that, these constraints are added only if the variable is called for when building the model. This eliminates adding unnecessary constraints to compute properties that are not needed in the model. 

List of Variables
-----------------
.. csv-table::
   :header: "Variable Name", "Description", "Units"

   "``flow_mol``", "Total molar flow rate", "mol/s"
   "``mole_frac_comp``", "Mixture mole fraction indexed by component", "None"
   "``temperature``", "Temperature", "K"
   "``pressure``", "Pressure", "Pa"
   "``flow_mol_phase``", "Molar flow rate indexed by phase", "mol/s"
   "``mole_frac_phase_comp``", "Mole fraction indexed by phase and component", "None"
   "``pressure_sat``", "Saturation or vapor pressure indexed by component", "Pa"
   "``dens_mol_phase``", "Molar density indexed by phase", "mol/m3"
   "``dens_mass_phase``", "Mass density indexed by phase", "kg/m3"
   "``enth_mol_phase``", "Molar enthalpy indexed by phase ", "J/mol"
   "``enth_mol``", "Molar enthalpy of mixture", "J/mol"
   "``entr_mol_phase``", "Molar entropy indexed by phase", "J/mol.K"
   "``entr_mol``", "Molar entropy of mixture", "J/mol.K"
   "``fug_phase_comp``", "Fugacity indexed by phase and component", "Pa"
   "``fug_coeff_phase_comp``", "Fugacity coefficient indexed by phase and component", "None"
   "``gibbs_mol_phase``", "Molar Gibbs energy indexed by phase", "J/mol"
   "``mw``", "Molecular weight of mixture", "kg/mol"
   "``mw_phase``", "Molecular weight by phase", "kg/mol"
   "``temperature_bubble``", "Bubble point temperature", "K"
   "``temperature_dew``", "Dew point temperature", "K"
   "``pressure_bubble``", "Bubble point pressure", "Pa"
   "``pressure_dew``", "Dew point pressure", "Pa"
   "``_teq``", "Temperature at which the VLE is calculated", "K"

List of Parameters
------------------
.. csv-table::
   :header: "Parameter Name", "Description", "Units"

   "``cubic_type``", "Type of cubic equation of state to use, from CubicEoS Enum", "None"
   "``pressure_ref``", "Reference pressure", "Pa"
   "``temperature_ref``", "Reference temperature", "K"
   "``omega``", "Pitzer acentricity factor", "None"
   "``kappa``", "Binary interaction paraters for EoS", "None"
   "``mw_comp``", "Component molecular weights", "kg/mol"
   "``cp_ig``", "Parameters for calculating component heat capacities", "varies"
   "``dh_form``", "Component standard heats of formation (used for enthalpy at reference state)", "J/mol"
   "``ds_form``", "Component standard entropies of formation (used for entropy at reference state)", "J/mol.K"
   "``antoine``", "Component Antoine coefficients (used to initialize bubble and dew point calculations)", "bar, K"

Config Block Documentation
--------------------------
.. module:: idaes.property_models.cubic_eos.cubic_prop_pack

.. autoclass:: CubicParameterData
   :members:

.. autoclass:: CubicStateBlock
   :members:
.. autoclass:: CubicStateBlockData
   :members:

