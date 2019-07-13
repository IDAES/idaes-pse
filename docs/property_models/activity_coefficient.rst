Vapor-Liquid Equilibrium Property Models (Ideal Gas - Non-ideal Liquids)
=========================================================================

This property package supports phase equilibrium calucations with a smooth phase transition formulation that makes it amenable for equation oriented optimization. The gas phase is assumed to be ideal and for the liquid phase,
the package supports an ideal liquid or a non-ideal liquid using an activity
coefficient model. To compute the activity coefficient, the package currently supports the Non Random Two Liquid Model (NRTL) or the
Wilson model. Therefore, this property package supports the following combinations for gas-liquid mixtures for VLE calculations:

1. Ideal - Ideal
2. Ideal - NRTL
3. Ideal - Wilson

**Flow basis**: Molar

**Units**: SI units

**State Variables**: 

.. raw:: html
	<head>
   	</head>
	
   	<body>
      	<ul>
 	<li>Total molar flow rate (mol/s) - <code> <font color="red"> flow_mol </font> </code>
 	<li>Temperature (K) - <code> <font color="red"> temperature </font> </code>
	<li>Presure (Pa) - <code> <font color="red"> pressure</font></code>
	<li>Mole fraction of the mixture - <code> <font color="red"> mole_frac</font></code> 
	</ul>
	</body>

Support for other combinations of state variables will be made available in the future.


Inputs
------

When instantiating the parameter block that uses this particular state block, 2 arguments can be passed:

.. raw:: html
	<head>

   	</head>
   	<body>
      	<ul>
 	<li><code> <font color="red"> valid_phase </font> </code> - <font color="green">"Liq"</font>  or <font color="green">"Vap"</font> or <font color="green">("Liq", "Vap") </font> or <font color="green">("Vap", "Liq") </font> </li></li>
 	<li><code> <font color="red"> activity_coeff_model </font> </code> - <font color="green">"Ideal"</font> or <font color="green">"NRTL"</font> or <font color="green">"Wilson" </font></li>
	</ul>
	</body>

The ``valid_phase`` argument denotes the phases a user expects for a given set of inlet conditions. For example, if the user knows with certainty that the valid phase will be only Liquid or Vapor, then it is best not include the complex flash equilibrium constraints in the model. If the user does not specify any option, then the package defaults to a 2 phase assumption meaning that the phase equilibrium will be computed.

The ``activity_coeff_model`` denotes the liquid phase assumption to be used. If the user does not specify any option, then the package defaults to asuming an ideal liquid assumption.


Degrees of Freedom
------------------
The number of degrees of freedom that need to be fixed to yield a square problem (i.e. degrees of freedom = 0) depends on the options selected. The following table provides a summary of the variables to be fixed and also the corresponding variable names in the model. 

.. csv-table::
   :header: "Property Model Type", "State variables", "Additional Variables", "Total number of variables"
   :widths: 25, 15, 10, 30

   "Ideal-Ideal", "``flow_mol``, ``temperature``, ``pressure``, ``mole_frac``", "None", "3 + :math:`N_{c}`"
   "Ideal-NRTL", "``flow_mol``, ``temperature``, ``pressure``, ``mole_frac``", "``alpha``, ``tau``", "3 + :math:`N_{c}` + :math:`2N_{c}^{2}`"
   "Ideal-Wilson", "``flow_mol``, ``temperature``, ``pressure``, ``mole_frac``", "``vol_mol_comp``, ``tau``", "3 + :math:`N_{c}` + :math:`2N_{c}^{2}`"

Please refer to reference 3 for recommended values for ``tau``.


VLE Model with Smooth Phase Transition
--------------------------------------

The flash equations consists of the total mass balance, component balances and the equilibrium condition. 

.. math:: F^{in} = F^{liq} + F^{vap}
.. math:: z_{i}^{in}F^{in} = x_{i}^{liq}F^{liq} + y_{i}^{vap}F^{vap}
.. math:: f_{i}^{vap} = f_{i}^{liq}

The fugacity of the vapor and liquid phase is defined as follows:

.. math:: f_{i}^{vap} = y_{i}\phi_{i}P
.. math:: f_{i}^{liq} = x_{i}p^{sat}_{i}\nu_{i}

The equilibrium constraint is written as a generic constraint such that it can be extended easily for non-ideal gases and liquids. As this property package only supports an ideal gas, the fugacity coefficient (:math:`\phi_{i}`) for the vapor phase is 1 and hence the expression reduces to :math:`y_{i}P`. For the liquid phase, if the ideal option is selected then the activity coefficient (:math:`\nu_{i}`) is 1. If an activity coefficient model is selected then corresponding constraints are added to compute the activity coefficient. 

Typically, the flash caluclations are computed at a given temperature, :math:`T`. However, the flash calculations become trivial if the given conditions do not fall in the two phase region. For simulation only studies, the user may know a priori the condition of the stream but when the same set of equations are used for optimization, there is a high probablity that the specifications can transcend the phase envelope and hence the flash equations included may be trivial in the single phase region (i.e. liquid or vapor only). To circumvent this problem, property packages in IDAES that support VLE will compute the flash calculations at an "equivalent" temperature :math:`T_{eq}`. The equivalent temperature is computed as follows:
.. math:: T_{1} = max(T_{bubble}, T) 
.. math:: T_{eq} = min(T_{1}, T_{dew})

where :math:`T_[eq}` is the equilibrium temperature at which flash calculations are computed, :math:`T` is the stream temperature, :math:`T_{1}` is the intermediate temperature variable, :math:`T_{bubble}` is the bubble point temperature of mixture, and :math:`T_{dew}` is the dew point temperature of the mixture. Note that, in the above equations, approximations are used for the max and min functions as follows:
.. math:: T_{1} = 0.5(T + T_{bubble} + \sqrt((T-T_{bubble})^2 + \epsilon_{1}^2))
.. math:: T_{eq} = 0.5(T_{1} + T_{dew} - \sqrt((T-T_{dew})^2 + \epsilon_{2}^2))

where :math:`\epsilon_1` and :math:`\epsilon_2` are smoothing parameters(mutable). The default values are 0.01 and 0.0005 respectively. It is recommended that :math:`\epsilon_1` > :math:`\epsilon_2`. Please refer to reference 4 for more details. 

Additional constraints are included in the model to compute the thermodynamic properties such as component saturation pressure, enthalpy, specific heat capacity. Please note that, these constraints are added only if the variable is called for when building the model. This eliminates adding unnecessary constraints to compute properties that are not needed in the model. 

The saturation or vapor pressure (``pressure_sat``) for component :math:`i` is computed using the following correlation[1]:

.. math:: \log_e\frac{P^{sat}}{P_c} = \frac{Ax+Bx^{1.5}+Cx^3+Dx^6}{1-x}
.. math:: x=1-\frac{T}{T_c}

where :math:`P_c` is the critical pressure, :math:`T_c` is the critical temperature of the component and :math:`T` is the tempertaure at which the saturation pressure is computed. Please note that when using this expression, :math:`T<T_{c}` is required and when violated it results in a negative number raised to the power of a fraction. When using the smooth flash formulation, the saturation pressure is computed at the "equilibrium" temperature. This ensures that the VLE constraints are satusfied even when a single phase is detected. 

The specific enthalpy (``enthalpy_comp_liq``) for component :math:`i` is computed using the following expression for the liquid phase:

.. math:: h_{i}^{liq} =  \int_{298.15}^{T}(A+BT+CT^2+DT^3+ET^4)dT

The specific enthalpy (``enthalpy_comp_vap``) for component :math:`i` is computed using the following expression for the vapor phase:

.. math:: h_{i}^{vap} = \Delta H_{vap,298.15} + \int_{298.15}^{T}(A+BT+CT^2+DT^3+ET^4)dT

The mixture specific enthapies (``enthalpy_liq`` & ``enthalpy_vap``) are computed using the following expressions for the liquid and vapor phase respectively:

.. math:: H^{liq} =  \sum_i{h_{i}^{liq}x_{i}}
.. math:: H^{vap} =  \sum_i{h_{i}^{vap}y_{i}}

Activity Coefficient Model - NRTL
----------------------------------
The activity coefficient for component :math:`i` is computed using the following equations when using the Non-Random Two Liquid model [3]:

.. math:: log_e{\gamma_i} = \frac{\sum_j{x_j\tau_jG_{ji}}}{\sum_kx_kG_{ki}} + \sum_j\frac{x_jG_{ij}}{\sum_kx_kG_{kj}}\lbrack\tau_{ij} - \frac{\sum_mx_m\tau_{mj}G_{mj}}{\sum_kx_kG_{kj}}\rbrack
.. math:: G_{ij}=\exp({-\alpha_{ij}\tau_{ij}})

where :math:`\alpha_{ij}` is the non-randomness parameter and :math:`\tau_{ij}` is the binary interaction parameter for the NRTL model. Note that in the IDAES implementation, these are declared as variables that allows for more flexibility and the ability to use these in a parameter estimation problem. These NRTL model specific variables need to be either fixed for a given component set or need to be estimated from VLE data.  

The bubble point is computed by enforcing the following condition:

.. math:: \sum_i{\lbrack z_{i}p^{sat}_{i}(T_{bubble})\nu_{i}}\rbrack-P=0

Activity Coefficient Model - Wilson
-----------------------------------
The activity coefficient for component :math:`i` is computed using the following equations when using the Wilson model [3]:

.. math:: log_e{\gamma_i} = \frac{\sum_j{x_j\tau_jG_{ji}}}{\sum_kx_kG_{ki}} + \sum_j\frac{x_jG_{ij}}{\sum_kx_kG_{kj}}\lbrack\tau_{ij} - \frac{\sum_mx_m\tau_{mj}G_{mj}}{\sum_kx_kG_{kj}}\rbrack
.. math:: G_{ij}=\exp({-\alpha_{ij}\tau_{ij}})


where :math:`\alpha_{ij}` and :math:`\tau_{ij}` are NRTL model specific variables that either need to be fixed for a given component set or need to be estimated from VLE data.  

The bubble point is computed by enforcing the following condition:

.. math:: \sum_i{\lbrack z_{i}p^{sat}_{i}(T_{bubble})\nu_{i}}\rbrack-P=0

List of Variables
-----------------
.. csv-table::
   :header: "Variable Name", "Description", "Units"

   "``flow_mol``", "Total molar flow rate", "mol/s"
   "``mole_frac``", "Mixture mole fraction indexed by component", "None"
   "``temperature``", "Temperature", "K"
   "``pressure``", "Pressure", "Pa"
   "``flow_mol_phase``", "Molar flow rate indexed by phase", "mol/s"
   "``mole_frac_phase``", "Mole fraction indexed by phase and component", "None"
   "``pressure_sat``", "Saturation or vapor pressure indexed by component", "Pa"
   "``density_mol_phase``", "Molar density indexed by phase", "mol/m3"
   "``enthalpy_comp_liq``", "Liquid molar enthalpy indexed by component", "J/mol"
   "``enthalpy_comp_vap``", "Vapor molar enthalpy indexed by component", "J/mol"
   "``enthalpy_liq``", "Liquid phase enthalpy", "J/mol"
   "``enthalpy_vap``", "Vapor phase enthalpy", "J/mol"
   "``temperature_bubble``", "Bubble point temperature", "K"
   "``temperature_dew``", "Dew point temperature", "K"
   "``_temperature_equilibrium``", "Temperature at which the VLE is calculated", "K"

.. csv-table:: NRTL model specific variables
   :header: "Variable Name", "Description", "Units"

   "``alpha``", "Non-randomness parameter indexed by component and component", "None"
   "``tau``", "Binary interaction parameter indexed by component and component", "None"
   "``activity_coeff_comp``", "Activity coefficient indexed by component", "None"

.. csv-table:: Wilson model specific variables
   :header: "Variable Name", "Description", "Units"

   "``vol_mol_comp``", "Molar volume of component indexed by component", "None"
   "``tau``", "Binary interaction parameter indexed by component and component", "None"
   "``activity_coeff_comp``", "Activity coefficient indexed by component", "None"

Initialization
--------------

Config Block Documentation
--------------------------
.. module:: idaes.property_models.activity_coeff_models.activity_coeff_prop_pack

.. autoclass:: ActivityCoeffParameterBlock
   :members:

.. autoclass:: ActivityCoeffStateBlock
   :members:
.. autoclass:: ActivityCoeffStateBlockData
   :members:

References
------------

1. "The properties of gases and liquids by Robert C. Reid"
2. "Perry's Chemical Engineers Handbook by Robert H. Perry".
3. H. Renon and J.M. Prausnitz, "Local compositions in thermodynamic excess
   functions for liquid mixtures.", AIChE Journal Vol. 14, No.1, 1968.
4. AP Burgard, JP Eason, JC Eslick, JH Ghouse, A Lee, LT Biegler, DC Miller. "A Smooth, Square Flash Formulation for Equation Oriented Flowsheet Optimization",
   Computer Aided Chemical Engineering 44, 871-876, 2018


