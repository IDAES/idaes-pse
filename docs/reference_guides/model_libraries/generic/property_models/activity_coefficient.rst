Vapor-Liquid Equilibrium Property Models (Ideal Gas - Non-ideal Liquids)
=========================================================================

This property package supports phase equilibrium calucations with a smooth phase transition formulation that makes it amenable for equation oriented optimization. The gas phase is assumed to be ideal and for the liquid phase,
the package supports an ideal liquid or a non-ideal liquid using an activity
coefficient model. To compute the activity coefficient, the package currently supports the Non Random Two Liquid Model (NRTL) or the
Wilson model. Therefore, this property package supports the following combinations for gas-liquid mixtures for VLE calculations:

1. Ideal (vapor) - Ideal (liquid)
2. Ideal (vapor) - NRTL (liquid)
3. Ideal (vapor) - Wilson (liquid)

**Flow basis**: Molar

**Units**: SI units

**State Variables**: 

The state block supports the following two sets of state variables:

Option 1 - "FTPz":

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

Option 2 - "FcTP":

.. raw:: html
	<head>
   	</head>
	
   	<body>
      	<ul>
 	<li>Component molar flow rate (mol/s) - <code> <font color="red"> flow_mol_comp </font> </code>
 	<li>Temperature (K) - <code> <font color="red"> temperature </font> </code>
	<li>Presure (Pa) - <code> <font color="red"> pressure</font></code>
	</ul>
	</body>


The user can specify the choice of state variables while instantiating the parameter block. See the Inputs section for more details. 

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
	<li><code> <font color="red"> state_vars </font> </code> - <font color="green">"FTPz"</font> or <font color="green">"FcTP"</font> </li>
	</ul>
	</body>

The ``valid_phase`` argument denotes the valid phases for a given set of inlet conditions. For example, if the user knows a priori that the it will only be a single phase (for example liquid only), then it is best not to include the complex flash equilibrium constraints in the model. If the user does not specify any option, then the package defaults to a 2 phase assumption meaning that the constraints to compute the phase equilibrium will be computed.

The ``activity_coeff_model`` denotes the liquid phase assumption to be used. If the user does not specify any option, then the package defaults to assuming an ideal liquid assumption.

The ``state_vars`` denotes the preferred set of state variables to be used. If the user does not specify any option, then the package defaults to using the total flow, mixture mole fraction, temperature and pressure as the state variables.



Degrees of Freedom
------------------
The number of degrees of freedom that need to be fixed to yield a square problem (i.e. degrees of freedom = 0) depends on the options selected. The following table provides a summary of the variables to be fixed and also the corresponding variable names in the model. 

.. csv-table::
   :header: "Property Model Type", "State variables", "Additional Variables", "Total number of variables"
   :widths: 25, 15, 10, 30

   "Ideal (vapor) - Ideal (liquid)", "``flow_mol``, ``temperature``, ``pressure``, ``mole_frac_comp``", "None", "3 + :math:`N_{c}`"
   "Ideal (vapor) - NRTL (liquid)", "``flow_mol``, ``temperature``, ``pressure``, ``mole_frac_comp``", "``alpha``, ``tau``", "3 + :math:`N_{c}` + :math:`2N_{c}^{2}`"
   "Ideal (vapor) - Wilson (liquid)", "``flow_mol``, ``temperature``, ``pressure``, ``mole_frac_comp``", "``vol_mol_comp``, ``tau``", "3 + :math:`N_{c}` + :math:`2N_{c}^{2}`"

Please refer to reference 3 for recommended values for ``tau``.


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

The equilibrium constraint is written as a generic constraint such that it can be extended easily for non-ideal gases and liquids. As this property package only supports an ideal gas, the fugacity coefficient (:math:`\phi_{i}`) for the vapor phase is 1 and hence the expression reduces to :math:`y_{i}P`. For the liquid phase, if the ideal option is selected then the activity coefficient (:math:`\nu_{i}`) is 1. If an activity coefficient model is selected then corresponding constraints are added to compute the activity coefficient. 

Typically, the flash calculations are computed at a given temperature, :math:`T`. However, the flash calculations become trivial if the given conditions do not fall in the two phase region. For simulation only studies, the user may know a priori the condition of the stream but when the same set of equations are used for optimization, there is a high probablity that the specifications can transcend the phase envelope and hence the flash equations included may be trivial in the single phase region (i.e. liquid or vapor only). To circumvent this problem, property packages in IDAES that support VLE will compute the flash calculations at an "equilibrium" temperature :math:`T_{eq}`. The equilibrium temperature is computed as follows:

.. math:: T_{1} = max(T_{bubble}, T) 
.. math:: T_{eq} = min(T_{1}, T_{dew})

where :math:`T_{eq}` is the equilibrium temperature at which flash calculations are computed, :math:`T` is the stream temperature, :math:`T_{1}` is the intermediate temperature variable, :math:`T_{bubble}` is the bubble point temperature of mixture, and :math:`T_{dew}` is the dew point temperature of the mixture. Note that, in the above equations, approximations are used for the max and min functions as follows:

.. math:: T_{1} = 0.5{[T + T_{bubble} + \sqrt{(T-T_{bubble})^2 + \epsilon_{1}^2}]}
.. math:: T_{eq} = 0.5{[T_{1} + T_{dew} - \sqrt{(T-T_{dew})^2 + \epsilon_{2}^2}]}

where :math:`\epsilon_1` and :math:`\epsilon_2` are smoothing parameters(mutable). The default values are 0.01 and 0.0005 respectively. It is recommended that :math:`\epsilon_1` > :math:`\epsilon_2`. Please refer to reference 4 for more details. Therefore, it can be seen that if the stream temperature is less than that of the bubble point temperature, the VLE calucalations will be computed at the bubble point. Similarly, if the stream temperature is greater than the dew point temperature, then the VLE calculations are computed at the dew point temperature. For all other conditions, the equilibrium calculations will be computed at the actual temperature. 

Additional constraints are included in the model to compute the thermodynamic properties such as component saturation pressure, enthalpy, specific heat capacity. Please note that, these constraints are added only if the variable is called for when building the model. This eliminates adding unnecessary constraints to compute properties that are not needed in the model. 

The saturation or vapor pressure (``pressure_sat``) for component :math:`i` is computed using the following correlation[1]:

.. math:: \log\frac{P^{sat}}{P_c} = \frac{Ax+Bx^{1.5}+Cx^3+Dx^6}{1-x}
.. math:: x=1-\frac{T_{eq}}{T_c}

where :math:`P_c` is the critical pressure, :math:`T_c` is the critical temperature of the component and :math:`T_{eq}` is the equilibrium temperature at which the saturation pressure is computed. Please note that when using this expression, :math:`T_{eq}<T_{c}` is required and when violated it results in a negative number raised to the power of a fraction. 

The specific enthalpy (``enthalpy_comp_liq``) for component :math:`i` is computed using the following expression for the liquid phase:

.. math:: h_{i}^{liq} =  \Delta h_{form,Liq,i} + \int_{298.15}^{T}(A+BT+CT^2+DT^3+ET^4)dT

The specific enthalpy (``enthalpy_comp_vap``) for component :math:`i` is computed using the following expression for the vapor phase:

.. math:: h_{i}^{vap} = \Delta h_{form,Vap,i} + \int_{298.15}^{T}(A+BT+CT^2+DT^3+ET^4)dT

The mixture specific enthapies (``enthalpy_liq`` & ``enthalpy_vap``) are computed using the following expressions for the liquid and vapor phase respectively:

.. math:: H^{liq} =  \sum_i{h_{i}^{liq}x_{i}}
.. math:: H^{vap} =  \sum_i{h_{i}^{vap}y_{i}}

Similarly, specific entropies are calculated as follows. The specific entropy (``entropy_comp_liq``) for component :math:`i` is computed using the following expression for the liquid phase:

.. math:: s_{i}^{liq} =  \Delta s_{form,Liq,i} + \int_{298.15}^{T}(A/T+B+CT+DT^2+ET^3)dT

The specific entropy (``entropy_comp_vap``) for component :math:`i` is computed using the following expression for the vapor phase:

.. math:: s_{i}^{vap} = \Delta s_{form,Vap,i} + \int_{298.15}^{T}(A/T+B+CT+DT^2+ET^3)dT

Please refer to references 1 and 2 to get parameters for different components. 

Activity Coefficient Model - NRTL
----------------------------------
The activity coefficient for component :math:`i` is computed using the following equations when using the Non-Random Two Liquid model [3]:

.. math:: \log{\gamma_i} = \frac{\sum_j{x_j\tau_jG_{ji}}}{\sum_kx_kG_{ki}} + \sum_j\frac{x_jG_{ij}}{\sum_kx_kG_{kj}}\lbrack\tau_{ij} - \frac{\sum_mx_m\tau_{mj}G_{mj}}{\sum_kx_kG_{kj}}\rbrack
.. math:: G_{ij}=\exp({-\alpha_{ij}\tau_{ij}})

where :math:`\alpha_{ij}` is the non-randomness parameter and :math:`\tau_{ij}` is the binary interaction parameter for the NRTL model. Note that in the IDAES implementation, these are declared as variables that allows for more flexibility and the ability to use these in a parameter estimation problem. These NRTL model specific variables need to be either fixed for a given component set or need to be estimated from VLE data.  

The bubble point is computed by enforcing the following condition:

.. math:: \sum_i{\lbrack z_{i}p^{sat}_{i}(T_{bubble})\nu_{i}}\rbrack-P=0

Activity Coefficient Model - Wilson
-----------------------------------
The activity coefficient for component :math:`i` is computed using the following equations when using the Wilson model [3]:

.. math:: \log{\gamma_i} = 1 - \log{\sum_jx_jG_{ji}} - \sum_j\frac{x_jG_{ij}}{\sum_kx_kG_{kj}}
.. math:: G_{ij}=(v_{i}/v_{j})\exp(-\tau_{ij})


where :math:`v_{i}` is the molar volume of component :math:`i` and :math:`\tau_{ij}` is the binary interaction parameter. These are Wilson model specific variables that either need to be fixed for a given component set or need to be estimated from VLE data.   

The bubble point is computed by enforcing the following condition:

.. math:: \sum_i{\lbrack z_{i}p^{sat}_{i}(T_{bubble})\nu_{i}}\rbrack-P=0

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
   "``density_mol_phase``", "Molar density indexed by phase", "mol/m3"
   "``ds_vap``", "Molar entropy of vaporization", "J/mol.K" 
   "``enthalpy_comp_liq``", "Liquid molar enthalpy indexed by component", "J/mol"
   "``enthalpy_comp_vap``", "Vapor molar enthalpy indexed by component", "J/mol"
   "``enthalpy_liq``", "Liquid phase enthalpy", "J/mol"
   "``enthalpy_vap``", "Vapor phase enthalpy", "J/mol"
   "``entropy_comp_liq``", "Liquid molar entropy indexed by component", "J/mol"
   "``entropy_comp_vap``", "Vapor molar entropy indexed by component", "J/mol"
   "``entrolpy_liq``", "Liquid phase entropy", "J/mol"
   "``entropy_vap``", "Vapor phase entropy", "J/mol"
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
.. module:: idaes.models.properties.activity_coeff_models.activity_coeff_prop_pack

.. autoclass:: ActivityCoeffInitializer
   :members: initialization_routine

Config Block Documentation
--------------------------

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


