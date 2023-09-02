Cubic Equations of State
========================

This property package implements a general form of a cubic equation of state which can be used for most cubic-type equations of state. This package supports phase equilibrium calculations with a smooth phase transition formulation that makes it amenable for equation oriented optimization. The following equations of state are currently supported:

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
	<li>Pressure (Pa) - <code> <font color="red"> pressure</font></code>
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
In general, the general cubic equation of state has a number of degrees of freedom equal to 2 + the number of components in the system (total flow rate, temperature, pressure and N-1 mole fractions). In some cases (primarily inlets to units), this is increased by 1 due to the removal of a constraint on the sum of mole fractions. 

General Cubic Equation of State
-------------------------------
All equations come from "The Properties of Gases and Liquids, 4th Edition" by Reid, Prausnitz and Poling. The general cubic equation of state is represented by the following equations:

.. math:: 0 = Z^3 - (1+B-uB)Z^2 + (A-uB-(u-w)B^2)Z - AB-wB^2-wB^3
.. math:: A = \frac{a_mP}{R^2T^2}
.. math:: B = \frac{b_mP}{RT}

where :math:`Z` is the compressibility factor of the mixture, :math:`a_m` and :math:`b_m` are properties of the mixture and :math:`u` and :math:`w` are parameters which depend on the specific equation of state being used as show in the table below.

.. csv-table::
   :header: "Equation", ":math:`u`", ":math:`w`", ":math:`\Omega_A`", ":math:`\Omega_B`", ":math:`\kappa_j`"

   "Peng-Robinson", "2", "-1", "0.45724", "0.07780", ":math:`(1+(1-T_r^2)(0.37464+1.54226\omega_j-0.26992\omega_j^2))^2`"
   "Soave-Redlich-Kwong", "1", "0", "0.42748", "0.08664", ":math:`(1+(1-T_r^2)(0.48+1.574\omega_j-0.176\omega_j^2))^2`"

The properties :math:`a_m` and :math:`b_m` are calculated from component specific properties :math:`a_j` and :math:`b_j` as shown below:

.. math:: a_j = \frac{\Omega_AR^2T_{c,j}^2}{P_{c, j}}\kappa_j
.. math:: b_j = \frac{\Omega_BRT_{c,j}}{P_{c,j}}
.. math:: a_m = \sum_i{\sum_j{y_iy_j(a_ia_j)^{1/2}(1-k_{ij})}}
.. math:: b_m = \sum_i{y_ib_i}

where :math:`P_{c,j}` and :math:`T_{c,j}` are the component critical pressures and temperatures, :math:`y_j` is the mole fraction of component :math:`j`, :math:`k_{ij}` are a set of binary interaction parameters which are specific to the equation of state and :math:`\Omega_A`, :math:`\Omega_B` and :math:`\kappa_j` are taken from the table above. :math:`\omega_j` is the Pitzer acentric factor of each component.

The cubic equation of state is solved for each phase via a call to an external function which automatically identifies the correct root of the cubic and returns the value of :math:`Z` as a function of :math:`A` and :math:`B` along with the first and second partial derivatives.

VLE Model with Smooth Phase Transition
--------------------------------------

The flash equations consists of the following equations:

.. math:: F^{in} = F^{liq} + F^{vap}
.. math:: z_{i}^{in}F^{in} = x_{i}^{liq}F^{liq} + y_{i}^{vap}F^{vap}

At the equilibrium condition, the fugacity of the vapor and liquid phase are defined as follows:

.. math:: \ln{f_{i}^{vap}} = \ln{f_{i}^{liq}}
.. math:: f_{i}^{phase} = y_{i}^{phase}\phi_{i}^{phase}P
.. math:: \ln{\phi_{i}} = \frac{b_i}{b_m}(Z-1) - \ln{(Z-B)} + \frac{A}{B\sqrt{u^2-4w}}\left(\frac{b_i}{b_m}-\delta_i\right)\ln{\left(\frac{2Z+B(u+\sqrt{u^2-4w})}{2Z+B(u-\sqrt{u^2-4w})}\right)}
.. math:: \delta_i = \frac{2a_i^{1/2}}{a_m}\sum_j{x_ja_j^{1/2}(1-k_{ij})}

The cubic equation of state is solved to find :math:`Z` for each phase subject to the composition of that phase. Typically, the flash calculations are computed at a given temperature, :math:`T`. However, the flash calculations become trivial if the given conditions do not fall in the two phase region. For simulation only studies, the user may know a priori the condition of the stream but when the same set of equations are used for optimization, there is a high probability that the specifications can transcend the phase envelope and hence the flash equations included may be trivial in the single phase region (i.e. liquid or vapor only). To circumvent this problem, property packages in IDAES that support VLE will compute the flash calculations at an "equilibrium" temperature :math:`T_{eq}`. The equilibrium temperature is computed as follows:

.. math:: T_{1} = max(T_{bubble}, T) 
.. math:: T_{eq} = min(T_{1}, T_{dew})

where :math:`T_{eq}` is the equilibrium temperature at which flash calculations are computed, :math:`T` is the stream temperature, :math:`T_{1}` is the intermediate temperature variable, :math:`T_{bubble}` is the bubble point temperature of mixture, and :math:`T_{dew}` is the dew point temperature of the mixture. Note that, in the above equations, approximations are used for the max and min functions as follows:

.. math:: T_{1} = 0.5{[T + T_{bubble} + \sqrt{(T-T_{bubble})^2 + \epsilon_{1}^2}]}
.. math:: T_{eq} = 0.5{[T_{1} + T_{dew} - \sqrt{(T-T_{dew})^2 + \epsilon_{2}^2}]}

where :math:`\epsilon_1` and :math:`\epsilon_2` are smoothing parameters (mutable). The default values are 0.01 and 0.0005 respectively. It is recommended that :math:`\epsilon_1` > :math:`\epsilon_2`. Please refer to reference 4 for more details. Therefore, it can be seen that if the stream temperature is less than that of the bubble point temperature, the VLE calculations will be computed at the bubble point. Similarly, if the stream temperature is greater than the dew point temperature, then the VLE calculations are computed at the dew point temperature. For all other conditions, the equilibrium calculations will be computed at the actual temperature.

Other Constraints
-----------------

Additional constraints are included in the model to compute the thermodynamic properties based on the cubic equation of state, such as enthalpies and entropies. Please note that, these constraints are added only if the variable is called for when building the model. This eliminates adding unnecessary constraints to compute properties that are not needed in the model.

All thermophysical properties are calculated using an ideal and residual term, such that:

.. math:: p = p^0 + p^r

The residual term is derived from the partial derivatives of the cubic equation of state, whilst the ideal term is determined using empirical correlations.

Enthalpy
^^^^^^^^

The ideal enthalpy term is given by:

.. math:: h_{i}^{0} =  \int_{298.15}^{T}(A+BT+CT^2+DT^3)dT + \Delta h_{form}^{298.15K}

The residual enthalpy term is given by:

.. math:: h_{i}^{r}b_m\sqrt{u^2-4w} = \left(T\frac{da}{dT}-a_m\right)\ln{\left(\frac{2Z+B(u+\sqrt{u^2-4w})}{2Z+B(u-\sqrt{u^2-4w})}\right)} +RT(Z-1)b_m\sqrt{u^2-4w}

.. math:: \frac{da}{dT}\sqrt{T} = -\frac{R}{2}\sqrt{\Omega_A}\sum_i{\sum_j{y_iy_j(1-k_{ij})\left(f_{w,j}\sqrt{a_i\frac{T_{c,j}}{P_{c,j}}}+f_{w,i}\sqrt{a_j\frac{T_{c,i}}{P_{c,i}}}\right)}}

Entropy
^^^^^^^^

The ideal entropy term is given by:

.. math:: s_{i}^{0} =  \int_{298.15}^{T}\frac{(A+BT+CT^2+DT^3)}{T}dT + \Delta s_{form}^{298.15K}

The residual entropy term is given by:

.. math:: s_{i}^{r}b_m\sqrt{u^2-4w} = R\ln{\frac{Z-B}{Z}}b_m\sqrt{u^2-4w} + R\ln{\frac{ZP^{ref}}{P}}b_m\sqrt{u^2-4w} + \frac{da}{dT}\ln{\left(\frac{2Z+B(u+\sqrt{u^2-4w})}{2Z+B(u-\sqrt{u^2-4w})}\right)}

Fugacity
^^^^^^^^

Fugacity is calculated from the system pressure, mole fractions and fugacity coefficients as follows:

.. math :: f_{i, p} = x_{i, p} \phi_{i, p} P

Fugacity Coefficient
^^^^^^^^^^^^^^^^^^^^

The fugacity coefficient is calculated from the departure function of the cubic equation of state as shown below:

.. math:: \ln{\phi_{i}} = \frac{b_i}{b_m}(Z-1) - \ln{(Z-B)} + \frac{A}{B\sqrt{u^2-4w}}\left(\frac{b_i}{b_m}-\delta_i\right)\ln{\left(\frac{2Z+B(u+\sqrt{u^2-4w})}{2Z+B(u-\sqrt{u^2-4w})}\right)}

.. math:: \delta_i = \frac{2a_i^{1/2}}{a_m} \sum_j{x_j a_j^{1/2}(1-k_{ij})}

Gibbs Energy
^^^^^^^^^^^^

The Gibbs energy of the system is calculated using the definition of Gibbs energy:

.. math:: g_i = h_i - T \Delta s_i

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
   "``cp_mol_phase``", "Isobaric molar heat capacity by phase", "J/mol/K"
   "``cv_mol_phase``", "Isochoric molar heat capacity by phase", "J/mol/K"
   "``dens_mol_phase``", "Molar density indexed by phase", "mol/m3"
   "``dens_mass_phase``", "Mass density indexed by phase", "kg/m3"
   "``enth_mol_phase``", "Molar enthalpy indexed by phase ", "J/mol"
   "``enth_mol``", "Molar enthalpy of mixture", "J/mol"
   "``entr_mol_phase``", "Molar entropy indexed by phase", "J/mol.K"
   "``entr_mol``", "Molar entropy of mixture", "J/mol.K"
   "``fug_phase_comp``", "Fugacity indexed by phase and component", "Pa"
   "``fug_coeff_phase_comp``", "Fugacity coefficient indexed by phase and component", "None"
   "``gibbs_mol_phase``", "Molar Gibbs energy indexed by phase", "J/mol"
   "``heat_capacity_ratio_phase``", "Heat capcity ratio by phase", "-"
   "``isothermal_speed_sound_phase``", "Isothermal speed of sound by phase", "m/s"
   "``isentropic_speed_sound_phase``", "Isentropic speed of sound by phase", "m/s"
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
   "``kappa``", "Binary interaction parameters for EoS (note that parameters are specific for a given EoS", "None"
   "``mw_comp``", "Component molecular weights", "kg/mol"
   "``cp_ig``", "Parameters for calculating component heat capacities", "varies"
   "``dh_form``", "Component standard heats of formation (used for enthalpy at reference state)", "J/mol"
   "``ds_form``", "Component standard entropies of formation (used for entropy at reference state)", "J/mol.K"
   "``antoine``", "Component Antoine coefficients (used to initialize bubble and dew point calculations)", "bar, K"

Initialization
--------------

.. module:: idaes.models.properties.cubic_eos.cubic_prop_pack

.. autoclass:: CubicEoSInitializer
   :members: initialization_routine

Config Block Documentation
--------------------------

.. autoclass:: CubicParameterData
   :members:

.. autoclass:: CubicStateBlock
   :members:
.. autoclass:: CubicStateBlockData
   :members:

