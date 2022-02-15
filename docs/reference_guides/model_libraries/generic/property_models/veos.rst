Virial Equation of State (VEOS)
===============================

This property package implements the B-truncated (second virial coefficient) volume explicit virial equation of state for gases where the generalized second virial coefficient is adequately represented by the Abbort equations. The B-truncated VEOS has the numerical advantage of having a single root and adequately represents the vapor phase at low to moderate pressure or where the vapor phase reduced density is less than one-half. The VEOS representation for pure species or gas mixtures is given by the following equation:

.. math:: Z = 1 + \frac{BP}{RT}

where :math:`Z` is the compressibility factor and :math:`B` is the second virial coefficient which is substance dependent and a function of temperature.
B is defined for both pure species and mixtures as follows:

.. math:: \mathrm{Pure\:Species}: B = B_{ii}
.. math:: \mathrm{Mixture}     : B = \sum_i{\sum_j{y_i y_j B_{ij}}}

The general form of the second virial coeffcient and its temperature derivatives are defined as shown in the table below:

.. csv-table::
   :header: "Generic Second Virial Coefficient and its Temperature Derivatives "

   ":math:`B_{ij} = \frac{RT_{c,ij}}{P_{c,ij}} \left( B^0_{ij} + \omega_{ij} B^1_{ij}\right)`"
   ":math:`\frac{dB_{ij}}{dT} = \frac{R}{P_{c,ij}} \left(\frac{dB_{ij}^0}{dT_{r,ij}} + \omega_{ij} \frac{dB^1_{ij}}{dT_{r,ij}} \right)`"
   ":math:`B^0_{ij} = 0.083 - \frac{0.422}{{T_{r,ij}}^{1.6}}`"
   ":math:`B^1_{ij} = 0.139 - \frac{0.172}{{T_{r,ij}}^{4.2}}`"
   ":math:`\frac{dB^0_{ij}}{dT_{r,ij}} = \frac{0.6752}{{T_{r,ij}}^{2.6}}`"
   ":math:`\frac{dB^1_{ij}}{dT_{r,ij}} = \frac{0.7224}{{T_{r,ij}}^{5.2}}`"
   ":math:`\frac{d}{dT_{r,ij}}\left( \frac{dB^0_{ij}}{dT_{r,ij}} \right) = \frac{-1.7552}{{T_{r,ij}}^{3.6}}`"
   ":math:`\frac{d}{dT_{r,ij}}\left( \frac{dB^1_{ij}}{dT_{r,ij}} \right) = \frac{-3.75648}{{T_{r,ij}}^{6.2}}`"

Where the like molecular terms :math:`(i=j)` indicate pure component terms , i.e. :math:`B_{ii}` is the second virial coefficient of pure substance :math:`i` and the unlike molecular terms :math:`(i \ne j)` indicate interaction between different molecules for extension to mixtures. The cross coefficients are evaluated by using the combining rules shown in the table below:

.. csv-table::
   :header: "Combining Rules for Cross Coefficients"

   ":math:`\omega_{ij} = 0.5(\omega_i +\omega_j)`"
   ":math:`Z_{c,ij} = 0.5 (Z_{c,i} + Z_{c,j})`"
   ":math:`V_{c,ij} = 0.125 \biggl(V_{c,i}^{1/3} + V_{c,j}^{1/3} \biggr)^3`"
   ":math:`T_{c,ij} = (1-\kappa_{ij})(T_{c,i}T_{c,j})^{1/2}`"
   ":math:`P_{c,ij} = \frac{Z_{c,ij}RT_{c,ij}}{V_{c,ij}}`"

The VEOS property package also has the option of treating a gas mixture as a pseudocritical component by passing the option: `use_pseudocritical_rules` in the configuration dictionary. This option builds pseudocritical parameters that treats the mixture as a pure component as shown in the table below

.. csv-table::
   :header: "Pseudocrital Parameters"

   ":math:`\omega_{p} = \sum_i{ y_i \omega_i}`"
   ":math:`T_{c,p} = \sum_i{ y_i T_{c,i}}`"
   ":math:`P_{c,p} = \sum_i{ y_i P_{c,i}}`"
   ":math:`T_{r,p} = \frac{\mathrm{T}}{T_{c,p}}`"
   ":math:`P_{r,p} = \frac{\mathrm{P}}{P_{c,p}}`"


Other Constraints
-----------------

Additional constraints are included in the model to compute thermodynamic properties based on heat capacity data of an ideal gas state and the virial equation of state. The isobaric heat capacity :math:`\biggl(Cp_i^{ig}\biggr)`, molar enthalpy :math:`\biggl(H_i^{ig}\biggr)` and molar entropy  :math:`\biggl(S_i^{ig}\biggr)` of an ideal gas are described below:

.. math:: Cp_i^{ig}(T,P) = \int_{T_{ref}}^{T}(A+BT+CT^2+DT^3)dT
.. math:: H_i^{ig}(T,P) = H_{i,ref}\;\bigl(T_{ref},P_{ref}\bigr) +  \int_{T_{ref}}^{T}Cp_i^{ig}dT
.. math:: S_i^{ig}(T,P) = S_{i,ref}\;\bigl(T_{ref},P_{ref}\bigr) +  \int_{T_{ref}}^{T}\frac{Cp_i^{ig}}{T}dT -
                          R\;\ln\biggl(\frac{P}{P_{ref}}\biggr)


The thermodynamic properties are shown in the table below:

.. csv-table::
   :header: "Property", "Description", "Residual Term"

   ":math:`\mathbf{Pure\; Properties}`"
   "Log fugacity coefficient",":math:`\ln\phi_i = \frac{P_r}{T_r}\biggl(B^0+\omega_i B^1\biggr)`",
   "Fugacity",":math:`f_i=\phi_i P`",
   "Molar Volume",":math:`V_i = \frac{Z_i R T}{P}`",
   "Molar Enthalpy",":math:`H_i = H_i^{ig} + H_i^{R}`",":math:`H_i^{R}=RT_cP_c \biggl[B^0 -T_r\frac{dB^0}{dT_r} + \omega*\biggl(B^1-T_r\frac{dB^1}{dT_r}\biggr)\biggr]`"
   "Molar Entropy",":math:`S_i = S_i^{ig} + S_i^{R}`",":math:`S_i^{R}=RP_r \biggl(\frac{dB^0}{dT_r} + \omega*\frac{dB^1}{dT_r}\biggr)`"
   "Molar Gibbs Energy",":math:`G_i = H_i -TS_i`",
   "Molar Internal Energy",":math:`U_i = H_i -PV_i`",
   ":math:`\mathbf{Mixture\; Properties}`"
   "Molar Volume",":math:`V = \frac{Z R T}{P}`",
   "Molar Enthalpy",":math:`H = \sum_i{y_i H_i^{ig}} + H^{R}`",":math:`H^{R}=PT\;\biggl(\frac{B}{T} -\frac{dB}{dT}\biggr)`"
   "Molar Entropy",":math:`S = \sum_i{y_i S_i^{ig}} + S^{R}`",":math:`S^{R}=-P\;\frac{dB}{dT}`"
   "Molar Gibbs Energy",":math:`G = H -TS`",
   "Molar Internal Energy",":math:`U = H -PV`",
   "Isobaric molar heat capacity",":math:`Cp = \sum_i{Cp_i^{ig}}+ Cp^{R}`",":math:`Cp^{R}=-\left(\frac{PT}{T_{c}^2}\;\frac{d^2B}{dT_r^2}\right)`"
   "Isochoric molar heat capacity",":math:`Cv = Cp - \frac{1}{R}\biggl(R + P\frac{dB}{dT}\biggr)^2`"
   ":math:`\mathbf{Partial\;Molar\; Properties}`"
   "Log fugacity coefficient",":math:`\ln\widehat \phi_j = \frac{P}{RT}\;\biggl(2\sum_i{y_i B_{ij}} - B\biggr)`",
   "Partial molar Gibbs energy",":math:`\bar{G_j} = H_j^{ig} -TS_j^{ig}+RT\;\ln y_j + \bar G_j^R`",":math:`\bar G_j^R=RT\;\ln\widehat \phi_j`"
   "Partial molar entropy",":math:`\bar{S_j} = -\biggl(\frac{\bar{G_j}}{dT}\biggr)_{P,x}`"
   "Partial molar enthalpy",":math:`\bar{H_j} = \bar{G_j} + T\bar{S_j}`"
   "Partial molar volume",":math:`\bar{V_j} = \frac{RT}{P}\;\biggl(\ln\widehat \phi_j+1)`",
   "Partial molar internal energy",":math:`\bar{U_j} = \bar{H_j} - P\bar{V_j}`"

List of Variables
-----------------
.. csv-table::
   :header: "Variable Name", "Description", "Units"

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

List of Parameters
------------------
.. csv-table::
   :header: "Parameter Name", "Description", "SI Units"

   "``use_pseudocritical_rules``", "Option to treat mixture as a pseudocritical component", "None"
   "``pressure_ref``", "Reference pressure", "Pa"
   "``temperature_ref``", "Reference temperature", "K"
   "``temperature_crit``", "Component critical temperature", "K"
   "``pressure_crit``", "Component critical pressure", "Pa"
   "``volume_crit``", "Component critical volume", ":math:`\mathrm{m^3/mol}`"
   "``omega``", "Component acentricity factor", "None"
   "``kappa``", "Binary interaction parameters for EoS", "None"
   "``mw_comp``", "Component molecular weights", "kg/mol"
   "``dh_form``", "Component standard heats of formation (used for enthalpy at reference state)", "J/mol"
   "``ds_form``", "Component standard entropies of formation (used for entropy at reference state)", "J/mol.K"
   "``dg_form``", "Component standard Gibbs energy of formation (used to compute entropy at reference state)", "J/mol"

Config Block Documentation
--------------------------
.. module:: idaes.generic_models.properties.core.eos.virial

.. autoclass:: VirialParameterData
   :members:

.. autoclass:: VirialStateBlock
   :members:
.. autoclass:: VirialStateBlockData
   :members:

