.. _standards:

IDAES Modeling Standards
========================

.. contents:: Contents 
    :depth: 3


Model Formatting and General Standards
--------------------------------------
The section describes the recommended formatting used within the IDAES framework. Users are strongly encouraged to follow these standards in developing their models in order to improve readability of their code.

Headers and Meta-data
^^^^^^^^^^^^^^^^^^^^^
Model developers are encouraged to include some documentation in the header of their model files which provides a brief description of the purpose of the model and how it was developed. Some suggested information to include is:

* Model name,
* Model publication date,
* Model author
* Any necessary licensing and disclaimer information (see below).
* Any additional information the modeler feels should be included.

Coding Standard
^^^^^^^^^^^^^^^
All code developed as part of IDAES should conform to the PEP-8 standard.

Model Organization
^^^^^^^^^^^^^^^^^^
Whilst the overall IDAES modeling framework enforces a hierarchical structure on models, model developers are still encouraged to arrange their models in a logical fashion to aid other users in understanding the model. Model constraints should be grouped with similar constraints, and each grouping of constraints should be clearly commented. 

For property packages, it is recommended that all the equations necessary for calculating a given property be grouped together, clearly separated and identified by using comments.

Additionally, model developers are encouraged to consider breaking their model up into a number of smaller methods where this makes sense. This can facilitate modification of the code by allowing future users to inherit from the base model and selectively overloading sub-methods where desired.

Commenting
^^^^^^^^^^
To help other modelers and users understand the how a model works, model builders are strongly encouraged to comment their code. It is suggested that every constraint should be commented with a description of the purpose of the constraint, and if possible/necessary a reference to a source or more detailed explanation. Any deviations from standard units or formatting should be clearly identified here. Any initialization procedures, or other procedures required to get the model to converge should be clearly commented and explained where they appear in the code. Additionally, modelers are strongly encouraged to add additional comments explaining how their model works to aid others in understanding the model.

Units of Measurement and Reference States
-----------------------------------------
Due to the flexibility provided by the IDAES modeling framework, there is no standard set of units of measurement or standard reference state that should be used in models. This places the onus on the user to understand the units of measurement being used within their models and to ensure that they are consistent.

The IDAES developers have generally used SI units without prefixes (i.e. Pa, not kPa) within models developed by the institute, with a default thermodynamic reference state of 298.15 K and 101325 Pa. Supercritical fluids have been consider to be part of the liquid phase, as they will be handled via pumps rather than compressors.

Standard Variable Names
-----------------------
In order for different models to communicate information effectively, it is necessary to have a standard naming convention for any variable that may need to be shared between different models. Within the IDAES modeling framework, this occurs most frequently with information regarding the state and properties of the material within the system, which is calculated in specialized property blocks, and then used in others parts of the model. This section of the documentation discusses the standard naming conventions used within the IDAES modeling framework.

Standard Naming Format
^^^^^^^^^^^^^^^^^^^^^^
There are a wide range of different variables which may be of interest to modelers, and a number of different ways in which these quantities can be expressed. In order to facilitate communication between different parts of models, a naming convention has been established to standardize the naming of variables across models. Variable names within IDAES follow to format below::

    {property_name}_{basis}_{state}_{condition}

Here, property_name is the name of the quantity in question, and should be drawn from the list of standard variable names given later in this document. If a particular quantity is not included in the list of standard names, users are encouraged to contact the IDAES developers so that it can be included in a future release. This is followed by a number of qualifiers which further indicate the specific conditions under which the quantity is being calculated. These qualifiers are described below, and some examples are given at the end of this document.

Basis Qualifier
"""""""""""""""
Many properties of interest to modelers are most conveniently represented on an intensive basis, that is quantity per unit amount of material. There are a number of different bases that can be used when expressing intensive quantities, and a list of standard basis qualifiers are given below.

============ =============
Basis        Standard Name
============ =============
Mass Basis   mass
Molar Basis  mol
Volume Basis vol
============ =============

State Qualifier
"""""""""""""""
Many quantities can be calculated either for the whole or a part of a mixture. In these cases, a qualifier is added to the quantity to indicate which part of the mixture the quantity applies to. In these cases, quantities may also be indexed by a Pyomo Set.

================= ============= ===================================
Basis             Standard Name Comments
================= ============= ===================================
Component         comp          Indexed by component list
Phase             phase         Indexed by phase list
Phase & Component phase_comp    Indexed by phase and component list
Total Mixture                   No state qualifier
================= ============= ===================================

=================== =============
Phase               Standard Name
=================== =============
Supercritical Fluid liq
Ionic Species       ion
Liquid Phase        liq
Solid Phase         sol
Vapor Phase         vap
Multiple Phases     e.g. liq1
=================== =============

Condition Qualifier
"""""""""""""""""""
There are also cases where a modeler may want to calculate a quantity at some state other than the actual state of the system (e.g. at the critical point, or at equilibrium).

================== =============
Basis              Standard Name
================== =============
Critical Point     crit
Equilibrium State  equil
Ideal Gas          ideal
Reduced Properties red
Reference State    ref
================== =============

Constants
^^^^^^^^^

=============================== =============
Constant                        Standard Name
=============================== =============
Gas Constant	                gas_const
=============================== =============

Thermophysical and Transport Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Below is a list of all the thermophysical properties which currently have a standard name associated with them in the IDAES framework.

=============================== =====================
Variable                        Standard Name
=============================== =====================
Activity	                act
Activity Coefficient	        act_coeff
Bubble Temperature	        t_bub
Compressibility Factor	        compress_fact
Concentration	                conc
Density	                        dens
Dew Temperature	                temperature_dew
Diffusivity	                diffus
Diffusion Coefficient (binary)	diffus_binary
Enthalpy	                enth
Entropy	                        entr
Fugacity	                fug
Fugacity Coefficient	        fug_coeff
Gibbs Energy	                energy_gibbs
Heat Capacity (const. P)	cp
Heat Capacity (const. V)	cv
Heat Capacity Ratio	        heat_capacity_ratio
Helmholtz Energy	        energy_helmholtz
Henry's Constant	        henry
Mass Fraction	                mass_frac
Material Flow	                flow
Molecular Weight	        mw
Mole Fraction	                mole_frac
pH	                        pH
Pressure	                pressure
Speed of Sound                  speed_sound
Surface Tension	                surf_tens
Temperature	                temperature
Thermal Conductivity	        therm_cond
Vapor Pressure	                pressure_sat
Viscosity (dynamic)	        visc_d
Viscosity (kinematic)	        visc_k
Vapor Fraction                  vap_frac
Volume Fraction	                vol_frac
=============================== =====================

Reaction Properties
^^^^^^^^^^^^^^^^^^^
Below is a list of all the reaction properties which currently have a standard name associated with them in the IDAES framework.

======================= =================
Variable                Standard Name
======================= =================
Activation Energy       energy_activation
Arrhenius Coefficient   arrhenius
Heat of Reaction        dh_rxn
Entropy of Reaction     ds_rxn
Equilibrium Constant    k_eq
Reaction Rate           reaction_rate
Rate constant           k_rxn
Solubility Constant     k_sol
======================= =================

Solid Properties
^^^^^^^^^^^^^^^^
Below is a list of all the properties of solid materials which currently have a standard name associated with them in the IDAES framework.

============================ =================
Variable                     Standard Name
============================ =================
Min. Fluidization Velocity   velocity_mf
Min. Fluidization Voidage    voidage_mf
Particle Size	             particle_dia
Pore Size	             pore_dia
Porosity	             particle_porosity
Specific Surface Area	     area_{basis}
Sphericity	             sphericity
Tortuosity	             tort
Voidage	                     bulk_voidage
============================ =================

Naming Examples
^^^^^^^^^^^^^^^
Below are some examples of the IDAES naming convention in use.

============================== ===========================================================
Variable Name                  Meaning
============================== ===========================================================
enth                           Specific enthalpy of the entire mixture (across all phases)
flow_comp["H2O"]               Total flow of H2O (across all phases)
entr_phase["liq"]              Specific entropy of the liquid phase mixture
conc_phase_comp["liq", "H2O"]  Concentration of H2O in the liquid phase
temperature_red                Reduced temperature
pressure_crit                  Critical pressure
============================== ===========================================================


