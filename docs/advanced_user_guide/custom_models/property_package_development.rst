Custom Property Model
=====================

.. contents:: :local:

.. warning:: This section is currently being developed

Property Packages
    • 3 classes
    • Build on demand
    • ParameterBlock
        ◦ base class
        ◦ config block
        ◦ decorator not necessary
        ◦ build
        ◦ Components
        ◦ Phases
        ◦ Parameters
        ◦ Metadata
            ▪ Units
            ▪ Add Properties
    • PropertyBlock (need better name)
        ◦ no decorator
        ◦ base class
        ◦ initialize
        ◦ release_state
    • PropertyBlockData
        ◦ special use of decorator
        ◦ base class
        ◦ config block
        ◦ build
        ◦ state vars
        ◦ properties
        ◦ Framework methods
            ▪ get_material_flow_basis
            ▪ get_material_flow_terms
            ▪ get_material_density_terms
            ▪ get_material_diffusion_terms
            ▪ get_enthalpy_flow_terms
            ▪ get_energy_density_terms
            ▪ get_energy_diffusion_terms
            ▪ default_material_balance_type
            ▪ default_energy_balance_type
            ▪ define_state_vars
            ▪ define_port_members
            ▪ define_display_vars

Tutorials
---------
Tutorials demonstrating how to create custom property models are found
:ref:`here<advanced_user_guide/learning_materials/property_tutorials/index:Property Model Tutorials>`.       
    
