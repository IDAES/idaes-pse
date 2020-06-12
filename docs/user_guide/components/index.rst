Components
==========

The purpose of this section of the documentation is to provide a general introduction to the
components of the IDAES modeling framework. Each component is described in greater detail with 
a link in their description.

.. toctree::
    :glob:
    :hidden:
    
    *
    */index

.. contents:: :local:


High-level components
---------------------
Flowsheet
^^^^^^^^^
:ref:`Flowsheet models<user_guide/components/flowsheet:Flowsheet>`
are the top level of the IDAES modeling framework. Flowsheet models represent 
traditional process flowsheets, containing a number of unit models connected together into a 
flow network and the property packages.

Unit Models
^^^^^^^^^^^
:ref:`Unit models<user_guide/components/unit_model:Unit Model>` 
represent individual pieces of equipment and their processes on the IDAES platform.

Property Packages
^^^^^^^^^^^^^^^^^
:ref:`Property packages<user_guide/components/property_model/index:Property Package>`  
represent the physical, thermodynamic, and reactive properties of the process streams on the 
IDAES platform.

Time Domain
^^^^^^^^^^^
In order to support dynamic modeling, the IDAES modeling framework has a 
:ref:`time domain<user_guide/components/time_domain:Time Domain>`. 
This time domain exists even for steady state models (single point in time).

Data Management Framework
^^^^^^^^^^^^^^^^^^^^^^^^^
The :ref:`Data Management Framework <user_guide/components/dmf/index:Data Management Framework>`
is used to manage all the data needed by the IDAES framework, including flowsheets, models, 
and results. It stores metadata and data in persistent storage.

Low-level components
---------------------
Control Volumes
^^^^^^^^^^^^^^^
:ref:`Control Volumes<user_guide/components/control_volume:Control Volume>` 
are a core component of the IDAES modeling framework and are an integral part
of representing Unit models. 

Physical Parameter Blocks
^^^^^^^^^^^^^^^^^^^^^^^^^
:ref:`Physical Parameter Blocks<user_guide/components/property_model/physical_param:Physical Parameter Block>`
are a core component of physical property packages. They serve as a central location for 
containing all the parameters and indexing sets used by a given physical property package.

State Blocks
^^^^^^^^^^^^
:ref:`State Blocks<user_guide/components/property_model/state_block:State Block>`
are a core component of physical property packages. They are used within all IDAES Unit models 
(generally within ControlVolume Blocks) in order to calculate physical properties given the 
state of the material. State Blocks are notably different to other types of Blocks within IDAES 
as they are always indexed by time (and possibly space as well).

Reaction Parameter Blocks
^^^^^^^^^^^^^^^^^^^^^^^^^
:ref:`Reaction Parameter Blocks<user_guide/components/property_model/reaction_param:Reaction Parameter Block>`
are a core component of reaction property packages. They serve as a central location for 
containing all the parameters and indexing sets used by a given reaction property package.

Reaction Blocks
^^^^^^^^^^^^^^^
:ref:`Reaction Blocks<user_guide/components/property_model/reaction_block:Reaction Block>`
are a core component of reaction property packages. They are used 
within IDAES Unit models (generally within ControlVolume Blocks) in order to calculate reaction 
properties given the state of the material (provided by an associated StateBlock). Reaction 
Blocks are notably different to other types of Blocks within IDAES as they are always indexed by time 
(and possibly space as well), and are also not fully self contained 
(in that they depend upon the associated state block for certain variables). 

Component Objects
^^^^^^^^^^^^^^^^^
:ref:`Component Objects<user_guide/components/comp:Component Object>` 
are used in the IDAES Process Modeling Framework to identify the chemical species of interest 
in a property package and to contain information describing the behavior of that component 
(such as properties of that component).

Phase objects
^^^^^^^^^^^^^
:ref:`Phase Objects<user_guide/components/phase:Phase Object>`
are used in the IDAES Process Modeling Framework to identify the thermodynamic phases of 
interest in a property package and to contain information describing the behavior of that phase 
(for example the equation of state which describes that phase).

Pyomo components
----------------
Pyomo provides the basis for the IDAES modeling framework. Pyomo components include: Sets, 
Parameters, Variables, Objectives, Constraints, Expressions, Suffixes as detailed
`here <https://pyomo.readthedocs.io/en/stable/pyomo_modeling_components/index.html>`_. Pyomo also
includes other important components such as Arcs and Ports, which are used to connect unit models
in flowsheets. A new user of IDAES would benefit from reading the 
`Pyomo documentation <https://pyomo.readthedocs.io/en/stable/index.html>`_.