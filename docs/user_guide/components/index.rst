Components
==========

The purpose of this section of the documentation is to provide a general introduction to the
components of the IDAES modeling framework. Each component is described in greater detail with 
a link in their description.

.. toctree::
    :glob:
    :hidden:
    
    *
    dmf/index

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
:ref:`Property packages<user_guide/components/property_package:Property Package>`  
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
Property Blocks
^^^^^^^^^^^^^^^

Control Volumes
^^^^^^^^^^^^^^^
:ref:`Control Volumes<user_guide/components/control_volume:Control Volume>` 
are a core component of the IDAES modeling framework and are an integral part
of representing Unit models. 

State Blocks
^^^^^^^^^^^^

Reaction Blocks
^^^^^^^^^^^^^^^

Component objects
^^^^^^^^^^^^^^^^^

Phase objects
^^^^^^^^^^^^^

Pyomo components
----------------
