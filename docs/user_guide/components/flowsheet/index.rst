Flowsheet
=========

.. toctree::
    :glob:
    :hidden:
    
    *

Flowsheet models are the top level of the IDAES modeling hierarchy. The 
flowsheet is implemented with a 
:ref:`FlowsheetBlock <technical_specs/core/flowsheet_block:Flowsheet Block>`,
which provides a container for other components. Flowsheet models generally contain 
three types of components:

1. :ref:`Unit models<user_guide/components/unit_model/index:Unit Model>`, representing unit operations
2. :ref:`Property packages<user_guide/components/property_package/index:Property Package>`, representing the parameters and relationships for property calculations
3. Arcs, representing connections between unit models

The FlowsheetBlock is also where the
:ref:`time domain<user_guide/components/flowsheet/time_domain:Time Domain>`
is implemented. While the time domain is essential for dynamic modeling,
the time domain exists even for steady state models (single point in time).

Flowsheet models may also contain additional constraints relating to how different unit models 
behave and interact, such as control and operational constraints. Generally speaking, if a 
constraint is purely internal to a single unit, and does not depend on information from other 
units in the flowsheet, then the constraint should be placed inside the relevant unit model. 
Otherwise, the constraint should be placed at the flowsheet level.





