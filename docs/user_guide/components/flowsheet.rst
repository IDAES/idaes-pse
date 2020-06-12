Flowsheet
=========

Flowsheet models are the top level of the IDAES modeling hierarchy. The 
flowsheet is implemented with a 
:ref:`FlowsheetBlock <technical_specs/core/flowsheet_block:Flowsheet Block>`,
which provides a container for other components. Flowsheet models generally contain 
three types of components:

1. :ref:`Unit models<user_guide/components/unit_model:Unit Model>`, representing unit operations
2. :ref:`Property packages<user_guide/components/property_model/index:Property Package>`, representing the parameters and relationships for property calculations
3. Arcs, representing connections between unit models

Flowsheet models are also an integral part of the time domain as described 
:ref:`here<user_guide/components/time_domain:Time Domain>`.

Flowsheet models may also contain additional constraints relating to how different Unit models 
behave and interact, such as control and operational constraints. Generally speaking, if a 
constraint is purely internal to a single unit, and does not depend on information from other 
units in the flowsheet, then the constraint should be placed inside the relevant unit model. 
Otherwise, the constraint should be placed at the Flowsheet level.





