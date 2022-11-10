.. index::
    pair: idaes.core.base.flowsheet_model;FlowsheetBlockData
    pair: idaes.core.base.flowsheet_model;FlowsheetBlock

Flowsheet Block
===============

Default Property Packages
-------------------------

Flowsheet Blocks may assign a property package to use as a default for all UnitModels within the Flowsheet. If a specific property package is not provided as an argument when constructing a UnitModel, the UnitModel will search up the model tree until it finds a default property package declared. The UnitModel will use the first default property package it finds during the search, and will return an error if no default is found.

Flowsheet Configuration Arguments
---------------------------------

Flowsheet blocks have three configuration arguments which are stored within a Config block (flowsheet.config). These arguments can be set by passing arguments when instantiating the class, and are described below:

* dynamic - indicates whether the flowsheet should be dynamic or steady-state. If dynamic = True, the flowsheet is declared to be a dynamic flowsheet, and the time domain will be a Pyomo ContinuousSet. If dynamic = False, the flowsheet is declared to be steady-state, and the time domain will be an ordered Pyomo Set. For top level Flowsheets, dynamic defaults to False if not provided. For lower level Flowsheets, the dynamic will take the same value as that of the parent model if not provided. It is possible to declare steady-state sub-Flowsheets as part of dynamic Flowsheets if desired, however the reverse is not true (cannot have dynamic Flowsheets within steady-state Flowsheets).
* time - a reference to the time domain for the flowsheet. During flowsheet creation, users may provide a Set or ContinuousSet that the flowsheet should use as the time domain. If not provided, then the flowsheet will look for a parent flowsheet and set this equal to the parent's time domain, otherwise a new time domain will be created and assigned here.
* time_units - used to specify the units of the time domain, and must be a Pyomo Unit object (cannot be a compound unit). This is necessary for dynamic flowsheets, but can be neglected in steady-state cases. In cases where the time domain is inherited from a parent flowsheet, the time units will also be inherited.
* time_set - used to initialize the time domain in top-level Flowsheets. When constructing the time domain in top-level Flowsheets, time_set is used to initialize the ContinuousSet or Set created. This can be used to set start and end times, and to establish points of interest in time (e.g. times when disturbances will occur). If dynamic = True, time_set defaults to [0.0, 1.0] if not provided, if dynamic = False time_set defaults to [0.0]. time_set is not used in sub-Flowsheets and will be ignored.
* default_property_package - can be used to assign the default property package for a Flowsheet. Defaults to None if not provided.

Flowsheet Classes
-----------------

.. module:: idaes.core.base.flowsheet_model

.. autoclass:: FlowsheetBlockData
    :members:

.. autoclass:: FlowsheetBlock
    :members:
