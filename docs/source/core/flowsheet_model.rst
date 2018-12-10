Flowsheet Models
================

.. contents:: Contents 
    :depth: 2

Flowsheet models make up the top level of the IDAES modeling framework, and represent the flow of material and energy through a process. Flowsheets will generally contain a number of UnitModels to represent unit operations within the process, and will contain one or more Property Packages which represent the thermophysical and transport properties of material within the process.

Flowsheet models are responsible for establishing and maintaining the time domain of the model, including declaring whether the process model will be dynamic or steady-state. This time domain is passed on to all models attached to the flowsheet (such as Unit Models and sub-Flowsheets). The Flowsheet model also serves as a centralized location for organizing property packages, and can set one property package to use as a default throughout the flowsheet.

Flowsheet Blocks may contain other Flowsheet Blocks in order to create nested flowsheets and to better organize large, complex process configurations. In these cases, the top-level Flowsheet Block creates the time domain, and each sub-flowsheet creates a reference this time domain. Sub-flowsheets may make use of any property package declared at a higher level, or declare new property package for use within itself - any of these may be set as the default property package for a sub-Flowsheet.

Default Property Packages
-------------------------

Flowsheet Blocks may assign a property package to use as a default for all UnitModels within the Flowsheet. If a specific property package is not provided as an argument when constructing a UnitModel, the UnitModel will search up the model tree until it finds a default property package declared. The UnitModel will use the first default property package it finds during the search, and will return an error if no default is found.

Flowsheet Configuration Arguments
---------------------------------

Flowsheet blocks have three configuration arguments which are stored within a Config block (flowsheet.config). These arguments can be set by passing arguments when instantiating the class, and are described below:

* dynamic - indicates whether the flowsheet should be dynamic or steady-state. If dynamic = True, the flowsheet is declared to be a dynamic flowsheet, and the time domain will be a Pyomo ContunuousSet. If dynamic = False, the flowsheet is declared to be steady-state, and the time domain will be an ordered Pyomo Set. For top level Flowsheets, dynamic defaults to False if not provided. For lower level Flowsheets, the dynamic will take the same value as that of the parent model if not provided. It is possible to declare steady-state sub-Flowsheets as part of dynamic Flowsheets if desired, however the reverse is not true (cannot have dynamic Flowsheets within steady-state Flowsheets).
* time_set - use to initialize the time domain in top-level Flowsheets. When constructing the time domain in top-level Flowsheets, time_set is used to initialize the ContinuousSet or Set created. This can be used to set start and end times, and to establish points of interest in time (e.g. times when disturbances will occur). If dynamic = True, time_set defaults to [0.0, 1.0] if not provided, if dynamic = False time_set defaults to [0.0]. time_set is not used in sub-Flowsheets and will be ignored.
* default_property_package - can be used to assign the default property package for a Flowsheet. Defaults to None if not provided.

Flowsheet Classes
-----------------

.. module:: idaes.core.flowsheet_model

.. autoclass:: FlowsheetBlockData
    :members:

.. autoclass:: FlowsheetBlock
    :members:
