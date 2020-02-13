Property Interrogator Tool
==========================

When preparing to model a process flowsheet, it is necessary to specify models for all the thermophysical and kinetic properties that will be required by the differnt unit operations to simulate the process. However, it is often difficult to know what properties will be required *a priory*. The IDAES Property Interrogator tool allows a user to define a general flowsheet structure and interrogate it for the full list of properties that will be required, thus informing them of what methods they will need to define in their property package(s).

Tool Usage
----------

The IDAES Properties Interrogator tool consists of two classes; a `PropertiesInterrogatorBlock` and a `ReactionInterrogatorBlock`.
