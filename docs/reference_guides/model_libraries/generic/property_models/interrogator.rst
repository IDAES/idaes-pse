Property Interrogator Tool
==========================

When preparing to model a process flowsheet, it is necessary to specify models for all the thermophysical and kinetic properties that will be required by the different unit operations to simulate the process. However, it is often difficult to know what properties will be required *a priori*. The IDAES Property Interrogator tool allows a user to define a general flowsheet structure and interrogate it for the full list of properties that will be required, thus informing them of what methods they will need to define in their property package(s).

Tool Usage
----------

The IDAES Properties Interrogator tool consists of two classes; a `PropertiesInterrogatorBlock` and a `ReactionInterrogatorBlock`. These blocks are used in place of the normal `PhysicalParameterBlock` and `ReactionParameterBlock` whilst declaring a flowsheet, however rather than constructing a solvable flowsheet they record all calls for properties made whilst constructing the flowsheet. These Blocks then contain a number of methods for reporting the logged property calls for the user.

An example of how Property Interrogator tool is used is shown below:

.. testcode::

    import pyomo.environ as pyo  # Pyomo environment
    from idaes.core import FlowsheetBlock
    from idaes.models.unit_models import CSTR
    from idaes.models.properties.interrogator import PropertyInterrogatorBlock, ReactionInterrogatorBlock

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": True, "time_units": pyo.units.s})

    m.fs.params = PropertyInterrogatorBlock()
    m.fs.rxn_params = ReactionInterrogatorBlock(
            default={"property_package": m.fs.params})

    m.fs.R01 = CSTR(default={"property_package": m.fs.params,
                             "reaction_package": m.fs.rxn_params,
                             "has_heat_of_reaction": True})

.. note::
    Flowsheets constructed using the Property Interrogator tools are not solvable flowsheets, and will result in errors if sent to a solver.

Output and Display Methods
--------------------------

Both the `PropertiesInterrogatorBlock` and `ReactionInterrogatorBlock` support the following methods for reporting the results of the flowsheet interrogation. The `PropertiesInterrogatorBlock` will contain a summary of all thermophysical properties expected of a `StateBlock` in the flowsheet, whilst the `ReactionInterrogatorBlock` will contain a summary of all reaction related properties required of a `ReactionBlock`.

* list_required_properties()  - returns a list containing all properties called for by the flowsheet.
* print_required_properties() - prints a summary of the required properties
* list_models_requiring_property(property)  - returns a list of unit models within the flowsheet that require the given property
* print_models_requiring_property(property)  - prints the name of all unit models within the flowsheet that require the given property
* list_properties_required_by_model(model)  - returns a list of all properties required by a given unit model in the flowsheet
* print_properties_required_by_model(model)  - prints a summary of all properties required by a given unit model in the flowsheet

For more details on these methods, see the detailed class documentation below.

Additionally, the `PropertiesInterrogatorBlock` and `ReactionInterrogatorBlock` contain a `dict` named `required_properties` which stores the data regarding the properties required by the model. The keys of this `dict` are the names of all the properties required (as strings) and the values are a list of names for the unit models requiring the given property.


Class Documentation
-------------------

.. currentmodule:: idaes.models.properties.interrogator.properties_interrogator

.. autoclass:: PropertyInterrogatorBlock
   :members:

.. autoclass:: PropertyInterrogatorData
   :members:

.. currentmodule:: idaes.models.properties.interrogator.reactions_interrogator

.. autoclass:: ReactionInterrogatorBlock
   :members:

.. autoclass:: ReactionInterrogatorData
   :members:

