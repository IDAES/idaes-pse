Unit Model Class
================

The UnitModelBlock is class is designed to form the basis of all IDAES Unit Models, and contains a number of methods which are common to all Unit Models.

UnitModelBlock Construction Arguments
--------------------------------------

The UnitModelBlock class by default has only one construction argument, which is listed below. However, most models inheriting from UnitModelBlock should declare their own set of configuration arguments which contain more information on how the model should be constructed.

* dynamic - indicates whether the Unit model should be dynamic or steady-state, and if dynamic = True, the unit is declared to be a dynamic model. dynamic defaults to useDefault if not provided when instantiating the Unit model (see below for more details). It is possible to declare steady-state Unit models as part of dynamic Flowsheets if desired, however the reverse is not true (cannot have dynamic Unit models within steady-state Flowsheets).

Collecting Time Domain
----------------------

The next task of the UnitModelBlock class is to establish the time domain for the unit by collecting the necessary information from the parent Flowsheet model. If the dynamic construction argument is set to `useDefault` then the Unit model looks to its parent model for the dynamic argument, otherwise the value provided at construction is used.

Finally, if the Unit model has a construction argument named "has_holdup" (not part of the base class), then this is checked to ensure that if dynamic = True then has_holdup is also True. If this check fails then a ConfigurationError exception will be thrown.

Modeling Support Methods
------------------------

The UnitModelBlock class also contains a number of methods designed to facilitate the construction of common components of a model, and these are described below.

Build Inlets Method
^^^^^^^^^^^^^^^^^^^

All (or almost all) Unit Models will have inlets and outlets which allow material to flow in and out of the unit being modeled. In order to save the model developer from having to write the code for each inlet themselves, UnitModelBlock contains a method named `build_inlet_port` which can automatically create an inlet to a specified ControlVolume block (or linked to a specified StateBlock). The `build_inlet_port` method is described in more detail in the documentation below.

Build Outlets Method
^^^^^^^^^^^^^^^^^^^^

Similar to `build_inlet_port`, UnitModelBlock also has a method named `build_outlet_port` for constructing outlets from Unit models. The `build_outlet_port` method is described in more detail in the documentation below.

Model Check Method
^^^^^^^^^^^^^^^^^^

In order to support the IDAES Model Check tools, UnitModelBlock contains a simple model_check method which assumes a single Holdup block and calls the `model_check` method on this block. Model developers are encouraged to create their own `model_check` methods for their particular applications.

Initialization Routine
^^^^^^^^^^^^^^^^^^^^^^

All Unit Models need to have an initialization routine, which should be customized for each Unit model, In order to ensure that all Unit models have at least a basic initialization routine, UnitModelBlock contains a generic initialization procedure which may be sufficient for simple models with only one Holdup Block. Model developers are strongly encouraged to write their own initialization routines rather than relying on the default method.

UnitModelBlock Classes
-----------------------

.. module:: idaes.core.base.unit_model

.. autoclass:: UnitModelBlockData
    :members:

.. autoclass:: UnitModelBlock
    :members:
