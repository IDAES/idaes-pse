Developing Custom Models
========================

It is difficult to build a single model library that will suit all modeling needs, and thus users will inevitably encounter situations where they need to create new models to represent their processes. The IDAES Process Modeling Framework has been developed with this in mind, and all components of the framework have been designed to be fully accessible and modifiable. This section of the documentation will explain how users can develop new models, or modify existing models, for use within the IDAES modeling environment.

.. contents:: :local:

Creating New Modeling Components
--------------------------------

All models within the IDAES Process Modeling Framework, be they models of unit operations of thermodynamic properties, are constructed in the same way and user defined models follow the same structure. Each model component is a set of instructions on how to assemble a Pyomo ‘Block’ containing the necessary variables, expressions and constraints to describe the desired process. These instructions are contained within Python `classes` which can be written and modified by the user.

Details on what is required when constructing custom models of different types will be provided in subsequent sections of this documentation, however there are some steps common to all types of models which will be discussed here.

Defining New Model Classes
--------------------------

There is a significant amount of work that has to be done behind the scenes in order to create a new model component and integrate it into the IDAES modeling framework. Rather than force the user to understand and implement this themselves, IDAES provides a number of standard tools which users can use to automate this when creating new models. These tools are used when declaring new components, and should be used in most cases when a user is creating a new model (exceptions are mentioned in the detailed documentation for each type of model).

An example of declaring  new model (named NewModel) is shown below:

.. code-block:: python

    @declare_process_block_class("NewModel")
    class NewModelData(BaseClass):

The first part of the declaration is the `declare_process_block_class` decorator, which automates all the code required to create a new type of Pyomo `Block`. This decorator needs to be provided with a name for the new model (NewModel). Understanding the details of class decoration within Python and the function of the `declare_process_block_class` decorator are not necessary for developing new models, however users who wish to read more can see the :ref:`technical specifications<technical_specs/core/process_block:ProcessBlock Class>`.

The second part of the declaration creates a `NewModelData` class which inherits from an existing `BaseClass`. The `NewModelData` class needs to contain the instructions necessary for building the desired model and must be populated by the user. In practice, each type of model (unit operation, thermophysical properties, etc.) has a set of common tasks which appear in most models of that type. To assist users with developing new model and reduce the need for them to rewrite these common tasks, IDAES provides a set of base classes which contain many of these common instructions. Model developers can use these base classes as a foundation for their new models using what is referred to as “inheritance”, as shown in the example above. In doing so, the new model class automatically gains access to all of the common instructions in the base class which can be used in constructing the new model.

Inheriting from Existing Models
-------------------------------

In addition to inheriting from the IDAES base classes, it is also possible to inherit from existing models, which provides an easy way to build off an existing model instead of starting afresh. When inheriting from an existing model, you gain access to all the instructions written for constructing that model and can then add to or modify those instructions as needed. All of the models in the core IDAES model library were written with this in mind, and were designed to provide a core model representing the simplest representation of each piece of process equipment possible to allow users to easily build upon these as a foundation.

A useful concept when modifying existing models through inheritance is “overloading” of methods. Any method defined by the inherited class can be overloaded and replaced by a new method of the same name defined in the new class. Thus, it is possible to selectively modify and replace parts of the existing model if they were defined using methods. For example, suppose there is an existing model that meets most of a user's needs, but the user would like to use a different equation for efficiency. If the existing model defined a method specifically for writing the efficiency constraint, then this can be replaced by inheriting the existing model and writing a new method for efficiency with the desired equation. This will overload the method in the original model, creating a new model which uses the desired equation. This requires little effort on the part of the user, but does require the original model to use modular methods for each performance equation however.

Config Blocks
-------------

Whilst the model class contains the instructions necessary to build a model object, it is often necessary to provide additional information when creating an instance of a model. One example of this is informing a unit model of which property package to use for a given instance. When creating a new model class, it is necessary to define the information that a user may pass to the class when creating an instance of the new model, which is done using configuration blocks (config blocks for short) – this is where the information in the “default” keyword is sent when an instance of a model is created.

Configuration blocks are defined by declaring a `CONFIG` object for each new model data class, as shown in the example below. The `CONFIG` object should be an instance of a Pyomo `ConfigBlock`.

.. code-block:: python

    @declare_process_block_class("NewModel")
    class NewModelData(BaseClass):

    CONFIG = BaseClass.CONFIG()

Each type of model has a set of expected inputs (or arguments) which are determined by the type of model and can be inherited from the appropriate base class (as shown above). Users may also add custom configuration arguments to their models as needed by declaring new entries to the `CONFIG` block as shown below:

.. code-block:: python

    from pyomo.common.config import ConfigValue

    @declare_process_block_class("NewModel")
    class NewModelData(BaseClass):
        CONFIG = BaseClass.CONFIG()
        CONFIG.declare("new_argument", ConfigValue(
            default =  # default value for argument,
            domain =  # condition input must satisfy,
            description = "short description of argument",
            doc = "longer description of argument"))

.. note::

    Configuration arguments are set when an instance of a model is created and are generally only used at build-time. That is, once a model has been constructed changing a configuration argument has no effect on the model structure.

The `build` Method
------------------

Finally, the core of any IDAES model class is the `build` method, which contains the set of instructions to be executed when a model is created. The `build` method acts as the rule for constructing the resulting Pyomo `Block`, and needs to contain the instructions necessary for constructing the variable, expressions and constraints which describe the model. The `build` method is written in Python code and should construct the necessary Pyomo components, and may make use of sub-methods to modularize the model construction.

In almost all cases, the first instruction in a `build` method should be to call the `build` method of the inherited (base) class. This is necessary to execute the instructions in the base class, and can be done with the following line of code:

.. code-block:: python

    super().build()

Types of Models
---------------

.. toctree::
    :maxdepth: 1
    
    unit_model_development
    property_package_development
    reaction_package_development
  

