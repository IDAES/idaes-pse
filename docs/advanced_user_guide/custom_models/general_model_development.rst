Developing Custom Models
========================

.. contents:: :local:

It is difficult to build a single model library that will suit all modeling needs, and thus users will inevitably encounter situations where they need to create new models to represent their processes. The IDAES Process Modeling Framework has been developed with this in mind, and all components of the framework have been designed to be fully accessible and modifiable. This section of the documentation will explain how users can develop new models, or modify existing models, for use within the IDAES modeling environment.

Creating New Modeling Components
--------------------------------

All models within the IDAES Process Modeling Framework, be they models of unit operations of thermodynamic properties, are constructed in the same way and user defined models follow the same structure. Each model component is a set of instructions on how to assemble a Pyomo ‘Block’ containing the necessary variables, expressions and constraints to describe the desired process. These instructions are contained within Python `classes` which can be written and modified by the user.

Details on what is required when constructing custom models of different types will be provided in subsequent sections of this documentation, whoever there are some steps common to all types of models which will be discussed here.

Defining New Model Classes
--------------------------

There is a significant amount of work that has to be done behind the scenes in order to create a new model component and integrate it into the IDAES modeling framework. Rather than force the user to understand and implement this themselves, IDAES provides a number of standard tools which users can use to automate this when creating new models. These tools are used when declaring new components, and should be used in most cases when a user is creating a new model (exceptions are mentioned in the detailed documentation for each type of model).

An example of declaring  new model (named NewModel) is shown below:

.. code-block:: python

    @declare_process_block_class("NewModel")
    class NewModelData(BaseClass):

The first part of the declaration is the `declare_process_block_class` decorator, which automates all the code required to create a new type of Pyomo `Block`. This decorator needs to be provided with a name for the new model (NewModel). Understanding the details of class decoration within Python and the function of the `declare_process_block_class` decorator are not necessary for developing new models, however users who wish to read more can see the :ref:`technical specifications<technical_specs/core/process_block:ProcessBlock Class>`.

Inheriting from Existing Models
-------------------------------

In addition to inheriting from the IDAES base classes, it is also possible to inherit from existing models, which provides an easy way to build off an existing model instead of starting afresh. When inheriting from an existing model, you gain access to all the instructions written for constructing that model and can then add to or modify those instructions as needed. All of the models in the core IDAES model library were written with this in mind, and were designed to provide a core model representing the simplest representation of each piece of process equipment possible to allow users to easily build upon these as a foundation.

A useful concept when modifying existing models through inheritance is “overloading” of methods. Any method defined by the inherited class can be overloaded and replaced by a new method of the same name defined in the new class. Thus, it is possible to selectively modify and replace parts of the existing model if they were defined using methods. A case where this may be useful is for cases where an existing model meets most of a user’s needs, but the user wishes to use a difference equation for efficiency (for example). If the existing model defined a method specifically for writing the efficiency constraint, then this can easily be replaced by inheriting the existing model and writing a new method for efficiency with the desired equation. This will overload the method in the original model resulting in the new model using the desired equation instead whilst requiring little effort on the part of the user (this does require the original model to use modular methods for each performance equation however).

Config Blocks
-------------

The `build` Method
------------------

.. toctree::
    :maxdepth: 1
    
    unit_model_development
    property_package_development
    reaction_package_development
  

