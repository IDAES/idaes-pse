Property Packages
=================

.. toctree::
    :maxdepth: 1

    state_block
    reaction_block

IDAES Property Packages
-----------------------

The IDAES process modeling framework divides property calculations into two parts;

* physical and transport properties
* chemical reaction properties

Defining the calculations to be used when calculating properties is done via "property packages", which contain a set of related calculations for a number of properties of interest. Property packages may be general in purpose, such as ideal gas equations, or specific to a certain application.

As Needed Properties
--------------------

Process flow sheets often require a large number of properties to be calculate, but not all of these are required in every unit operation. Calculating additional properties that are not required is undesirable, as it leads to larger problem sizes and unnecessary complexity of the resulting model.

To address this, the IDAES modeling framework supports "as needed" construction of properties, where the variables and constraints required to calculate a given quantity are not added to a model unless the model calls for this quantity. To designate a property as an "as needed" quantity, a method can be declared in the associated property BlockData class (StateBlockData or ReactionBlockData) which contains the instructions for constructing the variables and constraints associated with the quantity (rather than declaring these within the BlockData's build method). The name of this method can then be associated with the property via the add_properties metadata in the property packages ParameterBlock, which indicates to the framework that when this property is called for, the associated method should be run.

The add_properties metadata can also indicate that a property should always be present (i.e. constructed in the BlockData's build method) by setting the method to None, or that it is not supported by setting the method to False.
