Property Package
================

.. contents:: :local:

Overview
--------

Property packages provide the relationships and parameters necessary to determine the 
properties of process streams. Property packages may be general in purpose, such as ideal gas 
equations, or specific to a certain application. The IDAES modeling framework divides property 
packages into two parts:

* physical and transport properties
* chemical reaction properties

While the IDAES modeling framework provides several standard property packages, many process
modeling applications will require specific property packages. Information on developing custom
property packages is provided in the 
:ref:`advanced user guide<advanced_user_guide/property_development:Property Model Development>`.
Since the effort to develop a custom property package is substantial, the IDAES modeling
framework provides a 
:ref:`Generic Property Package Framework<user_guide/components/property_model/general/index:Generic Property Package Framework>` 
and :ref:`Generic Property Package Framework<user_guide/components/property_model/general_reactions/index:Generic Reaction Package Framework>`
to make it easier to create a package for common property and reaction models.

Physical properties
-------------------

Almost all process models depend on physical properties to some extent, such as 
calculation of specific enthalpy or internal energy for energy balances. These properties 
only depend on the material being considered and are independent of the unit operations in 
which they are used. As such, physical property calculations can be separated from the 
unit model calculations and treated as a separate submodel which is called by the unit model. 
Each unit model can then create instances of these submodels as required to calculate those 
properties required by each unit.

Within IDAES, this is handled by StateBlock objects – these are self-contained submodels 
containing the calculations for all necessary thermophysical properties for a given material 
at a given point in space and time. IDAES Unit models create instances of these StateBlocks 
wherever they need to calculate physical properties and link to variables within the 
StateBlock within the unit model constraints.

However, physical property calculations depend on a set of parameters which are specific 
to a given material or mixture. Thus, each instance of a StateBlock for a material use the same 
set of parameters. To avoid duplicating these parameters in every instance of a StateBlock for 
a given material, these parameters are instead grouped in a PhysicalParameterBlock for that 
material which the StateBlocks link to. In this way, there is a single common location for all 
parameters.

In summary, physical property packages consist of two parts:

* PhysicalParameterBlocks, which contain a set of parameters associated with the specific material(s) being modeled
* StateBlocks, which contain the actual calculations of the state variables and functions

Reaction properties
-------------------

Reaction property packages represent a collection of calculations necessary to determine the 
reaction behavior of a mixture at a given state. Reaction properties depend upon the state and 
physical properties of the material, and thus must be linked to a StateBlock which provides the 
necessary state and physical property information.

Reaction property packages consist of two parts:

* ReactionParameterBlocks, which contain a set of parameters associated with the specific reaction(s) being modeled, and
* ReactionBlocks, which contain the actual calculations of the reaction behavior.

As Needed Properties
--------------------

Process flow sheets often require a large number of properties to be calculate, but not all of 
these are required in every unit operation. Calculating additional properties that are not 
required is undesirable, as it leads to larger problem sizes and unnecessary complexity of the 
resulting model.

To address this, the IDAES modeling framework supports "as needed" construction of properties, 
where the variables and constraints required to calculate a given quantity are not added to a 
model unless the model calls for this quantity. To designate a property as an "as needed" 
quantity, a method can be declared in the associated property BlockData class (StateBlockData or 
ReactionBlockData) which contains the instructions for constructing the variables and 
constraints associated with the quantity (rather than declaring these within the BlockData's 
build method). The name of this method can then be associated with the property via the 
add_properties metadata in the property packages ParameterBlock, which indicates to the 
framework that when this property is called for, the associated method should be run.

The add_properties metadata can also indicate that a property should always be present 
(i.e. constructed in the BlockData's build method) by setting the method to None, or that it is 
not supported by setting the method to False.

Generic Property Package Framework
----------------------------------
Property packages represent the core of any process model, and having a suitable property 
package is key to successfully modeling any process system. However, developing property 
packages is a significant challenge even for experienced modelers as they involve large numbers 
of tightly coupled constraints and parameters. The IDAES modeling framework provides 
a Generic Property Package Framework to provide a flexible platform on which users can build 
property packages for common types of systems by calling upon libraries of modular sub-models 
to build up complex property calculations with the least effort possible.




