Property Package
================

.. contents:: :local:

Overview
--------

.. toctree::
    :glob:
    :hidden:

    *
    */index

Property packages provide the relationships and parameters necessary to determine the
properties of process streams. They may be general in purpose, such as ideal gas
equations, or specific to a certain application. Property packages are separated into two categories:

* physical and transport properties
* chemical reaction properties

While several standard property packages are provided in the IDAES model libraries, many process
modeling applications will require specific property packages. Information on developing custom
property packages is provided in the
:ref:`advanced user guide<how_to_guides/custom_models/property_package_development:Custom Property Packages>`.

Since the effort to develop a custom property package is substantial, the IDAES modeling
framework provides a
:ref:`Modular Property Package Framework<explanations/components/property_package/general/index:Modular Property Package Framework>`
and :ref:`Generic Reaction Package Framework<explanations/components/property_package/general_reactions/index:Generic Reaction Package Framework>`
to make it easier to create a package for common property and reaction models.

Physical properties
-------------------

Almost all process models depend on physical properties to some extent, such as
calculation of specific enthalpy or internal energy for energy balances. These properties
only depend on the material being considered and are independent of the unit operations in
that they are used. As such, physical property calculations can be separated from the
unit model calculations and treated as a separate submodel that is called by the unit model.
Each unit model can then create instances of these submodels as required to calculate those
properties required by each unit.

Within IDAES, this is handled by StateBlock objects – these are self-contained submodels
containing the calculations for all necessary thermophysical properties for a given material
at a given point in space and time. IDAES UnitModels create instances of these StateBlocks
wherever they need to calculate physical properties and link to variables within the
StateBlock within the unit model constraints.

However, physical property calculations depend on a set of parameters that are specific
to a given material or mixture. Thus, each instance of a StateBlock for a material use the same
set of parameters. To avoid duplicating these parameters in every instance of a StateBlock for
a given material, these parameters are instead grouped in a PhysicalParameterBlock for that
material that the StateBlocks link to. In this way, there is a single common location for all
parameters.

In summary, physical property packages consist of two parts:

* :ref:`PhysicalParameterBlocks<explanations/components/property_package/physical_param:Physical Parameter Block>`, which contain a set of parameters associated with the specific material(s) being modeled
* :ref:`StateBlocks<explanations/components/property_package/state_block:State Block>`, which contain the actual calculations of the state variables and functions

Reaction properties
-------------------

Reaction property packages represent a collection of calculations necessary to determine the
reaction behavior of a mixture at a given state. Reaction properties depend upon the state and
physical properties of the material, and thus must be linked to a StateBlock that provides the
necessary state and physical property information.

Reaction property packages consist of two parts:

* :ref:`ReactionParameterBlocks<explanations/components/property_package/reaction_param:Reaction Parameter Block>`, which contain a set of parameters associated with the specific reaction(s) being modeled, and
* :ref:`ReactionBlocks<explanations/components/property_package/reaction_block:Reaction Block>`, which contain the actual calculations of the reaction behavior.

Property Metadata
-----------------

One of the goals of the IDAES CMF is to provide users with the maximum degree of flexibility when defining their models. One of the most important roles property packages play within the modeling framework is to define the units of measurement that will be used for those models that use the property packages. Any variable that is created in a unit model will derive its units of measurement from those defined in the associated property package in order to ensure consistency of units. Secondly, the IDAES CMF aims to allow users to only define calculations for those thermophysical properties that will actually be used in their process model, rather than needing to define a comprehensive set of calculations for all possible properties.

In order to do this, each property package is expected to contain a set of metadata that defines the units of measurement used by that package and to list which thermophysical properties are supported by the package (and conversely, which properties are not supported). A discussion of how metadata is defined and handled can be :ref:`found here<reference_guides/core/property_metadata:Property Metadata Classes>`.

Component and Phase Objects
---------------------------

Property packages also rely on component and phase objects.

:ref:`Component Objects<explanations/components/property_package/comp:Component Object>`
are used to identify the chemical species of interest
in a property package and to contain information describing the behavior of that component
(such as properties of that component).

:ref:`Phase Objects<explanations/components/property_package/phase:Phase Object>`
are used to identify the thermodynamic phases of
interest in a property package and to contain information describing the behavior of that phase
(for example the equation of state that describes that phase).

As Needed Properties
--------------------

Process flow sheets often require a large number of properties to be calculate, but not all of
these are required in every unit operation. Calculating additional properties that are not
required is undesirable, as it leads to larger problem sizes and unnecessary complexity of the
resulting model.

To address this, IDAES supports "as needed" construction of properties,
where the variables and constraints required to calculate a given quantity are not added to a
model unless the model calls for this quantity. To designate a property as an "as needed"
quantity, a method can be declared in the associated property BlockData class (StateBlockData or
ReactionBlockData) that contains the instructions for constructing the variables and
constraints associated with the quantity (rather than declaring these within the BlockData's
build method). The name of this method can then be associated with the property via the
add_properties metadata in the property packages ParameterBlock, which indicates that when
this property is called for, the associated method should be run.

The add_properties metadata can also indicate that a property should always be present
(i.e. constructed in the BlockData's build method) by setting the method to `None`, or that it is
not supported by setting the method to `False`.

Modular Property Package Framework
----------------------------------

Property packages represent the core of any process model, and having a suitable property
package is key to successfully modeling any process system. However, developing property
packages is a significant challenge even for experienced modelers as they involve large numbers
of tightly coupled constraints and parameters. The
:ref:`Modular Property Package Framework<explanations/components/property_package/general/index:Modular Property Package Framework>`
was designed to help users build property packages with the least effort possible by levarging libraries
of modular sub-models that include common types of property calculations.


Generic Reaction Package Framework
----------------------------------

Similar to the Generic Property Package Framework, the
:ref:`Generic Reaction Package Framework<explanations/components/property_package/general_reactions/index:Generic Reaction Package Framework>`
helps users create reaction property packages for common systems.

