Property Sets
=============

Due to the large number of potential thermophysical properties that are used in process applications and the complexity associated with these calculations, one of the key aspects of the IDAES CMF is to be able to only specify calculations for those properties that are actually required for a given problem. In order to support this, and allow users to understand what is supported by a given property package, each property package is expected to define a `PropertySet` that lists what properties are supported and how these are constructed.

`PropertySets` are application specific, and are intended to contain a comprehensive list of all properties of interest for the given application. Each property package should select a `PropertySet` appropriate for its application, and can then specify the sub-set of properties in the `PropertySet` are supported by that package, and how they will be constructed. The `PropertySet` will also define the expected units of measurement for each property by linking to the `UnitSet` defined for the property package.

The code below shows an example of how to set up a property package using the `StandardPropertySet` (suitable for most liquid and vapor phase property packages) and identify which properties are supported.

.. code-block:: python

    @classmethod
    def define_metadata(cls, obj):
        # Select StandardPropertySet
        obj.define_property_set(StandardPropertySet)
        # Update the flow mol property to indicate that it is supported
        obj.properties.FLOW_MOL.update_property(supported=True)

.. note::

  Any property defined in a `PropertySet` that is not explicitly marked as `supported` is assumed to not be supported (i.e., the default value of `supported` is `False`).

Core Property Sets
------------------

A list of the existing `PropertySets` defined in the IDAES CMF is shown below:

.. toctree::
  :maxdepth: 1

  standard_property_set
  electrolyte_property_set

Property Metadata Object
------------------------

Each individual property within a `PropertySet` is defined using a `PropertyMetadata` object, which records information such as whether the property is supported by the given property package, and what method is used to construct this property (if applicable).

.. module:: idaes.core.base.property_set

.. autoclass:: PropertyMetadata
  :members:
  :undoc-members:

Property Set Base Class
-----------------------

The `PropertySetBase` class defines common methods associated with `PropertySets`, and should be used as the base class for defining custom `PropertySets`.

.. autoclass:: PropertySetBase
  :members:
  :undoc-members:
