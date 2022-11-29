Property Metadata Classes
=========================

All property packages (both for thermophysical and reaction properties) are expected to define a set of metadata that contains the following:

1. A `UnitSet` object that records the units of measurement used within the property package.
2. A `PropertySet` object that describes what properties are supported by the property package, how these are constructed, and what dependencies this property package may have (primarily for reaction packages that depend on thermophysical properties from another package).

More information on `UnitSets` and `PropertySets` can be found through the links below.

.. toctree::
    :maxdepth: 1
    
    uom
    property_set

Property Class Metadata
-----------------------

Each property package should create an instance of a `PropertyClassMetadata` object to contain the necessary metadata. This class contains methods that are used hold and interact with the `UnitSet` and `PropertySet` for a specific property package.

.. module:: idaes.core.base.property_meta

.. autoclass:: PropertyClassMetadata
  :members:
  :undoc-members:
 
