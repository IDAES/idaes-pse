Unit Sets
=========

All property packages within IDAES are expected to define a metadata class as part of the package's ParameterBlock, which among other things contains a definition of the base units of measurement used by that property package.  The definition of units of measurement is held within a `UnitSet` object that defines a number of common quantities of interest in process applications.

An example of defining the default units for a property package is shown below.

.. code-block:: python

    from pyomo.environ import units

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units({'time': units.s,
                               'length': units.m,
                               'mass': units.kg,
                               'amount': units.mol,
                               'temperature': units.K})

Each property package can define a default units for 7 base quantities listed below along with the default units assumed if not provided by the property package developer.

* time (seconds)
* length (meters)
* mass (kilograms)
* amount of substance (mole)
* temperature (Kelvin)
* current (Ampere)
* luminous intensity (candela)

Units must be defined using Pyomo's Units container (`from pyomo.environ import units`).

Units of measurement are stored in a `UnitSet` object that uses these base units to derive the units of all other quantities that will be used in the property model and any unit operation that uses that property package. In order to avoid the need for unit conversion in all expressions, the IDAES Framework assumes that all quantities within the property package are based on the chosen set of base units. Parameters and correlations may be based on different sets of unit as necessary (e.g., from literature sources using different base units), however the final quantity must be converted to the set of base units defined in the metadata.

Retrieving Units of Measurement
-------------------------------

Modelers may retrieve the units of measurement for any defined quantity by looking these up in the `UnitSet` associated with the property package. The `UnitSet` may be accessed using the following code:

.. code-block:: python

    unit_set = property_package.get_metadata().default_units

Modelers can then retrieve the units for a given quantity by looking up the quantity on the `UnitSet` as shown below (using time as an example):

.. code-block:: python

    time_units = unit_set.TIME

.. note::
    For backwards-compatibility, units can also be retrieved using dict-like notation, i.e., units_set["time"]. However, this feature may be
    deprecated in the future, thus modelers are encouraged to use the attribute-based approach to access units.

A full list of defined quantities is shown below.

UnitSet Class
-------------

.. module:: idaes.core.base.property_meta
  :noindex:

.. autoclass:: UnitSet
  :members:
  :undoc-members:
  
