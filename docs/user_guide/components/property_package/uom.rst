Defining Units of Measurement
=============================

All property packages within IDAES are expected to define a metadata class as part of the package's ParameterBlock, which amongst other things contains a definition of the base units of measurement used by that property package. An example of defining the default units for a property package is shown below.

.. code-block:: python

    from pyomo.environ import units

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units({'time': units.s,
                               'length': units.m,
                               'mass': units.kg,
                               'amount': units.mol,
                               'temperature': units.K})

Each property package should define a default units for 7 base quantities listed below:

* time
* length
* mass
* amount of substance
* temperature
* current (optional)
* luminous intensity (optional)

Units must be defined using Pyomo's Units container (`from pyomo.environ import units`), and all quantities within the property package must be based on the chosen set of base units. Parameters and correlations may be based on different sets of unit as necessary (e.g. from literature sources using different base units), however the final quantity must be converted to the set of base units defined in the metadata.

