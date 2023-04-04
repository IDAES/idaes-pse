No Method
=========

.. index::
   pair: idaes.models.properties.modular_properties.transport_properties.no_method; NoMethod

.. currentmodule:: idaes.models.properties.modular_properties.transport_properties.no_method

If a transport property, like dynamic viscosity or thermal conductivity, is specified in any phase, the modular
property framework requires it be defined in every phase. However, these properties are often not necessary for
calculations in current unit models. A constant, dummy value should not be used, because then future users may
mistakenly believe that an appropriate method and value have been chosen and obtain erroneous results.
The :code:`NoMethod` method for transport properties has been implemented as a solution to this problem. It returns
an `Expression.Skip`, so a user who attempts to call for the (previously unused) property is returned a
`KeyError` during model construction. The method can be imported from :code:`idaes.models.properties.modular_properties.transport_properties`.

Example
-------

The code snippet below demonstrates how to specify use of the :code:`NoMethod` method. Note that users need only use this method
when they are using the same property (with an actual method of calculation) in another phase.

.. code-block:: python

  from idaes.models.properties.modular_properties.transport_properties import ThermalConductivityWMS, NoMethod

  configuration = {
    "phases":{
      "Vap": {
        "type": VaporPhase,
        ...
        "therm_cond_phase": ThermalConductivityWMS,
        ...
      },
      "Liq": {
        "type": LiquidPhase,
        ...
        "therm_cond_phase": NoMethod,
        ...
      },
      ...
    }
    ...
  }
