Heat Exchangers (0D)
====================

.. module:: idaes.unit_models.heat_exchanger

The :code:`idaes.unit_models.heat_exchanger` module contains 0D models for heat
transfer. There are two models available Heater and HeatExchanger.  Heater is a
simple model to add or remove heat from a material.  HeatExchanger models
transfer of heat from one material to another.

The Heater and HeatExchanger models can be imported from :code:`idaes.unit_models`,
while additional rules and utility functions can be imported from
:code:`idaes.unit_models.heat_exchanger`.

.. index::
  pair: idaes.unit_models.heat_exchanger;Heater

Heater
------

The Heater model is a simple 0D model that removes heat from a material stream.
A simple example of how to use a Heater model is given below. A heater model
has a standard inlet and outlet port and a :code:`heat_duty` variable which
specifies a heat transfer rate. The Heater contains a ControlVolume0D, see the
ControlVolume documentation for more information.

The following code provides a simple Heater example.

.. code-block:: python

  import pyomo.environ as pe # Pyomo environment
  from idaes.core import FlowsheetBlock, StateBlockBase
  from idaes.unit_models import Heater, HeatExchanger
  from idaes.unit_models.heat_exchanger import delta_temperature_amtd_rule
  from idaes.property_models import iapws95_ph

  if __name__ == "__main__":
    # Create an empyty flowsheet and steam property parameter block.
    model = pe.ConcreteModel()
    model.fs = FlowsheetBlock(default={"dynamic": False})
    model.fs.properties = iapws95_ph.Iapws95ParameterBlock()

    # Add a Heater model to the flowsheet.
    model.fs.heater = Heater(default={"property_package": model.fs.properties})

    # Setup the heater model by fixing the inputs and heat duty
    model.fs.heater.inlet[:].enth_mol.fix(4000)
    model.fs.heater.inlet[:].flow_mol.fix(100)
    model.fs.heater.inlet[:].pressure.fix(101325)
    model.fs.heater.heat_duty[:].fix(100*20000)

    # Initialize the model.
    model.fs.heater.initialize()

Heater class reference documentation is provided below.

.. autoclass:: Heater
  :members:

.. autoclass:: HeaterData
  :members:

.. index::
  pair: idaes.unit_models.heat_exchanger;HeatExchanger

HeatExchanger
-------------

The HeatExchanger model contains two ControlVolume0D blocks (side_1 and side_2),
which are configured the same as the ControlVolume0D in the Heater model. The
HeatExchanger model contains additional constraints that calculate the amount of
heat transferred from side_1 to side_2. By default the following equation
is used to calculate heat transfer:

.. math::
  Q = UA\Delta T.

Where:

:|Q|:
  Exchanger heat duty or heat transferred from side_1 to side_2,
  **attribute name:** duty, **type:** Var, reference to side_2.heat
:|A|:
  Heat exchange area, the model attribute is "area" with the type Pyomo Var
:|U|:
  Heat transfer coefficient the model attribute is "heat_transfer_coefficient"
  with the type Pyomo Var
:|deltaT|:
  Temperature difference the model attribute is "delta_temperature" with the type
  Pyomo Expression

.. |Q| replace:: :math:`Q`
.. |A| replace:: :math:`A`
.. |U| replace:: :math:`U`
.. |deltaT| replace:: :math:`\Delta T`

By default :math:`\Delta T` is calculated as the log-mean temperature difference
and :math:`U` is fixed.  The equations for these three quantities can be changed
by providing alternate rules to the configuration.

Class Documentation
~~~~~~~~~~~~~~~~~~~

.. autoclass:: HeatExchanger
  :members:

.. autoclass:: HeatExchangerData
  :members:

Rules
~~~~~

A selection of functions for Pyomo rules are provided in the
:code:`heat_exchanger` module, with provide options for different calculation
methods. Users can also provide their own rule functions. See the source for
the rules below for examples.


**Rules for the :code:`delta_temperature_rule` Option**

These rules provide expressions for the temperature difference used in the
heat transfer equations.

.. autofunction:: delta_temperature_lmtd_rule

.. autofunction:: delta_temperature_amtd_rule


**Rules for the :code:`heat_transfer_rule` Option**

These rules provide constraints for the heat transfer rate.

.. autofunction:: heat_transfer_rule


**Rules for the :code:`heat_transfer_coefficient_rule` Option**

There are currently no rules provided for heat transfer coefficient calculation.
When the rule is set to None :code:`heat_transfer_coefficient` is a fixed Var.
User provided heat transfer coefficient rules should return a constraint.
