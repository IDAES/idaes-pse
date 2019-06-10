HeatExchanger (0D)
==================

.. index::
  pair: idaes.unit_models.heat_exchanger;HeatExchanger

.. module:: idaes.unit_models.heat_exchanger

The HeatExchanger model can be imported from :code:`idaes.unit_models`,
while additional rules and utility functions can be imported from
``idaes.unit_models.heat_exchanger``.

Example
-------

The example below demonstrates how to initialize the HeatExchanger model, and
override the default temperature difference calculation.

.. testcode::

  import pyomo.environ as pe # Pyomo environment
  from idaes.core import FlowsheetBlock, StateBlock
  from idaes.unit_models import HeatExchanger
  from idaes.unit_models.heat_exchanger import delta_temperature_amtd_rule
  from idaes.property_models import iapws95_ph

  # Create an empty flowsheet and steam property parameter block.
  model = pe.ConcreteModel()
  model.fs = FlowsheetBlock(default={"dynamic": False})
  model.fs.properties = iapws95_ph.Iapws95ParameterBlock()

  # Add a Heater model to the flowsheet.
  model.fs.heat_exchanger = HeatExchanger(default={
          "delta_temperature_rule":delta_temperature_amtd_rule,
          "side_1":{"property_package": model.fs.properties},
          "side_2":{"property_package": model.fs.properties}})

  model.fs.heat_exchanger.area.fix(1000)
  model.fs.heat_exchanger.overall_heat_transfer_coefficient[0].fix(100)
  model.fs.heat_exchanger.inlet_1.flow_mol.fix(100)
  model.fs.heat_exchanger.inlet_1.pressure.fix(101325)
  model.fs.heat_exchanger.inlet_1.enth_mol.fix(4000)
  model.fs.heat_exchanger.inlet_2.flow_mol.fix(100)
  model.fs.heat_exchanger.inlet_2.pressure.fix(101325)
  model.fs.heat_exchanger.inlet_2.enth_mol.fix(3000)

  # Initialize the model
  model.fs.heat_exchanger.initialize()

Degrees of Freedom
------------------

Aside from the inlet conditions, a heat exchanger model usually has two degrees
of freedom, which can be fixed for it to be fully specified:

* heat transfer area
* heat transfer coefficient.

The user may also provide constants to calculate the heat transfer coefficient.

Model Structure
---------------

The ``HeatExchanger`` model contains two ``ControlVolume0DBlock`` blocks (side_1 and side_2),
which are configured the same as the ``ControlVolume0DBlock`` in the
:ref:`Heater model <models/heater:Heater>`. The ``HeatExchanger`` model contains additional
constraints that calculate the amount of heat transferred from side_1 to side_2.

The ``HeatExchanger`` has two inlet ports inlet_1 (inlet for side_1) and inlet_2
(outlet for side_2), and two outlet ports inlet ports inlet_1 (outlet for side_1)
and outlet_2 (outlet for side_2).

If the :math:`\Delta T` calculation method requires one side to be hotter than
the other, side_1 is assumed to be the hot side.

Variables
---------

=========================== ================== =========== ======================================================================
Variable                    Symbol             Index Sets  Doc
=========================== ================== =========== ======================================================================
heat_duty                   :math:`Q`          t           Heat transferred from side_1 to side_2 a reference to side_2.heat
area                        :math:`A`          None        Heat transfer area
heat_transfer_coefficient   :math:`U`          t           Heat transfer coefficient
delta_temperature           :math:`\Delta T`   t           Temperature difference for heat transfer calculations defaults to LMTD
=========================== ================== =========== ======================================================================


Constraints
-----------

The default constants can be overridden by providing :ref:`alternative rules <models/heat_exchanger:Rules>` for
the heat transfer equation, temperature difference, and heat transfer coefficient. The section
describes the default constraints.

Heat transfer from side_1 to side_2:

.. math::
  Q = UA\Delta T


Temperature difference is an expression:

.. math::
  \Delta T = \frac{\Delta T_1 - \Delta T_2}{\log_e \Delta T_1 - \log_e \Delta T_2}

The heat transfer coefficient is a variable with no associated constraints by default.


See :class:`HeatExchanger` and :class:`HeatExchangerData` for full class details.

.. .. autoclass:: HeatExchanger
..   :members:
..
.. .. autoclass:: HeatExchangerData
..   :members:

Rules
-----

A selection of functions for Pyomo rules are provided in the
``idaes.unit_models.heat_exchanger`` module, which provide options for different
calculation methods. Users can also provide their own rule functions. See the
source code for the rules below for examples.

Rules for the ``delta_temperature_rule`` Option
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These rules provide expressions for the temperature difference used in the
heat transfer equations.

.. autofunction:: delta_temperature_lmtd_rule

.. autofunction:: delta_temperature_amtd_rule


Rules for the ``heat_transfer_rule`` Option
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These rules provide constraints for the heat transfer rate.

.. autofunction:: _heat_transfer_rule


Rules for the ``heat_transfer_coefficient_rule`` Option
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are currently no rules provided for heat transfer coefficient calculation.
When the rule is set to `None`, ``heat_transfer_coefficient`` is a fixed Var.
User provided heat transfer coefficient rules should return a constraint.
