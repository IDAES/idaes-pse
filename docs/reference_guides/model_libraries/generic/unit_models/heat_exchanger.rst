HeatExchanger (0D)
==================

.. index::
   pair: idaes.models.unit_models.heat_exchanger;HeatExchanger

.. currentmodule:: idaes.models.unit_models.heat_exchanger

The HeatExchanger model can be imported from :code:`idaes.models.unit_models`,
while additional rules and utility functions can be imported from
``idaes.models.unit_models.heat_exchanger``.

Example
-------

The example below demonstrates how to initialize the HeatExchanger model, and
override the default temperature difference calculation.

.. testcode::

  import pyomo.environ as pe # Pyomo environment
  from idaes.core import FlowsheetBlock, StateBlock
  from idaes.models.unit_models import HeatExchanger
  from idaes.models.unit_models.heat_exchanger import delta_temperature_amtd_callback
  from idaes.models.properties import iapws95

  # Create an empty flowsheet and steam property parameter block.
  model = pe.ConcreteModel()
  model.fs = FlowsheetBlock(default={"dynamic": False})
  model.fs.properties = iapws95.Iapws95ParameterBlock()

  # Add a Heater model to the flowsheet.
  model.fs.heat_exchanger = HeatExchanger(default={
          "delta_temperature_callback":delta_temperature_amtd_callback,
          "shell":{"property_package": model.fs.properties},
          "tube":{"property_package": model.fs.properties}})

  model.fs.heat_exchanger.area.fix(1000)
  model.fs.heat_exchanger.overall_heat_transfer_coefficient[0].fix(100)
  model.fs.heat_exchanger.shell_inlet.flow_mol.fix(100)
  model.fs.heat_exchanger.shell_inlet.pressure.fix(101325)
  model.fs.heat_exchanger.shell_inlet.enth_mol.fix(4000)
  model.fs.heat_exchanger.tube_inlet.flow_mol.fix(100)
  model.fs.heat_exchanger.tube_inlet.pressure.fix(101325)
  model.fs.heat_exchanger.tube_inlet.enth_mol.fix(3000)

  # Initialize the model
  model.fs.heat_exchanger.initialize()

Degrees of Freedom
------------------

Aside from the inlet conditions, a heat exchanger model usually has two degrees
of freedom, which can be fixed for it to be fully specified. Things that are
frequently fixed are two of:

* heat transfer area,
* heat transfer coefficient, or
* temperature approach.

The user may also provide constraints to calculate the heat transfer coefficient.

Model Structure
---------------

The ``HeatExchanger`` model contains two ``ControlVolume0DBlock`` blocks. By default the
hot side is named ``shell`` and the cold side is named ``tube``. These names are configurable.
The sign convention is that duty is positive for heat flowing from the hot side to the cold
side.  Aside from the sign convention there is no requirement that the hot side be hotter
than the cold side.

The control volumes are configured the same as the ``ControlVolume0DBlock`` in the
:ref:`Heater model <reference_guides/model_libraries/generic/unit_models/heater:Heater>`. The ``HeatExchanger`` model contains additional
constraints that calculate the amount of heat transferred from the hot side to the cold side.

The ``HeatExchanger`` has two inlet ports and two outlet ports. By default these are
``shell_inlet``, ``tube_inlet``, ``shell_outlet``, and ``tube_outlet``. If the user
supplies different hot and cold side names the inlet and outlets are named accordingly.

Variables
---------

=========================== ================== =========== =============================================================================
Variable                    Symbol             Index Sets  Doc
=========================== ================== =========== =============================================================================
heat_duty                   :math:`Q`          t           Heat transferred from hot side to the cold side
area                        :math:`A`          None        Heat transfer area
heat_transfer_coefficient   :math:`U`          t           Heat transfer coefficient
delta_temperature           :math:`\Delta T`   t           Temperature difference, defaults to LMTD
=========================== ================== =========== =============================================================================

Note: ``delta_temperature`` may be either a variable or expression depending on the callback used.  If the specified cold side is hotter
than the specified hot side this value will be negative.

Constraints
-----------

The default constants can be overridden by providing :ref:`alternative rules <reference_guides/model_libraries/generic/unit_models/heat_exchanger:Callbacks>` for
the heat transfer equation, temperature difference, and heat transfer coefficient. The section
describes the default constraints.

Heat transfer from shell to tube:

.. math::
  Q = UA\Delta T


Temperature difference is an expression:

.. math::
  \Delta T = \frac{\Delta T_1 - \Delta T_2}{\log_e\left(\frac{\Delta T_1}{\Delta T_2}\right)}

The heat transfer coefficient is a variable with no associated constraints by default.

Class Documentation
-------------------

.. Note::
  The ``hot_side_config`` and ``cold_side_config`` can also be supplied using the name of
  the hot and cold sides (``shell`` and ``tube`` by default) as in :ref:`the example <reference_guides/model_libraries/generic/unit_models/heat_exchanger:Example>`.

.. autoclass:: HeatExchanger
   :members:

.. autoclass:: HeatExchangerData
   :members:

Callbacks
---------

A selection of functions for constructing the ``delta_temperature`` variable or
expression are provided in the ``idaes.models.unit_models.heat_exchanger`` module.
The user may also provide their own function. These callbacks should all take
one argument (the HeatExchanger block). With the block argument, the function
can add any additional variables, constraints, and expressions needed.  The only
requirement is that either a variable or expression called ``delta_temperature``
must be added to the block.

Defined Callbacks for the ``delta_temperature_callback`` Option
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These callbacks provide expressions for the temperature difference used in the
heat transfer equations. There is a choice of three forms for the LMTD
calculation.  Depending on the expected approach temperatures, one form may be
more favorable.  The first two forms do not require one side or the other to
be the hot side, while the third form requires the hot side to be the hot side
to avoid an evaluation error in the log function.

.. autofunction:: delta_temperature_lmtd_callback

.. autofunction:: delta_temperature_lmtd2_callback

.. autofunction:: delta_temperature_lmtd3_callback

.. autofunction:: delta_temperature_amtd_callback

.. autofunction:: delta_temperature_underwood_callback
