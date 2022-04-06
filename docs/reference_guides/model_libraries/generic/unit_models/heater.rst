Heater
======

.. index::
  pair: idaes.models.unit_models.heater;Heater

.. module:: idaes.models.unit_models.heater

The Heater model is a simple 0D model that adds or removes heat from a
material stream.

Example
-------

.. testcode::

  import pyomo.environ as pe # Pyomo environment
  from idaes.core import FlowsheetBlock, StateBlock
  from idaes.models.unit_models import Heater
  from idaes.models.properties import iapws95

  # Create an empty flowsheet and steam property parameter block.
  model = pe.ConcreteModel()
  model.fs = FlowsheetBlock(default={"dynamic": False})
  model.fs.properties = iapws95.Iapws95ParameterBlock()

  # Add a Heater model to the flowsheet.
  model.fs.heater = Heater(default={"property_package": model.fs.properties})

  # Setup the heater model by fixing the inputs and heat duty
  model.fs.heater.inlet[:].enth_mol.fix(4000)
  model.fs.heater.inlet[:].flow_mol.fix(100)
  model.fs.heater.inlet[:].pressure.fix(101325)
  model.fs.heater.heat_duty[:].fix(100*20000)

  # Initialize the model.
  model.fs.heater.initialize()

Degrees of Freedom
------------------

Aside from the inlet conditions, a heater model usually has one degree of
freedom, which is the heat duty.

Model Structure
---------------

A heater model contains one ControlVolume0DBlock block.

Variables
---------

The ``heat_duty`` variable is a reference to ``control_volume.heat``.

Constraints
-----------

A heater model contains no additional constraints beyond what are contained in
a ``ControlVolume0DBlock`` model.

Heater Class
------------

.. autoclass:: Heater
  :members:

HeaterData Class
----------------

.. autoclass:: HeaterData
  :members:
