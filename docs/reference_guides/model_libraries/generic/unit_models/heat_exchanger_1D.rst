Heat Exchangers (1D)
====================

Heat Exchanger models represents a unit operation with two material streams which exchange heat. The IDAES 1-D Heat
Exchanger model is used for detailed modeling of heat exchanger units with variations in one spatial dimension.
For a simpler representation of a heat exchanger unit see
:ref:`Heat Exchanger (0-D) <reference_guides/model_libraries/generic/unit_models/heat_exchanger:HeatExchanger (0D)>`.

Degrees of Freedom
------------------

1-D Heat Exchangers generally have 2 + number of finite elements degrees of freedom.

Typical fixed variables are:

* heat transfer area,
* heat exchanger length,
* average heat transfer coefficients (at all spatial points).

For dynamic simulations (and cases where velocities are required), the cross-sectional areas of the hot and
cold sides need to be provided as well (``unit.hot_side.area`` and ``unit.hot_side.area``).

Model Structure
---------------

The core 1-D Heat Exchanger Model unit model consists of two ControlVolume1DBlock Blocks named ``hot_side``
and ``cold_side``, each with one Inlet Port (named ``hot_side_inlet`` and ``cold_side_inlet``) and one Outlet Port
(named ``hot_side_outlet`` and ``cold_side_outlet``). These names are configurable using the ``hot_side_name`` and
``cold_side_name`` configuration arguments, in which case aliases are assigned to the control volumes and associated
Ports using the names provided (note that ``hot_side`` and ``cold_side`` will always work).

Construction Arguments
----------------------

1-D Heat Exchanger units have construction arguments specific to the hot and cold sides and for the unit as a
whole.

Arguments that are applicable to the heat exchanger unit are as follows:

* flow_type - indicates the flow arrangement within the unit to be modeled. Options are:

    - 'co-current' - (default) shell and tube both flow in the same direction (from x=0 to x=1)
    - 'counter-current' - shell and tube flow in opposite directions (shell from x=0 to x=1 and tube from x=1 to x=0).

* finite_elements - sets the number of finite elements to use when discretizing the spatial domains (default = 20). This is used for both shell and tube side domains.
* collocation_points - sets the number of collocation points to use when discretizing the spatial domains (default = 5, collocation methods only). This is used for both shell and tube side domains.
* hot_side_name
* cold_side_name

Arguments that are applicable to the hot and cold sides:

* property_package - property package to use when constructing shell side Property Blocks (default = 'use_parent_value'). This is provided as a Physical Parameter Block by the Flowsheet when creating the model. If a value is not provided, the ControlVolume Block will try to use the default property package if one is defined.
* property_package_args - set of arguments to be passed to the shell side Property Blocks when they are created.
* transformation_method - argument to specify the DAE transformation method for the shell side; should be compatible with the Pyomo DAE TransformationFactory
* transformation_scheme - argument to specify the scheme to use for the selected DAE transformation method; should be compatible with the Pyomo DAE TransformationFactory

Additionally, 1-D Heat Exchanger units have the following construction arguments for each side which are passed
to the ControlVolume1DBlocks for determining which terms to construct in the balance equations for the
hot and cold sides.

========================= =================
Argument                  Default Value
========================= =================
dynamic                   useDefault
has_holdup                False
material_balance_type     'componentTotal'
energy_balance_type       'enthalpyTotal'
momentum_balance_type     'pressureTotal'
has_phase_equilibrium     False
has_heat_transfer         True
has_pressure_change       False
========================= =================

Variables
---------

1-D Heat Exchanger units add the following additional Variables beyond those created by the ControlVolume1DBlock Block.

================= ========================== ==================================
Variable          Name                       Notes
================= ========================== ==================================
:math:`L`         length                     Reference to hot_side.length
:math:`A`         area                       Overall heat transfer area
:math:`U_{t,x}`   heat_transfer_coefficient  Average heat transfer coefficient
================= ========================== ==================================

Additional Constraints
----------------------

1-D Heat Exchanger models write the following additional Constraints to describe the heat transfer between the
two sides of the heat exchanger. Firstly, overall heat transfer is calculated as:

.. math:: Q_{hot,t,x} = -U_{t,x} \times \frac{A}{L_{hot}} \times (T_{hot,t,x}-T_{cold,t,x}))

where :math:`Q_{hot,t,x}` is the hot-side heat duty at point :math:`x` and time :math:`t`, :math:`A` is the total heat
transfer area, :math:`U_{t,x}` is the average heat transfer coefficient, and :math:`T_{hot,t,x}` and
:math:`T_{cold,t,x}` are the hot and cold side temperatures respectively.

Next, overall heat conservation is enforced by the following constraint:

.. math:: Q_{cold,t,x} = -Q_{hot,t,x}

Finally, the following Constraints are written to describe the unit geometry:

.. math:: L_{hot} = L_{cold}

where :math:`L_{hot}` and :math:`L_{cold}` are the length of the hot and cold side respectively.

Initialization
--------------
.. module:: idaes.models.unit_models.heat_exchanger_1D

.. autoclass:: HX1DInitializer
   :members: initialization_routine

HeatExchanger1d Class
---------------------

.. autoclass:: HeatExchanger1D
  :members:

HeatExchanger1dData Class
-------------------------

.. autoclass:: HeatExchanger1DData
    :members:
