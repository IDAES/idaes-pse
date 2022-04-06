Heat Exchangers (1D)
====================================

Heat Exchanger models represents a unit operation with two material streams which exchange heat. The IDAES 1-D Heat Exchanger model is used for detailed modeling of heat exchanger units with variations in one spatial dimension. For a simpler representation of a heat exchanger unit see Heat Exchanger (0-D).

Degrees of Freedom
------------------

1-D Heat Exchangers generally have 7 degrees of freedom.

Typical fixed variables are:

* shell length and diameter,
* tube length and diameter,
* number of tubes,
* heat transfer coefficients (at all spatial points) for both shell and tube sides.

Model Structure
---------------

The core 1-D Heat Exchanger Model unit model consists of two ControlVolume1DBlock Blocks (named shell and tube), each with one Inlet Port (named shell_inlet and tube_inlet) and one Outlet Port (named shell_outlet and tube_outlet).

Construction Arguments
----------------------

1-D Heat Exchanger units have construction arguments specific to the shell side, tube side and for the unit as a whole.

Arguments that are applicable to the heat exchanger unit are as follows:

* flow_type - indicates the flow arrangement within the unit to be modeled. Options are:

    - 'co-current' - (default) shell and tube both flow in the same direction (from x=0 to x=1)
    - 'counter-current' - shell and tube flow in opposite directions (shell from x=0 to x=1 and tube from x=1 to x=0).

* finite_elements - sets the number of finite elements to use when discretizing the spatial domains (default = 20). This is used for both shell and tube side domains.
* collocation_points - sets the number of collocation points to use when discretizing the spatial domains (default = 5, collocation methods only). This is used for both shell and tube side domains.
* has_wall_conduction - option to enable a model for heat conduction across the tube wall:
    - 'none' - 0D wall model
    - '1D' - 1D heat conduction equation along the thickness of the tube wall
    - '2D' - 2D heat conduction equation along the length and thickness of the tube wall

Arguments that are applicable to the shell side:

* property_package - property package to use when constructing shell side Property Blocks (default = 'use_parent_value'). This is provided as a Physical Parameter Block by the Flowsheet when creating the model. If a value is not provided, the ControlVolume Block will try to use the default property package if one is defined.
* property_package_args - set of arguments to be passed to the shell side Property Blocks when they are created.
* transformation_method - argument to specify the DAE transformation method for the shell side; should be compatible with the Pyomo DAE TransformationFactory
* transformation_scheme - argument to specify the scheme to use for the selected DAE transformation method; should be compatible with the Pyomo DAE TransformationFactory

Arguments that are applicable to the tube side:

* property_package - property package to use when constructing tube side Property Blocks (default = 'use_parent_value'). This is provided as a Property Parameter Block by the Flowsheet when creating the model. If a value is not provided, the ControlVolume Block will try to use the default property package if one is defined.
* property_package_args - set of arguments to be passed to the tube side Property Blocks when they are created.
* transformation_method - argument to specify the DAE transformation method for the tube side; should be compatible with the Pyomo DAE TransformationFactory
* transformation_scheme - argument to specify the scheme to use for the selected DAE transformation method; should be compatible with the Pyomo DAE TransformationFactory


Additionally, 1-D Heat Exchanger units have the following construction arguments which are passed to the ControlVolume1DBlock Block for determining which terms to construct in the balance equations for the shell and tube side.

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

Additional Constraints
----------------------

1-D Heat Exchanger models write the following additional Constraints to describe the heat transfer between the two sides of the heat exchanger. Firstly, the shell- and tube-side heat transfer is calculated as:

.. math :: Q_{shell,t,x} = - N_{tubes} \times (\pi \times U_{shell,t,x} \times D_{tube,outer} \times (T_{shell,t,x}-T_{wall,t,x}))

where :math:`Q_{shell,t,x}` is the shell-side heat duty at point :math:`x` and time :math:`t`, :math:`N_{tubes}` :math:`D_{tube}` are the number of and diameter of the tubes in the heat exchanger, :math:`U_{shell,t,x}` is the shell-side heat transfer coefficient, and :math:`T_{shell,t,x}` and :math:`T_{wall,t,x}` are the shell-side and tube wall temperatures respectively.

.. math :: Q_{tube,t,x} = N_{tubes} \times (\pi \times U_{tube,t,x} \times D_{tube,inner} \times (T_{wall,t,x}-T_{tube,t,x}))

where :math:`Q_{tube,t,x}` is the tube-side heat duty at point :math:`x` and time :math:`t`, :math:`U_{tube,t,x}` is the tube-side heat transfer coefficient and :math:`T_{tube,t,x}` is the tube-side temperature.

If a OD wall model is used for the tube wall conduction, the following constraint is implemented to connect the heat terms on the shell and tube side:

.. math :: N_{tubes} \times Q_{tube,t,x} =  - Q_{shell,t,x}

Finally, the following Constraints are written to describe the unit geometry:

.. math:: 4 \times A_{tube} = \pi \times D_{tube}^2

.. math:: 4 \times A_{shell} = \pi \times (D_{shell}^2 - N_{tubes} \times D_{tube}^2)

where :math:`A_{shell}` and :math:`A_{tube}` are the shell and tube areas respectively and :math:`D_{shell}` and :math:`D_{tube}` are the shell and tube diameters.

Variables
---------

1-D Heat Exchanger units add the following additional Variables beyond those created by the ControlVolume1DBlock Block.

====================== =============================== =========================
Variable               Name                            Notes
====================== =============================== =========================
:math:`L_{shell}`      shell_length                    Reference to shell.length
:math:`A_{shell}`      shell_area                      Reference to shell.area
:math:`D_{shell}`      d_shell
:math:`L_{tube}`       tube_length                     Reference to tube.length
:math:`A_{tube}`       tube_area                       Reference to tube.area
:math:`D_{tube}`       d_tube
:math:`N_{tubes}`      N_tubes
:math:`T_{wall,t,x}`   temperature_wall
:math:`U_{shell,t,x}`  shell_heat_transfer_coefficient
:math:`U_{tube,t,x}`   tube_heat_transfer_coefficient
====================== =============================== =========================

HeatExchanger1dClass
--------------------

.. module:: idaes.models.unit_models.heat_exchanger_1D

.. autoclass:: HeatExchanger1D
  :members:

HeatExchanger1dDataClass
------------------------

.. autoclass:: HeatExchanger1DData
    :members:
