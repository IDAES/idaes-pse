1D Shell and Tube Heat Exchanger
================================

The IDAES Shell and Tube Heat Exchanger model is an extension of the standard
:ref:`1-D Heat Exchanger model <reference_guides/model_libraries/generic/unit_models/heat_exchanger_1D:Heat Exchangers (1D)>`
model which includes the effects of wall temperature and geometry of co- and counter-current flow shell and tube
type heat exchangers.

Degrees of Freedom
------------------

1-D Shell and Tube Heat Exchangers generally have 5 + 2 times number of finite elements degrees of freedom.

Typical fixed variables are:

* heat exchanger length,
* shell diameter,
* tube inner and outer diameters,
* number of tubes,
* hot and cold side heat transfer coefficients (at all spatial points).

Model Structure
---------------

The core 1-D Heat Exchanger Model unit model consists of two ControlVolume1DBlock Blocks named ``hot_side``
and ``cold_side``, each with one Inlet Port (named ``hot_side_inlet`` and ``cold_side_inlet``) and one Outlet Port
(named ``hot_side_outlet`` and ``cold_side_outlet``). These names are configurable using the ``hot_side_name`` and
``cold_side_name`` configuration arguments, in which case aliases are assigned to the control volumes and associated
Ports using the names provided (note that ``hot_side`` and ``cold_side`` will always work). If custom names are not
provided, then default names of ``shell`` and ``tube`` will be used based on the ``shell_is_hot`` configuration
argument.

Construction Arguments
----------------------

Shell and Tube Heat Exchanger units have construction arguments specific to the hot and cold sides and for the unit as a
whole.

Arguments that are applicable to the heat exchanger unit are as follows:

* flow_type - indicates the flow arrangement within the unit to be modeled. Options are:

    - 'co-current' - (default) shell and tube both flow in the same direction (from x=0 to x=1)
    - 'counter-current' - shell and tube flow in opposite directions (shell from x=0 to x=1 and tube from x=1 to x=0).

* finite_elements - sets the number of finite elements to use when discretizing the spatial domains (default = 20). This is used for both shell and tube side domains.
* collocation_points - sets the number of collocation points to use when discretizing the spatial domains (default = 5, collocation methods only). This is used for both shell and tube side domains.
* hot_side_name
* cold_side_name
* shell_is_hot - (``bool``) indicates whether the shell will be mapped to the hot side of the heat exchanger or the cold side.

Arguments that are applicable to the hot and cold sides:

* property_package - property package to use when constructing shell side Property Blocks (default = 'use_parent_value'). This is provided as a Physical Parameter Block by the Flowsheet when creating the model. If a value is not provided, the ControlVolume Block will try to use the default property package if one is defined.
* property_package_args - set of arguments to be passed to the shell side Property Blocks when they are created.
* transformation_method - argument to specify the DAE transformation method for the shell side; should be compatible with the Pyomo DAE TransformationFactory
* transformation_scheme - argument to specify the scheme to use for the selected DAE transformation method; should be compatible with the Pyomo DAE TransformationFactory

Additionally, Shell and Tube Heat Exchanger units have the following construction arguments for each side which are passed
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

1-D Shell and Tube Heat Exchanger units add the following additional Variables beyond those created by the
ControlVolume1DBlock Block.

======================== ==================================== ==================================
Variable                 Name                                 Notes
======================== ==================================== ==================================
:math:`L`                length                               Reference to hot_side.length
:math:`D_{shell}`        shell_diameter                       Diameter of shell
:math:`D_{tube, inner}`  tube_inner_diameter                  Inner diameter of tubes
:math:`D_{tube, outer}`  tube_outer_diameter                  Outer diameter of tubes
:math:`N_{tubes}`        number_of_tubes
:math:`T_{wall}`         temperature_wall                     Temperature of tube wall material
:math:`U_{hot,t,x}`      hot_side_heat_transfer_coefficient
:math:`U_{cold,t,x}`     cold_side_heat_transfer_coefficient
======================== ==================================== ==================================

Additional Constraints
----------------------

1-D Shell and TubeHeat Exchanger models write the following additional Constraints to describe the heat transfer between the
two sides of the heat exchanger.

Firstly, heat transfer for the hot and cold sides to the tube wall is calculated as:

.. math:: Q_{hot,t,x} = -U_{hot_t,x} \times N_{tubes} \times \pi \times D_{tube, outer} \times (T_{hot,t,x}-T_{wall,t,x}))

.. math:: Q_{cold,t,x} = U_{cold_t,x} \times N_{tubes} \times \pi \times D_{tube, inner} \times (T_{wall,t,x}-T_{cold,t,x}))

where :math:`Q_{hot,t,x}` and :math:`Q_{cold,t,x}` are the hot -and cold-side heat duties at point :math:`x` and time
:math:`t`.

Next, overall heat conservation is enforced by the following constraint:

.. math:: Q_{cold,t,x} = -Q_{hot,t,x}

Finally, the following Constraints are written to describe the unit geometry:

.. math:: L_{hot} = L_{cold}

where :math:`L_{hot}` and :math:`L_{cold}` are the length of the hot and cold side respectively.

.. math:: A_{shell} = \pi \times \frac{(D_{shell}^2 - N_{tubes} \times D_{tube, outer}^2)}{4}

.. math:: A_{tube} = N_{tubes} \times \pi \times \frac{D_{tube,inner}^2}{4}

where :math:`A_{shell}` and :math:`A_{tube}` are the cross-sectional areas of the shell and tube side control volumes
respectively.

Initialization
--------------
.. module:: idaes.models.unit_models.shell_and_tube_1d

.. autoclass:: ShellAndTubeInitializer
   :members: initialization_routine


ShellAndTube1D Class
--------------------



.. autoclass:: ShellAndTube1D
  :members:

ShellAndTube1Data Class
-----------------------

.. autoclass:: ShellAndTube1DData
    :members:
