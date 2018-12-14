Pressure Changer
================

The IDAES Pressure Changer model represents a unit operation with a single stream of material which undergoes a change in pressure due to the application of a work. The Pressure Changer model contains support for a number of different thermodynamic assumptions regarding the working fluid.

Degrees of Freedom
------------------

Pressure Changer units generally have one or more degrees of freedom, depending on the thermodynamic assumption used.

Typical fixed variables are:

* outlet pressure, :math:`P_{ratio}` or :math:`\Delta P`,
* unit efficiency (isentropic or pump assumption).

Model Structure
---------------

The core Pressure Changer unit model consists of a single Holdup0D (named holdup) with one Inlet Port (named inlet) and one Outlet Port (named outlet).

Construction Arguments
----------------------

Pressure Changers have the following construction arguments:

* compressor - argument indicates whether the unit should be considered a compressor (True, default) or an expander/turbine (False). This determines how unit efficiency is calculated.
* thermodynamic_assumption - indicates which thermodynamic assumption should be used when constructing the model. Options are:

    - 'isothermal' - (default) assumes no temperature change occurs between the inlet and outlet of the unit.
    - 'adiabatic' - assumes no heat loss occurs between the inlet and outlet of the unit.
    - 'isentropic' - assumes isentropic behavior. This requires an additional set of property calculations for the isentropic outlet conditions.
    - 'pump' - assumes that the fluid work is proportional to the pressure difference and flow rate of fluid. This is suitable for incompressible fluids.

* property_package - property package to use when constructing Property Blocks (default = 'use_parent_value'). This is provided as a Property Parameter Block by the Flowsheet when creating the model. If a value is not provided, the Holdup Block will try to use the default property package if one is defined.
* property_package_args - set of arguments to be passed to the Property Blocks when they are created.
* inlet_list - list of names to be passed to the build_inlets method (default = None).
* num_inlets - number of inlets argument to be passed to the build_inlets method (default = None).
* outlet_list - list of names to be passed to the build_outlets method (default = None).
* num_outlets - number of outlets argument to be passed to the build_outlets method (default = None).

Additionally, Pressure Changers have the following construction arguments which are passed to the Holdup Block for determining which terms to construct in the balance equations.

========================= =================
Argument                  Default Value
========================= =================
material_balance_type     'component_phase'
energy_balance_type       'total'
momentum_balance_type     'total'
dynamic                   False
include_holdup            False
has_rate_reactions        False
has_equilibrium_reactions False
has_phase_equilibrium     True
has_mass_transfer         False
has_heat_transfer         False
has_work_transfer         True
has_pressure_change       True
========================= =================

Additional Constraints
----------------------

In addition to the Constraints written by the Holdup Block, Pressure Changer writes additional Constraints which depend on the thermodynamic assumption chosen. All Pressure Changers add the following Constraint to calculate the pressure ratio:

.. math:: P_{ratio,t} \times P_{in,t} = P_{out,t}

Isothermal Assumption
^^^^^^^^^^^^^^^^^^^^^

The isothermal assumption writes one additional Constraint:

.. math:: T_{out} = T_{in}

Adiabatic Assumption
^^^^^^^^^^^^^^^^^^^^^

The isothermal assumption writes one additional Constraint:

.. math:: Q_{out} = Q_{in}

Isentropic Assumption
^^^^^^^^^^^^^^^^^^^^^

The isentropic assumption creates an additional set of Property Blocks (indexed by time) for the isentropic fluid calculations (named properties_isentropic). This requires a set of balance equations relating the inlet state to the isentropic conditions, which are shown below:

.. math:: F_{in,t,p,j} = F_{out,t,p,j}
.. math:: s_{in,t} = s_{isentropic,t}
.. math:: P_{in,t} \times P_{ratio,t} = P_{isentropic,t}

where :math:`F_{t,p,j}` is the flow of component :math:`j` in phase :math:`p` at time :math:`t` and :math:`s` is the specific entropy of the fluid at time :math:`t`.

Next, the isentropic work is calculated as follows:

.. math:: W_{isentropic,t} = \sum_p{H_{isentropic,t,p}} - \sum_p{H_{in,t,p}}

where :math:`H_{t,p}` is the total energy flow of phase :math:`p` at time :math:`t`. Finally, a constraint which relates the fluid work to the actual mechanical work via an efficiency term :math:`\eta`.

If compressor is True, :math:`W_{isentropic,t} = W_{mechanical,t} \times \eta_t`

If compressor is False, :math:`W_{isentropic,t} \times \eta_t = W_{mechanical,t}`

Pump (Incompressible Fluid) Assumption
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The incompressible fluid assumption writes two additional constraints. Firstly, a Constraint is written which relates fluid work to the pressure change of the fluid.

.. math:: W_{fluid,t} = (P_{out,t}-P_{in,t})\times F_{vol,t}

where :math:`F_{vol,t}` is the total volumetric flowrate of material at time :math:`t` (from the outlet Property Block). Secondly, a constraint which relates the fluid work to the actual mechanical work via an efficiency term :math:`\eta`.

If compressor is True, :math::`W_{fluid,t} = W_{mechanical,t} \times \eta_t`

If compressor is False, :math::`W_{fluid,t} \times \eta_t = W_{mechanical,t}`

Variables
---------

Pressure Changers contain the following Variables (not including those contained within the Holdup Block):

=========================== ===================== ===========================================================================
Variable                    Name                  Notes
=========================== ===================== ===========================================================================
:math:`P_{ratio}`           ratioP
:math:`V_t`                 volume                Only if has_rate_reactions = True, reference to holdup.rate_reaction_extent
:math:`W_{mechanical,t}`    work_mechanical       Reference to holdup.work
:math:`W_{fluid,t}`         work_fluid            Pump assumption only
:math:`\eta_{pump,t}`       efficiency_pump       Pump assumption only
:math:`W_{isentropic,t}`    work_isentropic       Isentropic assumption only
:math:`\eta_{isentropic,t}` efficiency_isentropic Isentropic assumption only
=========================== ===================== ===========================================================================

Isentropic Pressure Changers also have an additional Property Block named `properties_isentropic` (attached to the Unit Model, not the Holdup Block).

PressureChangerData Class
-------------------------
