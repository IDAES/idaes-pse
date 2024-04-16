CrossFlowHeatExchanger1D
========================

.. index::
  pair: idaes.models_extra.power_generation.unit_models.cross_flow_heat_exchanger_1D;CrossFlowHeatExchanger1D

.. module:: idaes.models_extra.power_generation.unit_models.cross_flow_heat_exchanger_1D

This model is for a cross flow heat exchanger between two gases. 

Example
-------

.. code-block:: python

  import pyomo.environ as pyo

  from idaes.core import FlowsheetBlock
  import idaes.core.util.scaling as iscale
  from idaes.models.unit_models import HeatExchangerFlowPattern
  from idaes.models.properties.modular_properties import GenericParameterBlock
  from idaes.models_extra.power_generation.properties.natural_gas_PR import (
      get_prop,
      EosType,
  )
  from idaes.models_extra.power_generation.unit_models import CrossFlowHeatExchanger1D
  import idaes.core.util.model_statistics as mstat
  from idaes.core.util.model_statistics import degrees_of_freedom
  from idaes.core.solvers import get_solver


  # Set up solver
  optarg = {
      "constr_viol_tol": 1e-8,
      "nlp_scaling_method": "user-scaling",
      "linear_solver": "ma57",
      "OF_ma57_automatic_scaling": "yes",
      "max_iter": 350,
      "tol": 1e-8,
      "halt_on_ampl_error": "no",
  }
  solver = get_solver("ipopt", options=optarg)

  m = pyo.ConcreteModel()
  m.fs = pyo.FlowsheetBlock(dynamic=False)
  m.fs.properties = GenericParameterBlock(
        **get_prop(["H2", "H2O", "Ar", "N2"], {"Vap"}, eos=EosType.IDEAL),
        doc="H2O + H2 gas property parameters",
    )
  m.fs.heat_exchanger = CrossFlowHeatExchanger1D(
        has_holdup=True,
        dynamic=False,
        cold_side={
            "property_package": m.fs.h2_side_prop_params,
            "has_holdup": False,
            "dynamic": False,
            "has_pressure_change": pressure_drop,
            "transformation_method": "dae.finite_difference",
            "transformation_scheme": "FORWARD",
        },
        hot_side={
            "property_package": m.fs.h2_side_prop_params,
            "has_holdup": False,
            "dynamic": False,
            "has_pressure_change": pressure_drop,
            "transformation_method": "dae.finite_difference",
            "transformation_scheme": "BACKWARD",
        },
        shell_is_hot=True,
        flow_type=HeatExchangerFlowPattern.countercurrent,
        finite_elements=12,
        tube_arrangement="staggered",
    )

  hx = m.fs.heat_exchanger

  hx.hot_side_inlet.flow_mol.fix(2619.7)
  hx.hot_side_inlet.temperature.fix(971.6)
  hx.hot_side_inlet.pressure.fix(1.2e5)
  hx.hot_side_inlet.mole_frac_comp[0, "H2"].fix(0.79715)
  hx.hot_side_inlet.mole_frac_comp[0, "H2O"].fix(0.20177)
  hx.hot_side_inlet.mole_frac_comp[0, "Ar"].fix(0.00086358)
  hx.hot_side_inlet.mole_frac_comp[0, "N2"].fix(0.00021589)

  hx.cold_side_inlet.flow_mol.fix(2619.7)
  hx.cold_side_inlet.temperature.fix(446.21)
  hx.cold_side_inlet.pressure.fix(1.2e5)
  hx.cold_side_inlet.mole_frac_comp[0, "H2"].fix(0.36203)
  hx.cold_side_inlet.mole_frac_comp[0, "H2O"].fix(0.63689)
  hx.cold_side_inlet.mole_frac_comp[0, "Ar"].fix(0.00086358)
  hx.cold_side_inlet.mole_frac_comp[0, "N2"].fix(0.00021589)

  hx.di_tube.fix(0.0525018)
  hx.thickness_tube.fix(0.0039116)
  hx.length_tube_seg.fix(4.3)
  hx.nseg_tube.fix(12)
  hx.ncol_tube.fix(50)
  hx.nrow_inlet.fix(25)

  hx.pitch_x.fix(0.1)
  hx.pitch_y.fix(0.1)
  hx.therm_cond_wall = 43.0
  hx.rfouling_tube = 0.0001
  hx.rfouling_shell = 0.0001
  hx.fcorrection_htc_tube.fix(1)
  hx.fcorrection_htc_shell.fix(1)
  if pressure_drop:
      hx.fcorrection_dp_tube.fix(1)
      hx.fcorrection_dp_shell.fix(1)

  hx.cp_wall.value = 502.4

  pp = m.fs.properties
  pp.set_default_scaling("enth_mol_phase", 1e-3)
  pp.set_default_scaling("pressure", 1e-5)
  pp.set_default_scaling("temperature", 1e-2)
  pp.set_default_scaling("flow_mol", 1e-3)

  _mf_scale = {
      "H2": 1,
      "H2O": 1,
      "N2": 10,
      "Ar": 10,
  }
  for comp, s in _mf_scale.items():
      pp.set_default_scaling("mole_frac_comp", s, index=comp)
      pp.set_default_scaling("mole_frac_phase_comp", s, index=("Vap", comp))
      pp.set_default_scaling("flow_mol_phase_comp", s * 1e-3, index=("Vap", comp))

  shell = hx.hot_side
  tube = hx.cold_side
  iscale.set_scaling_factor(shell.area, 1e-1)
  iscale.set_scaling_factor(shell.heat, 1e-6)
  iscale.set_scaling_factor(tube.area, 1)
  iscale.set_scaling_factor(tube.heat, 1e-6)
  iscale.set_scaling_factor(shell._enthalpy_flow, 1e-8)  # pylint: disable=W0212
  iscale.set_scaling_factor(tube._enthalpy_flow, 1e-8)  # pylint: disable=W0212
  iscale.set_scaling_factor(shell.enthalpy_flow_dx, 1e-7)
  iscale.set_scaling_factor(tube.enthalpy_flow_dx, 1e-7)
  iscale.set_scaling_factor(hx.heat_holdup, 1e-8)

  iscale.calculate_scaling_factors(m)

  initializer = CrossFlowHeatExchanger1DInitializer(
        solver="ipopt",
        solver_options=optarg
  )
  initializer.initialize(m.fs.heat_exchanger)

  
  

Heat Exchanger Geometry
-----------------------
=========================== =========== =================================================================================
Variable                    Index Sets  Doc
=========================== =========== =================================================================================
``number_columns_per_pass`` None        Number of columns of tube per pass
``number_rows_per_pass``    None        Number of rows of tube per pass
``number_passes``           None        Number of tube banks of `nrow_tube * ncol_inlet` tubes
``pitch_x``                 None        Distance between columns (TODO rows?) of tubes, measured from center-of-tube to center-of-tube
``pitch_y``                 None        Distance between rows (TODO columns?) of tubes, measured from center-of-tube to center-of-tube
``length_tube_seg``         None        Length of tube segment perpendicular to flow in each pass
``area_flow_shell_min``     None        Minimum flow area on shell side
``di_tube``                 None        Inner diameter of tubes
``thickness_tube``          None        Thickness of tube wall.
=========================== =========== =================================================================================

============================ =========== ===========================================================================
Expression                   Index Sets  Doc
============================ =========== ===========================================================================
``nrow_tube``                None        Total number of rows of tube
``do_tube``                  None        Outer diameter of tube (equal to `di_tube+2*thickness_tube`)
``pitch_x_to_do``            None        Ratio of `pitch_x` to `do_tube`
``pitch_y_to_do``            None        Ratio of `pitch_y` to `do_tube`
``area_wall_seg``            None        Total cross-sectional area of tube per pass
``total_heat_transfer_area`` None        Total heat transfer area, as measured on outer surface of tubes
============================ =========== ===========================================================================

=========================== =========== =================================================================================================
Constraint                  Index Sets  Doc
=========================== =========== =================================================================================================
``length_flow_shell_eqn``   None        Constrains shell flow length from control volume to equal value implied by geometry
``area_flow_shell_eqn``     None        Constrains shell flow cross-sectional area from control volume to equal value implied by geometry
``area_flow_shell_min_eqn`` None        Constraints `area_flow_shell_min` to equal value determined by geometry
``length_flow_tube_eqn``    None        Constrains tube flow length from control volume to equal value implied by geometry
``area_flow_tube_eqn``      None        Constrains tube flow cross-sectional area from control volume to equal value implied by geometry
=========================== =========== =================================================================================================

Performance Equations
-----------------------

================================== ============ =================================================================================
Variable                           Index Sets   Doc
================================== ============ =================================================================================
``fcorrection_htc_shell``          time, length Correction factor for shell convective heat transfer
``conv_heat_transfer_coeff_shell`` time, length Shell-side convective heat transfer coefficient
``temp_wall_shell``                time, length Shell-side wall temperature of tube
``temp_wall_center``               time, length Temperature at center of tube wall
``v_shell``                        time, length Flow velocity on shell side through minimum area
``N_Re_shell``                     time, length Reynolds number on shell side
``N_Nu_shell``                     time, length Nusselt number on shell side
``heat_transfer_coeff_tube``       time, length Tube-side heat transfer coefficient
``temp_wall_tube``                 time, length Tube-side wall temperature of tube
``v_tube``                         time, length Flow velocity of gas in tube
``N_Re_tube``                      time, length Reynolds number in tube
``N_Nu_tube``                      time, length Nusselt number on tube side
================================== ============ =================================================================================

=========================== =========== =================================================================================
Parameter                   Index Sets  Doc
=========================== =========== =================================================================================
``therm_cond_wall``         None        Thermal conductivity of tube wall
``density_wall``            None        Mass density of tube wall metal
``cp_wall``                 None        Tube wall heat capacity (mass basis)
``rfouling_shell``          None        Fouling resistance on shell side
``f_arrangement``           None        Adjustment factor depending on `tube_arrangement` in config
=========================== =========== =================================================================================

====================================== ============ =================================================================================
Constraint                             Index Sets   Doc
====================================== ============ =================================================================================
``v_shell_eqn``                        time, length Calculates velocity of flow through shell using `area_flow_shell_min`
``N_Re_shell_eqn``                     time, length Calculates the shell-side Reynolds number
``conv_heat_transfer_coeff_shell_eqn`` time, length Calculates the shell-side convective heat transfer coefficient 
``v_tube_eqn``                         time, length Calculates gas velocity in tube
``N_Re_tube_eqn``                      time, length Calculates the tube-side Reynolds number
``heat_transfer_coeff_tube_eqn``       time, length Calcualtes the tube-side heat transfer coefficient
====================================== ============ =================================================================================

====================================== ============ ===================================================================================
Expression                             Index Sets   Doc
====================================== ============ ===================================================================================
``total_heat_transfer_coeff_shell``    time, length Returns ``conv_heat_transfer_coeff_shell``. Could be extended to include radiation.
====================================== ============ ===================================================================================

Pressure Change Equations
-------------------------

=========================== ============ =================================================================================
Parameter                    Index Sets   Doc
=========================== ============ =================================================================================
``fcorrection_dp_shell``    None         Correction factor for shell side pressure drop
``kloss_uturn``             None         Loss coefficient of a tube u-turn
``fcorrection_dp_tube``     None         Correction factor for tube side pressure drop
=========================== ============ =================================================================================

=========================== ============ =================================================================================
Variable                    Index Sets   Doc
=========================== ============ =================================================================================
``fcorrection_dp_shell``    None         Correction factor for shell side pressure drop
``friction_factor_shell``   time, length Friction factor on shell side
``friction_factor_tube``    time, length Friction factor on tube side
``deltaP_tube_friction``    time, length Change of pressure in tube due to friction
``deltaP_tube_uturn``       time, length Change of pressure in tube due to U turns
=========================== ============ =================================================================================

================================== ============ =================================================================================
Constraint                         Index Sets   Doc
================================== ============ =================================================================================
``friction_factor_shell_eqn``      time, length Calculates the shell-side friction factor
``deltaP_shell_eqn``               time, length Sets `deltaP_shell` based on the friction factor and shell properties
``friction_factor_tube_eqn``       time, length Calculates the tube-side friction factor
``deltaP_tube_friction_eqn``       time, length Sets `deltaP_tube_friction` based on friction factor
``deltaP_tube_uturn_eqn``          time, length Sets `deltaP_tube_uturn` based on `kloss_uturn`
``deltaP_tube_eqn``                time, length Sets `deltaP_tube` by summing `deltaP_tube_friction` and `deltaP_tube_uturn`
================================== ============ =================================================================================


Holdup Equations
----------------
Created when `has_holdup=True` in the config.
=========================== ============ =================================================================================
Variable                    Index Sets   Doc
=========================== ============ =================================================================================
``heat_holdup``             time, length Energy holdup per unit length of shell flow path
=========================== ============ =================================================================================

=========================== ============ =================================================================================
Constraint                  Index Sets   Doc
=========================== ============ =================================================================================
``heat_holdup_eqn``         time, length Defines heat holdup in terms of geometry and physical properties
=========================== ============ =================================================================================

Dynamic Equations
-----------------
Created when `dynamic=True` in the config.
=========================== ============ =================================================================================
Derivative Variable         Index Sets   Doc
=========================== ============ =================================================================================
``heat_accumulation``       time, length Energy accumulation in tube wall per unit length of shell flow path per unit time
=========================== ============ =================================================================================


Constraints
-----------

The pressure flow relation is added to the inherited constraints from the :ref:`PressureChanger model
<reference_guides/model_libraries/generic/unit_models/pressure_changer:Pressure Changer>`.

If the ``phase`` option is set to ``"Liq"`` the following equation describes the pressure-flow relation.

.. math::

  \frac{1}{s_f^2}F^2 = \frac{1}{s_f^2}C_v^2\left(P_{in} - P_{out}\right)f(x)^2

If the ``phase`` option is set to ``"Vap"`` the following equation describes the pressure-flow relation.

.. math::

  \frac{1}{s_f^2}F^2 = \frac{1}{s_f^2}C_v^2\left(P_{in}^2 - P_{out}^2\right)f(x)^2


Initialization
--------------
First, the shell and tube control volumes are initialized without heat transfer. Next
the total possible heat transfer between streams is estimated based on heat capacity, 
flow rate, and inlet/outlet temperatures. The actual temperature change is set to be 
half the theoretical maximum, and the shell and tube are initalized with linear 
temperature profiles. Finally, temperatures besides the inlets are unfixed and 
the performance equations are activated before a full solve of the system model.


HelmValve Class
----------------

.. autoclass:: CrossFlowHeatExchanger1D
  :members:

HelmValveData Class
---------------------

.. autoclass:: CrossFlowHeatExchanger1DData
  :members:
