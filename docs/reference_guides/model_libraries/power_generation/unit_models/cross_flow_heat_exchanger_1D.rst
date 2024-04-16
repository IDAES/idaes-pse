========================
CrossFlowHeatExchanger1D
========================

.. index::
  pair: idaes.models_extra.power_generation.unit_models.cross_flow_heat_exchanger_1D;CrossFlowHeatExchanger1D

.. module:: idaes.models_extra.power_generation.unit_models.cross_flow_heat_exchanger_1D

This model is for a cross flow heat exchanger between two gases. 

Example
=======

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
=======================

Variables
---------

This model features many variables describing the geometry of the heat exchanger in 
order to accurately calculate heat transfer coefficients, pressure drop, and heat holdup.

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


Expressions
-----------

=========================== =========== ===========================================================================
Expression                  Index Sets  Doc
=========================== =========== ===========================================================================
``nrow_tube``               None        Total number of rows of tube
``do_tube``                 None        Outer diameter of tube (equal to `di_tube+2*thickness_tube`)
``pitch_x_to_do``           None        Ratio of `pitch_x` to `do_tube`
``pitch_y_to_do``           None        Ratio of `pitch_y` to `do_tube`
``area_wall_seg``           None        Total cross-sectional area of tube per pass
=========================== =========== ===========================================================================

Constraints
-----------

=========================== =========== =================================================================================================
Constraint                  Index Sets  Doc
=========================== =========== =================================================================================================
``length_flow_shell_eqn``   None        Constrains shell flow length from control volume to equal value implied by geometry
``area_flow_shell_eqn``     None        Constrains shell flow cross-sectional area from control volume to equal value implied by geometry

=========================== =========== =================================================================================================


===============

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

This just calls the initialization routine from PressureChanger, but it is wrapped in
a function to ensure the state after initialization is the same as before initialization.
The arguments to the initialization method are the same as PressureChanger.

HelmValve Class
----------------

.. autoclass:: HelmValve
  :members:

HelmValveData Class
---------------------

.. autoclass:: HelmValveData
  :members:
