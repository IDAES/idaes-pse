Heater1D
========

.. index::
  pair: idaes.models_extra.power_generation.unit_models.heater_1D;Heater1D

.. module:: idaes.models_extra.power_generation.unit_models.heater_1D

This model is for a gas trim heater modeled as gas being blown perpendicularly across banks of hollow tubes,
which are heated by resistive heating. Note that the ``finite_elements`` option in the control
volume config should be set to an integer factor of ``number_passes`` in order for the
discretization equations to make sense as a cross-flow heat exchanger.

Example
-------

.. code-block:: python

    import pyomo.environ as pyo

    from idaes.core import FlowsheetBlock
    import idaes.core.util.scaling as iscale
    from idaes.models.properties.modular_properties import GenericParameterBlock
    from idaes.models_extra.power_generation.properties.natural_gas_PR import (
        get_prop,
        EosType,
    )
    from idaes.models_extra.power_generation.unit_models import Heater1D
    from idaes.core.util.model_statistics import degrees_of_freedom

    optarg = {
        "constr_viol_tol": 1e-8,
        "nlp_scaling_method": "user-scaling",
        "linear_solver": "ma57",
        "OF_ma57_automatic_scaling": "yes",
        "max_iter": 350,
        "tol": 1e-8,
        "halt_on_ampl_error": "no",
    }

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.h2_side_prop_params = GenericParameterBlock(
        **get_prop(["H2", "H2O", "Ar", "N2"], {"Vap"}, eos=EosType.IDEAL),
        doc="H2O + H2 gas property parameters",
    )
    m.fs.heater = Heater1D(
        property_package=m.fs.h2_side_prop_params,
        has_holdup=True,
        dynamic=False,
        has_fluid_holdup=False,
        has_pressure_change=pressure_drop,
        finite_elements=4,
        tube_arrangement="in-line",
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
    )

    heater = m.fs.heater

    heater.inlet.flow_mol.fix(5102.5)
    heater.inlet.temperature.fix(938.83)
    heater.inlet.pressure.fix(1.2e5)
    heater.inlet.mole_frac_comp[0, "H2"].fix(0.57375)
    heater.inlet.mole_frac_comp[0, "H2O"].fix(0.42517)
    heater.inlet.mole_frac_comp[0, "Ar"].fix(0.00086358)
    heater.inlet.mole_frac_comp[0, "N2"].fix(0.00021589)

    heater.di_tube.fix(0.0525018)
    heater.thickness_tube.fix(0.0039116)
    heater.pitch_x.fix(0.1)
    heater.pitch_y.fix(0.1)
    heater.length_tube_seg.fix(10)
    heater.number_passes.fix(1)
    heater.rfouling = 0.0001
    heater.fcorrection_htc_shell.fix(1)
    heater.cp_wall = 502.4
    if pressure_drop:
        heater.fcorrection_dp_shell.fix(1)

    heater.number_columns_per_pass.fix(40)
    heater.number_rows_per_pass.fix(40)
    heater.electric_heat_duty.fix(3.6504e06)

    pp = m.fs.h2_side_prop_params
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

    shell = heater.control_volume
    iscale.set_scaling_factor(shell.area, 1e-1)
    iscale.set_scaling_factor(shell.heat, 1e-6)
    iscale.set_scaling_factor(shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(heater.heat_holdup, 1e-8)

    iscale.calculate_scaling_factors(m)

    initializer = m.fs.heat_exchanger.default_initializer(
            solver="ipopt",
            solver_options=optarg
    )
    initializer.initialize(m.fs.heat_exchanger)

  
  

Heater Geometry
---------------
=========================== =========== =============================================================================================
Variable                    Index Sets  Doc
=========================== =========== =============================================================================================
``number_columns_per_pass`` None        Number of columns of tube per pass
``number_rows_per_pass``    None        Number of rows of tube per pass
``number_passes``           None        Number of tube banks of ``nrow_tube * ncol_inlet`` tubes
``pitch_x``                 None        Distance between tubes parallel to flow, measured from center-of-tube to center-of-tube
``pitch_y``                 None        Distance between tubes perpendicular to flow, measured from center-of-tube to center-of-tube
``length_tube_seg``         None        Length of tube segment perpendicular to flow in each pass
``area_flow_shell``         None        Reference to flow area on control volume
``length_flow_shell``       None        Reference to flow length on control volume
``area_flow_shell_min``     None        Minimum flow area on shell side
``di_tube``                 None        Inner diameter of tubes
``thickness_tube``          None        Thickness of tube wall.
=========================== =========== =============================================================================================

============================ =========== ===========================================================================
Expression                   Index Sets  Doc
============================ =========== ===========================================================================
``nrow_tube``                None        Total number of rows of tube
``do_tube``                  None        Outer diameter of tube (equal to ``di_tube+2*thickness_tube``)
``pitch_x_to_do``            None        Ratio of ``pitch_x`` to ``do_tube``
``pitch_y_to_do``            None        Ratio of ``pitch_y`` to ``do_tube``
``area_wall_seg``            None        Total cross-sectional area of tube per pass
``total_heat_transfer_area`` None        Total heat transfer area, as measured on outer surface of tubes
============================ =========== ===========================================================================

=========================== =========== =================================================================================================
Constraint                  Index Sets  Doc
=========================== =========== =================================================================================================
``length_flow_shell_eqn``   None        Constrains flow length from control volume to equal value implied by geometry
``area_flow_shell_eqn``     None        Constrains flow cross-sectional area from control volume to equal value implied by geometry
``area_flow_shell_min_eqn`` None        Constraints ``area_flow_shell_min`` to equal value determined by geometry
=========================== =========== =================================================================================================

Performance Equations
-----------------------

================================== ============ =================================================================================
Variable                           Index Sets   Doc
================================== ============ =================================================================================
``electric_heat_duty``             time         Electric heat duty supplied to entire heater unit
``fcorrection_htc_shell``          time, length Correction factor for convective heat transfer
``conv_heat_transfer_coeff_shell`` time, length Convective heat transfer coefficient
``temp_wall_shell``                time, length Wall temperature of tube
``temp_wall_center``               time, length Temperature at center of tube wall
``v_shell``                        time, length Flow velocity through minimum area
``N_Re_shell``                     time, length Reynolds number
``N_Nu_shell``                     time, length Nusselt number
================================== ============ =================================================================================

=========================== =========== =================================================================================
Parameter                   Index Sets  Doc
=========================== =========== =================================================================================
``therm_cond_wall``         None        Thermal conductivity of tube wall
``density_wall``            None        Mass density of tube wall metal
``cp_wall``                 None        Tube wall heat capacity (mass basis)
``rfouling_shell``          None        Fouling resistance on shell side
``f_arrangement``           None        Adjustment factor depending on ``tube_arrangement`` in config
=========================== =========== =================================================================================

====================================== ============ =================================================================================
Constraint                             Index Sets   Doc
====================================== ============ =================================================================================
``v_shell_eqn``                        time, length Calculates velocity of flow through shell using ``area_flow_shell_min``
``N_Re_shell_eqn``                     time, length Calculates the Reynolds number
``conv_heat_transfer_coeff_shell_eqn`` time, length Calculates the convective heat transfer coefficient 
``N_Nu_shell_eqn``                     time, length Calculate the Nusselt number
``heat_shell_eqn``                     time, length Calculates heat transfer per unit length
``temp_wall_shell_eqn``                time, length Calculate the wall temperature of the outer tube
``temp_wall_center_eqn``               time, length Overall energy balance on tube metal
====================================== ============ =================================================================================

====================================== ============ ===================================================================================
Expression                             Index Sets   Doc
====================================== ============ ===================================================================================
``total_heat_transfer_coeff_shell``    time         Returns ``conv_heat_transfer_coeff_shell``. Could be extended to include radiation.
====================================== ============ ===================================================================================

Pressure Change Equations
-------------------------

=========================== ============ =================================================================================
Parameter                    Index Sets   Doc
=========================== ============ =================================================================================
``fcorrection_dp_shell``    None         Correction factor for pressure drop
=========================== ============ =================================================================================

=========================== ============ =================================================================================
Variable                    Index Sets   Doc
=========================== ============ =================================================================================
``fcorrection_dp_shell``    None         Correction factor for pressure drop
``friction_factor_shell``   time, length Friction factor
=========================== ============ =================================================================================

================================== ============ =================================================================================
Constraint                         Index Sets   Doc
================================== ============ =================================================================================
``friction_factor_shell_eqn``      time, length Calculates the friction factor
``deltaP_shell_eqn``               time, length Sets ``deltaP_shell`` based on the friction factor and physical properties
================================== ============ =================================================================================


Holdup Equations
----------------

Created when ``has_holdup=True`` in the config.

=========================== ============ =================================================================================
Variable                    Index Sets   Doc
=========================== ============ =================================================================================
``heat_holdup``             time, length Energy holdup per unit length of flow path
=========================== ============ =================================================================================

=========================== ============ =================================================================================
Constraint                  Index Sets   Doc
=========================== ============ =================================================================================
``heat_holdup_eqn``         time, length Defines heat holdup in terms of geometry and physical properties
=========================== ============ =================================================================================

Dynamic Equations
-----------------

Created when ``dynamic=True`` in the config.

=========================== ============ =================================================================================
Derivative Variable         Index Sets   Doc
=========================== ============ =================================================================================
``heat_accumulation``       time, length Energy accumulation in tube wall per unit length of shell flow path per unit time
=========================== ============ =================================================================================


Initialization
--------------

A simple initialization method that first initializes the control volume without heat transfer,
then adds heat transfer in and solves it again, then finally solves the entire model.


Heater1D Class
--------------

.. autoclass:: Heater1D
  :members:

Heater1DData Class
------------------

.. autoclass:: Heater1DData
  :members:

Heater1DInitializer Class
-------------------------

.. autoclass:: Heater1DInitializer
  :members: