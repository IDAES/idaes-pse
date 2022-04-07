Heat Exchanger (Lumped Capacitance)
===================================

.. index::
   pair: idaes.models.unit_models.heat_exchanger_lc;HeatExchangerLumpedCapacitance

.. currentmodule:: idaes.models.unit_models.heat_exchanger_lc

Author: `Rusty Gentile <https://github.com/rustygentile>`_

The ``HeatExchangerLumpedCapacitance`` model can be imported from :code:`idaes.models.unit_models`.
This model is an extension of ``idaes.models.unit_models.heat_exchanger``,
with wall temperature and heat holdup terms added for use in transient simulations. Using the electric 
circuit analogy, heat stored in the wall material is similar to charge in a capacitor.

Degrees of Freedom
------------------

``HeatExchangerLumpedCapacitance`` has three degrees of freedom. For most scenarios, these will be the wall 
temperature and exit temperatures on the hot and cold sides. The constraints for this model 
should be similar to those of the standard 
:ref:`Heat Exchanger <reference_guides/model_libraries/generic/unit_models/heat_exchanger:HeatExchanger (0D)>` 
model. In lieu of "heat_transfer_coefficient", however, the following should be fixed:

* ua_cold_side
* ua_hot_side

The user may also add constraints to make these functions of fluid state.

Model Structure
---------------

Structure is derived from the base 0D :ref:`Heat Exchanger <reference_guides/model_libraries/generic/unit_models/heat_exchanger:Model Structure>` 
model. Aside from adding new variables and constraints, the structure of this model is unchanged.

Variables
---------

All variables from the base 0D :ref:`Heat Exchanger <reference_guides/model_libraries/generic/unit_models/heat_exchanger:Variables>` 
model are included. This model contains the following variables in addition. Each of these is indexed in the time domain.

=========================== ============================= ==========================================================================================
Variable                    Symbol                        Doc
=========================== ============================= ==========================================================================================
temperature_wall            :math:`T_{wall}`              Average wall temperature
dT_wall_dt                  :math:`\frac{dT_{wall}}{dt}`  Derivative of wall temperature with respect to time 
ua_cold_side                :math:`UA_{cold}`             Overall heat transfer coefficient from the cold side
ua_hot_side                 :math:`UA_{hot}`              Overall heat transfer coefficient from the hot side
ua_hot_side_to_wall         :math:`UA_{hot, wall}`        Overall heat transfer coefficient between the hot side and the center of the wall material
hot_side_heat               :math:`Q_{hot}`               Heat entering the hot side control volume (value is typically negative)
cold_side_heat              :math:`Q_{cold}`              Heat entering the cold side control volume (value is typically positive)
=========================== ============================= ==========================================================================================

Note: "ua" variables are of the form :math:`hconv \times area`. So the units of measure of these should be :math:`\frac{energy}{time \times temperature}`

Parameters
----------

The model uses the following parameters. These are not indexed as we assume they are constant.

=========================== ====================== ==========================================================================================
Parameter                    Symbol                 Doc
=========================== ====================== ==========================================================================================
heat_capacity_wall          :math:`C_{wall}`       Total heat capacity of heat exchanger material
thermal_resistance_wall     :math:`R_{wall}`       Total thermal resistance of heat exchanger material
thermal_fouling_cold_side   :math:`R_{foul, cold}` Total thermal resistance due to fouling on the cold side
thermal_fouling_hot_side    :math:`R_{foul, hot}`  Total thermal resistance due to fouling on the hot side
=========================== ====================== ==========================================================================================

Constraints
-----------

Note: all constraints from the base 0D :ref:`Heat Exchanger <reference_guides/model_libraries/generic/unit_models/heat_exchanger:Constraints>` 
model are also included.

Total overall heat transfer coefficient:

.. math::
  \frac{1}{UA} = \frac{1}{UA_{hot}} + R_{foul, hot} + R_{wall} + R_{foul, cold} + \frac{1}{UA_{cold}}

Overall heat transfer coefficient from the hot side to the center of the wall:

.. math::
  \frac{1}{UA_{hot, wall}} = \frac{1}{UA_{hot}} + R_{foul, hot} + \frac{1}{2}R_{wall}

Wall temperature equation:

.. math::
  T_{wall} = \frac{1}{2}(T_{hot, in} + T_{hot, out}) + \frac{Q_{hot}}{UA_{hot, wall}}

Dynamic heat balance (if ``dynamic_heat_balance`` is set to ``True``):

.. math::
  Q_{hot} + Q_{cold} + \frac{dT_{wall}}{dt}C_{wall} = 0

Standard heat balance (if ``dynamic_heat_balance`` is set to ``False``):

.. math::
  Q_{hot} + Q_{cold} = 0

Example
-------

In this example, we run a transient simulation of an air-cooled SCO2 heat exchanger using a dynamic 
heat balance. The control volumes, however, are static. This is essentially an assumption that any 
mass holdup of the working fluid is negligible.

To use dynamic control volumes, constraints should be added that relate flow rate to 
pressure loss. In that case, the :ref:`volume <reference_guides/core/control_volume_0d:add_geometry>` 
variables of the hot and cold sides should also be fixed.

.. testcode::
  
  import pyomo.environ as pe
  from idaes.core import FlowsheetBlock
  from idaes.models.unit_models import HeatExchangerLumpedCapacitance, HeatExchangerFlowPattern
  from idaes.models.unit_models.heat_exchanger import delta_temperature_lmtd_callback
  from idaes.models.properties import swco2
  from idaes.models_extra.power_generation.properties import FlueGasParameterBlock
  from idaes.core.util.plot import plot_grid_dynamic
  
  m = pe.ConcreteModel()
  m.fs = FlowsheetBlock(default={"dynamic": True,
                                 "time_set": [0, 300, 600, 900, 1200, 1500],
                                 "time_units": pe.units.s})
  
  m.fs.prop_sco2 = swco2.SWCO2ParameterBlock()
  m.fs.prop_fluegas = FlueGasParameterBlock()
  
  m.fs.HE = HeatExchangerLumpedCapacitance(default={
      "delta_temperature_callback": delta_temperature_lmtd_callback,
      "cold_side_name": "shell",
      "hot_side_name": "tube",
      "shell": {"property_package": m.fs.prop_fluegas,
                "has_pressure_change": False},
      "tube": {"property_package": m.fs.prop_sco2,
               "has_pressure_change": True},
      "flow_pattern": HeatExchangerFlowPattern.crossflow,
      "dynamic_heat_balance": True,
      "dynamic": False})
  
  m.discretizer = pe.TransformationFactory('dae.finite_difference')
  m.discretizer.apply_to(m, nfe=200, wrt=m.fs.time, scheme="BACKWARD")
  
  # Cold-side boundary conditions
  shell_inlet_temperature = 288.15
  shell_flow = 44004.14222
  shell_area = 690073.9153
  shell_hconv = 22
  m.fs.HE.ua_cold_side[:].fix(shell_area * shell_hconv)
  m.fs.HE.shell_inlet.flow_mol_comp[:, "H2O"].fix(0.01027 * shell_flow)
  m.fs.HE.shell_inlet.flow_mol_comp[:, "CO2"].fix(0.000411592 * shell_flow)
  m.fs.HE.shell_inlet.flow_mol_comp[:, "N2"].fix(0.780066026 * shell_flow)
  m.fs.HE.shell_inlet.flow_mol_comp[:, "O2"].fix(0.209252382 * shell_flow)
  m.fs.HE.shell_inlet.flow_mol_comp[:, "NO"].fix(0)
  m.fs.HE.shell_inlet.flow_mol_comp[:, "SO2"].fix(0)
  m.fs.HE.shell_inlet.temperature[:].fix(shell_inlet_temperature)
  m.fs.HE.shell_inlet.pressure[:].fix(101325)
  
  # Hot-side boundary conditions
  tube_inlet_temperature = 384.35
  tube_inlet_pressure = 7653000
  tube_outlet_pressure = 7500000
  tube_inlet_enthalpy = swco2.htpx(T=tube_inlet_temperature * pe.units.K,
                                   P=tube_inlet_pressure * pe.units.Pa)
  tube_flow = 13896.84163
  tube_area = 19542.2771
  tube_hconv = 1000
  m.fs.HE.ua_hot_side[:].fix(tube_area * tube_hconv)
  m.fs.HE.tube_inlet.flow_mol[:].fix(tube_flow)
  m.fs.HE.tube_inlet.pressure[:].fix(tube_inlet_pressure)
  m.fs.HE.tube_inlet.enth_mol[:].fix(tube_inlet_enthalpy)
  m.fs.HE.tube_outlet.pressure[:].fix(tube_outlet_pressure)
  
  # number of tubes * tube mass * specific heat capacity
  m.fs.HE.heat_capacity_wall = 1160 * 322 * 466
  m.fs.HE.crossflow_factor.fix(0.8)
  
  # Area has no effect on the result so long as it isn't zero
  m.fs.HE.area.fix(1)
  
  # Initialize the model with steady-state boundary conditions
  m.fs.HE.initialize()
  solver = pe.SolverFactory('ipopt')
  
  # Activate the heat holdup term and add temperature disturbances
  m.fs.HE.activate_dynamic_heat_eq()
  for t in m.fs.time:
      if 300 <= t < 600:
          m.fs.HE.shell_inlet.temperature[t].fix(288.15 - 10)
      elif 600 <= t < 900:
          m.fs.HE.shell_inlet.temperature[t].fix(288.15)
      elif 900 <= t < 1200:
          m.fs.HE.shell_inlet.temperature[t].fix(288.15 + 10)
      elif t >= 1200:
          m.fs.HE.shell_inlet.temperature[t].fix(288.15)
  
  solver.solve(m)
  
  # Plot some results and save the figure to a file.
  plot_grid_dynamic(
      x=m.fs.time,
      xlabel="time (s)",
      y=[
          m.fs.HE.tube.properties_out[:].temperature,
          m.fs.HE.tube.properties_out[:].dens_mass,
          m.fs.HE.temperature_wall[:],
          m.fs.HE.shell.properties_in[:].temperature
      ],
      ylabel=[
          "sCO2 out (K)",
          "sCO2 out (kg/m^3)",
          "wall temperature (K)",
          "air in (K)"
      ],
      yunits=[
          pe.units.K,
          pe.units.kg / pe.units.m ** 3,
          pe.units.K,
          pe.units.K
      ],
      cols=2,
      rows=2,
      to_file="transient_sco2.pdf"
  )
  

Class Documentation
-------------------

.. autoclass:: HeatExchangerLumpedCapacitance
   :members:

.. autoclass:: HeatExchangerLumpedCapacitanceData
   :members:
