Turbine (Multistage)
====================

.. index::
  pair: idaes.power_generation.unit_models.helm.turbine_multistage;HelmTurbineMultistage

.. module:: idaes.power_generation.unit_models.helm.turbine_multistage

This is a composite model for a power plant turbine with high, intermediate and
low pressure sections. This model contains an inlet stage with throttle valves
for partial arc admission and optional splitters for steam extraction.

The figure below shows the layout of the mutistage turbine model.  Optional splitters
provide for steam extraction.  The splitters can have two or more outlets (one being
the main steam outlet).  The streams that connect one stage to the next can also be
omitted.  This allows for connecting additional unit models (usually reheaters) between
stages.

.. figure:: turbine_multistage.svg
  :width: 800
  :align: center

  MultiStage Turbine Model

Example
-------

This example sets up a turbine multistage turbine model similar to what could be
found in a power plant steam cycle.  There are 7 high-pressure stages, 14
intermediate-pressure stages, and 11 low-pressure stages. Steam extractions
are provided after stages hp4, hp7, ip5, ip14, lp4, lp7, lp9, lp11.  The
extraction at ip14 uses a splitter with three outlets, one for the main steam,
one for the boiler feed pump, and one for a feedwater heater.  There is a
disconnection between the HP and IP sections so that steam can be sent to a
reheater.  In this example, a heater block is a stand-in for a reheater model.

.. code-block:: python

  from pyomo.environ import (ConcreteModel, SolverFactory, TransformationFactory,
                             Constraint, value)
  from pyomo.network import Arc

  from idaes.core import FlowsheetBlock
  from idaes.unit_models import Heater
  from idaes.power_generation.unit_models.helm import HelmTurbineMultistage
  from idaes.generic_models.properties import iapws95

  solver = SolverFactory('ipopt')
  solver.options = {'tol': 1e-6}

  m = ConcreteModel()
  m.fs = FlowsheetBlock(default={"dynamic": False})
  m.fs.properties = iapws95.Iapws95ParameterBlock()
  m.fs.turb = HelmTurbineMultistage(default={
      "property_package": m.fs.properties,
      "num_hp": 7,
      "num_ip": 14,
      "num_lp": 11,
      "hp_split_locations": [4,7],
      "ip_split_locations": [5, 14],
      "lp_split_locations": [4,7,9,11],
      "hp_disconnect": [7], # 7 is last stage in hp so disconnect hp from ip
      "ip_split_num_outlets": {14:3}})
  # Add reheater (for example using a simple heater block)
  m.fs.reheat = Heater(default={"property_package": m.fs.properties})
  # Add Arcs (streams) to connect the HP and IP sections through reheater
  m.fs.hp_to_reheat = Arc(source=m.fs.turb.hp_split[7].outlet_1,
                          destination=m.fs.reheat.inlet)
  m.fs.reheat_to_ip = Arc(source=m.fs.reheat.outlet,
                          destination=m.fs.turb.ip_stages[1].inlet)
  # Set the turbine inlet conditions and an initial flow guess
  p = 2.4233e7
  hin = iapws95.htpx(T=880, P=p)
  m.fs.turb.inlet_split.inlet.enth_mol[0].fix(hin)
  m.fs.turb.inlet_split.inlet.flow_mol[0].fix(26000)
  m.fs.turb.inlet_split.inlet.pressure[0].fix(p)

  # Set the inlet of the ip section for initialization, since it is disconnected
  p = 7.802e+06
  hin = iapws95.htpx(T=880, P=p)
  m.fs.turb.ip_stages[1].inlet.enth_mol[0].value = hin
  m.fs.turb.ip_stages[1].inlet.flow_mol[0].value = 25220.0
  m.fs.turb.ip_stages[1].inlet.pressure[0].value = p
  # Set the efficency and pressure ratios of stages other than inlet and outlet
  for i, s in turb.hp_stages.items():
      s.ratioP[:] = 0.88
      s.efficiency_isentropic[:] = 0.9
  for i, s in turb.ip_stages.items():
      s.ratioP[:] = 0.85
      s.efficiency_isentropic[:] = 0.9
  for i, s in turb.lp_stages.items():
      s.ratioP[:] = 0.82
      s.efficiency_isentropic[:] = 0.9
  # Usually these fractions would be determined by the boiler feed water heater
  # network. Since this example doesn't include them, just fix split fractions
  turb.hp_split[4].split_fraction[0,"outlet_2"].fix(0.03)
  turb.hp_split[7].split_fraction[0,"outlet_2"].fix(0.03)
  turb.ip_split[5].split_fraction[0,"outlet_2"].fix(0.04)
  turb.ip_split[14].split_fraction[0,"outlet_2"].fix(0.04)
  turb.ip_split[14].split_fraction[0,"outlet_3"].fix(0.15)
  turb.lp_split[4].split_fraction[0,"outlet_2"].fix(0.04)
  turb.lp_split[7].split_fraction[0,"outlet_2"].fix(0.04)
  turb.lp_split[9].split_fraction[0,"outlet_2"].fix(0.04)
  turb.lp_split[11].split_fraction[0,"outlet_2"].fix(0.04)
  # unfix inlet flow for pressure driven simulation
  turb.inlet_split.inlet.flow_mol.unfix()
  # Set the inlet steam mixer to use the constraints that the pressures of all
  # inlet streams are equal
  turb.inlet_mix.use_equal_pressure_constraint()
  # Initialize turbine
  turb.initialize(outlvl=1)
  # Copy conditions out of turbine to initialize the reheater
  for t in m.fs.time:
      m.fs.reheat.inlet.flow_mol[t].value = \
          value(turb.hp_split[7].outlet_1_state[t].flow_mol)
      m.fs.reheat.inlet.enth_mol[t].value = \
          value(turb.hp_split[7].outlet_1_state[t].enth_mol)
      m.fs.reheat.inlet.pressure[t].value = \
          value(turb.hp_split[7].outlet_1_state[t].pressure)
  # initialize the reheater
  m.fs.reheat.initialize(outlvl=4)
  # Add constraint to the reheater to result in 880K outlet temperature
  def reheat_T_rule(b, t):
      return m.fs.reheat.control_volume.properties_out[t].temperature == 880
  m.fs.reheat.temperature_out_equation = Constraint(m.fs.reheat.time_ref,
      rule=reheat_T_rule)
  # Expand the Arcs connecting the turbine to the reheater
  TransformationFactory("network.expand_arcs").apply_to(m)
  # Fix the outlet pressure (usually determined by condenser)
  m.fs.turb.outlet_stage.control_volume.properties_out[0].pressure.fix()

  # Solve the pressure driven flow model with reheat
  solver.solve(m, tee=True)

Unit Models
-----------

The multistage turbine model contains the models in the table below.  The splitters for steam extraction are not present if a turbine section contains no steam extractions.

=========================== ==================== ======================================================================================================================================================
Unit                        Index Sets           Doc
=========================== ==================== ======================================================================================================================================================
``inlet_split``             None                 Splitter to split the main steam feed into steams for each arc (:ref:`Separator <reference_guides/model_libraries/generic/unit_models/separator:Separator>`)
``throttle_valve``          Admission Arcs       Throttle valves for each admission arc (:ref:`HelmValve <reference_guides/model_libraries/power_generation/unit_models/steam_valve:HelmValve>`)
``inlet_stage``             Admission Arcs       Parallel inlet turbine stages that represent admission arcs (:ref:`TurbineInlet <reference_guides/model_libraries/power_generation/unit_models/turbine_inlet:Turbine (Inlet Stage)>`)
``inlet_mix``               None                 Mixer to combine the streams from each arc back to one stream (:ref:`Mixer <reference_guides/model_libraries/generic/unit_models/mixer:Mixer>`)
``hp_stages``               HP stages            Turbine stages in the high-pressure section (:ref:`TurbineStage <reference_guides/model_libraries/power_generation/unit_models/turbine_stage:Turbine (Stage)>`)
``ip_stages``               IP stages            Turbine stages in the intermediate-pressure section (:ref:`TurbineStage <reference_guides/model_libraries/power_generation/unit_models/turbine_stage:Turbine (Stage)>`)
``lp_stages``               LP stages            Turbine stages in the low-pressure section (:ref:`TurbineStage <reference_guides/model_libraries/power_generation/unit_models/turbine_stage:Turbine (Stage)>`)
``hp_splits``               subset of HP stages  Extraction splitters in the high-pressure section (:ref:`Separator <reference_guides/model_libraries/generic/unit_models/separator:Separator>`)
``ip_splits``               subset of IP stages  Extraction splitters in the high-pressure section (:ref:`Separator <reference_guides/model_libraries/generic/unit_models/separator:Separator>`)
``lp_splits``               subset of LP stages  Extraction splitters in the high-pressure section (:ref:`Separator <reference_guides/model_libraries/generic/unit_models/separator:Separator>`)
``outlet_stage``            None                 The final stage in the turbine, which calculates exhaust losses (:ref:`TurbineOutlet <reference_guides/model_libraries/power_generation/unit_models/turbine_outlet:Turbine (Outlet Stage)>`)
=========================== ==================== ======================================================================================================================================================

Initialization
--------------

The initialization approach is to sequentially initialize each sub-unit using the outlet of the previous
model. Before initializing the model, the inlet of the turbine, and any stage that is disconnected should
be given a reasonable guess.  The efficiency and pressure ration of the stages in the HP, IP and LP
sections should be specified. For the inlet and outlet stages the flow coefficient should be specified.
Valve coefficients should also be specified.  A reasonable guess for split fractions should also be given
for any extraction splitters present. The most likely cause of initialization failure is flow coefficients
in inlet stage, outlet stage, or valves that do not pair well with the specified flow rates.

The flow coefficients for the inlet and outlet stage can be difficult to determine, therefore the
initialization arguments ``calculate_outlet_cf`` and ``calculate_outlet_cf`` are provided. If these are
True, the first stage flow coefficient is calculated from the flow and pressure ratio guesses, and the
outlet flow coefficient is calculated from the exhaust pressure and flow.


TurbineMultistage Class
-----------------------

.. autoclass:: HelmTurbineMultistage
  :members:

TurbineMultistageData Class
---------------------------

.. autoclass:: HelmTurbineMultistageData
  :members:
