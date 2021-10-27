##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
This is an example flowsheet for the discharge cycle of a natural gas fuel cell
(NGFC) power plant with carbon capture integrated with compressed air energy
storage. The model uses a reduced order model created by PNNL to calculate
the performance of the solid oxide fuel cell.
During the discharge cycle, the steam cycle is turned off.
A compressed gas tank model with fixed volume is used to model gas storage.
"""

# Import Pyomo libraries
from pyomo.environ import (
    TransformationFactory,
    value,
    Var,
    Constraint,
    Block)
from pyomo.network import Arc
from pyomo.opt import TerminationCondition

# IDAES Imports
from idaes.core import MaterialBalanceType
from idaes.core.util import get_solver
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock)
from idaes.generic_models.unit_models import (
    Heater,
    HeatExchanger,
    PressureChanger,
    Valve,
    ValveFunctionType)
from idaes.generic_models.unit_models.pressure_changer import \
    ThermodynamicAssumption
from idaes.generic_models.unit_models.heat_exchanger import (
    delta_temperature_underwood_callback)
import idaes.core.util.unit_costing as icost
import idaes.logger as idaeslog
from idaes.power_generation.flowsheets.\
    sofc.properties.natural_gas_PR_scaled_units import (
        get_NG_properties)
from idaes.power_generation.flowsheets.sofc import sofc as SOFC

# Import tank model
from idaes.power_generation.flowsheets.sofc.compressed_gas_tank import (
    CompressedGasTank
)


def build_model(m):

    # create property packages for air
    air_config = get_NG_properties(
        components=['H2O', 'CO2', 'N2', 'O2', 'Ar'])
    m.fs.air_props = GenericParameterBlock(default=air_config)

    # create unti model instances
    m.fs.air_turbine = PressureChanger(
        default={
            "property_package": m.fs.air_props,
            "compressor": False,
            "material_balance_type": MaterialBalanceType.componentTotal,
            "thermodynamic_assumption": ThermodynamicAssumption.isentropic,
        }
    )

    m.fs.hrsg_heater = Heater(
        default={
            "dynamic": False,
            "property_package": m.fs.air_props,
            "has_pressure_change": True
        }
    )

    m.fs.air_preheater = HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {
                "property_package": m.fs.air_props,
                "material_balance_type": MaterialBalanceType.componentTotal,
                "has_pressure_change": True,
            },
            "tube": {
                "property_package": m.fs.air_props,
                "material_balance_type": MaterialBalanceType.componentTotal,
                "has_pressure_change": True,
            },
        }
    )

    m.fs.storage_tank = CompressedGasTank(
        default={
            "property_package": m.fs.air_props,
            "dynamic": False
        }
    )

    m.fs.flow_valve = Valve(
        default={
            "valve_function_callback": ValveFunctionType.linear,
            "property_package": m.fs.air_props,
            }
    )
    # deactivate pressure-flow correlation to control pressure and flow
    m.fs.flow_valve.pressure_flow_equation.deactivate()

    # create arcs

    # tank to valve
    m.fs.tank_to_valve = Arc(
        source=m.fs.storage_tank.outlet,
        destination=m.fs.flow_valve.inlet
    )
    # valve to preheater tube
    m.fs.valve_to_preheater = Arc(
        source=m.fs.flow_valve.outlet,
        destination=m.fs.air_preheater.inlet_2
    )
    # preheater tube to HRSG
    m.fs.preheater_to_hrsg = Arc(
        source=m.fs.air_preheater.outlet_2,
        destination=m.fs.hrsg_heater.inlet
    )
    # HRSG to turbine
    m.fs.hrsg_to_turbine = Arc(
        source=m.fs.hrsg_heater.outlet,
        destination=m.fs.air_turbine.inlet
    )
    # Turbine to preheater shell
    m.fs.turbine_to_preheater = Arc(
        source=m.fs.air_turbine.outlet,
        destination=m.fs.air_preheater.inlet_1
    )

    TransformationFactory("network.expand_arcs").apply_to(m)


def set_inputs(m):

    # *********** Preheater INPUTS ***********
    m.fs.air_preheater.area.fix(5000)  # m2
    m.fs.air_preheater.overall_heat_transfer_coefficient.fix(0.025)  # kW/m2/K
    m.fs.air_preheater.tube.deltaP[0].fix(-10)  # 0.1 bar = 100 kPa
    m.fs.air_preheater.shell.deltaP[0].fix(-10)  # 0.1 bar = 100 kPa

    # *********** HRSG INPUTS ***********
    m.fs.hrsg_heater.control_volume.properties_out[0].\
        temperature.fix(1173)  # K
    m.fs.hrsg_heater.deltaP.fix(10)  # 0.1 bar = 10 kPa

    # *********** Turbine INPUTS ***********
    m.fs.air_turbine.control_volume.properties_out[0].\
        pressure.fix(150)  # kPa
    m.fs.air_turbine.efficiency_isentropic.fix(0.9)

    # *********** Tank INPUTS ***********
    # Fix tank geometry
    m.fs.storage_tank.tank_diameter.fix(8)
    m.fs.storage_tank.tank_length.fix(20)

    # Fix initial state of tank
    m.fs.storage_tank.previous_state[0].temperature.fix(323)  # K
    m.fs.storage_tank.previous_state[0].pressure.fix(7e3)  # kPa
    m.fs.storage_tank.previous_state[0].\
        mole_frac_comp['H2O'].fix(0.0104)
    m.fs.storage_tank.previous_state[0].\
        mole_frac_comp['N2'].fix(0.7722)
    m.fs.storage_tank.previous_state[0].\
        mole_frac_comp['O2'].fix(0.2077)
    m.fs.storage_tank.previous_state[0].\
        mole_frac_comp['Ar'].fix(0.0094)

    # Fix inlet state of tank
    m.fs.storage_tank.control_volume.properties_in[0].\
        flow_mol.fix(5)  # kmol/s
    m.fs.storage_tank.control_volume.properties_in[0].\
        temperature.fix(323)  # K
    m.fs.storage_tank.control_volume.properties_in[0].\
        pressure.fix(7e3)  # kPa
    m.fs.storage_tank.control_volume.properties_in[0].\
        mole_frac_comp['H2O'].fix(0.0104)
    m.fs.storage_tank.control_volume.properties_in[0].\
        mole_frac_comp['CO2'].fix(0.0003)
    m.fs.storage_tank.control_volume.properties_in[0].\
        mole_frac_comp['N2'].fix(0.7722)
    m.fs.storage_tank.control_volume.properties_in[0].\
        mole_frac_comp['O2'].fix(0.2077)
    m.fs.storage_tank.control_volume.properties_in[0].\
        mole_frac_comp['Ar'].fix(0.0094)

    # Fix Duration of Operation (Time Step =  1hr)
    m.fs.storage_tank.dt[0].fix(3600)
    m.fs.storage_tank.control_volume.properties_out[0].\
        flow_mol.fix(0.5)  # kmol/s

    # Fix outlet pressure for Valve
    m.fs.flow_valve.outlet.pressure[0].fix(4e3)


def add_bounds(m):

    # Add required bounds to the variables

    m.fs.air_turbine.control_volume.properties_in[0].\
        pressure.setub(1e15)
    m.fs.air_turbine.control_volume.properties_out[0].\
        pressure.setub(1e15)

    m.fs.hrsg_heater.control_volume.properties_in[0].\
        pressure.setub(1e15)
    m.fs.hrsg_heater.control_volume.properties_out[0].\
        pressure.setub(1e15)

    m.fs.air_preheater.tube.properties_in[0].\
        pressure.setub(1e15)
    m.fs.air_preheater.tube.properties_out[0].\
        pressure.setub(1e15)

    m.fs.flow_valve.control_volume.properties_in[0].\
        pressure.setub(1e15)
    m.fs.flow_valve.control_volume.properties_out[0].\
        pressure.setub(1e15)

    # Setting the bounds on the state variables
    m.fs.storage_tank.control_volume.properties_in[0].\
        pressure.setub(1e15)
    m.fs.storage_tank.control_volume.properties_out[0].\
        pressure.setub(1e15)
    m.fs.storage_tank.previous_state[0].pressure.setub(1e15)


def initialize(m, outlvl=idaeslog.NOTSET, solver=None, optarg=None):

    optarg = {
        "max_iter": 300,
        "halt_on_ampl_error": "yes",
        'bound_push': 1e-23,
        "mu_init": 1e-6,
    }
    solver = get_solver(solver, optarg)

    iscale.calculate_scaling_factors(m)

    # initialize tank
    m.fs.storage_tank.\
        initialize(outlvl=outlvl, optarg=solver.options)
    solver.solve(m.fs.storage_tank)
    m.fs.storage_tank.volume_cons.deactivate()
    m.fs.storage_tank.control_volume.volume[:].fix(300000)
    m.fs.storage_tank.control_volume.properties_in[0].\
        flow_mol.fix(0)
    m.fs.storage_tank.control_volume.properties_out[0].\
        flow_mol.fix(7)
    solver.solve(m.fs.storage_tank)

    # initialize valve
    propagate_state(m.fs.tank_to_valve)
    m.fs.flow_valve.control_volume.properties_in[0].\
        flow_mol.value = 5
    m.fs.flow_valve.\
        initialize(outlvl=outlvl, optarg=solver.options)

    # initialize preheater
    m.fs.air_preheater.shell.properties_in[0].\
        flow_mol.value = 1  # kmol/s
    m.fs.air_preheater.shell.properties_in[0].\
        temperature.value = 600  # K
    m.fs.air_preheater.shell.properties_in[0].\
        pressure.value = 1740  # kPa
    m.fs.air_preheater.shell.properties_in[0].\
        mole_frac_comp['H2O'].value = 0.0104
    m.fs.air_preheater.shell.properties_in[0].\
        mole_frac_comp['CO2'].value = 0.0003
    m.fs.air_preheater.shell.properties_in[0].\
        mole_frac_comp['N2'].value = 0.7722
    m.fs.air_preheater.shell.properties_in[0].\
        mole_frac_comp['O2'].value = 0.2077
    m.fs.air_preheater.shell.properties_in[0].\
        mole_frac_comp['Ar'].value = 0.0094
    propagate_state(m.fs.valve_to_preheater)
    m.fs.air_preheater.\
        initialize(outlvl=outlvl, optarg=solver.options)

    # initialize HRSG
    propagate_state(m.fs.preheater_to_hrsg)
    m.fs.hrsg_heater.\
        initialize(outlvl=outlvl, optarg=solver.options)

    # initialize turbine
    propagate_state(m.fs.hrsg_to_turbine)
    m.fs.air_turbine.\
        initialize(outlvl=outlvl, optarg=solver.options)

    # update conditions for discharge scenario
    m.fs.storage_tank.control_volume.properties_in[0].\
        flow_mol.fix(0)  # kmol/s

    res = solver.solve(m.fs)
    if res.solver.termination_condition == TerminationCondition.optimal:
        print('Flowsheet Initialization Successful.')
    else:
        print('Flowsheet Initialization Failed.')

    print("*************  Caes storage block Initialized   **************")
    print("      ")


def build_costing(m):

    # Chemical engineering cost index for 2019
    m.CE_index = 607.5  
    # Number of years for annulaizing capital cost
    m.number_of_years = 15

    m.fs.discharge_capital_cost = Var(
        initialize=1000000,
        doc="Annualized capital cost")

    m.fs.air_turbine.get_costing()
    m.fs.air_turbine.costing.CE_index = m.CE_index
    icost.initialize(m.fs.air_turbine.costing)

    m.fs.air_preheater.get_costing()
    m.fs.air_preheater.costing.CE_index = m.CE_index
    icost.initialize(m.fs.air_preheater.costing)

    # creating a dummy block to add capital costs associated with gas storage
    # this is a placeholder and should be updated with actual cost data
    m.fs.storage_tank.costing = Block()
    m.fs.storage_tank.purchase_cost = Var(
    initialize=100000,
    doc="Capital costs associated gas storage")
    m.fs.storage_tank.purchase_cost.fix()

    # No costing added for the interstage_cooler block as this is a dummy unit

    # Computing annualized capital cost
    def discharge_capital_cost_rule(b):
        return m.fs.discharge_capital_cost == ((
            m.fs.air_turbine.costing.purchase_cost +
            m.fs.air_preheater.costing.purchase_cost +
            m.fs.storage_tank.purchase_cost) / m.number_of_years
        )
    m.fs.discharge_capital_cost_eq = Constraint(rule=discharge_capital_cost_rule)


def build_caes_discharge(m):

    # build the storage model, set model inputs, and add bounds
    build_model(m)
    set_inputs(m)
    add_bounds(m)

    # check if the storage model is complete by asserting DOF = 0
    assert degrees_of_freedom(m) == 0

    # initialize the charge storage model
    initialize(m)
    assert degrees_of_freedom(m) == 0

    return m


def integrate_storage(m):

    # Deactivating constraints to turn off steams cycle and integrate storage
    m.fs.steam_cycle_heat_constraint.deactivate()
    m.fs.steam_cycle_power_constraint.deactivate()
    m.fs.steam_cycle_loss_constraint.deactivate()
    m.fs.gross_power_constraint.deactivate()
    m.fs.feedwater_pump_work_constraint.deactivate()
    m.fs.condensate_pump_work_constraint.deactivate()
    m.fs.steam_turbine_auxiliary_constraint.deactivate()
    m.fs.net_power_constraint.deactivate()

    m.fs.steam_cycle_heat.fix()
    m.fs.steam_cycle_power.fix()
    m.fs.steam_cycle_loss.fix()
    m.fs.feedwater_pump_work.fix()
    m.fs.condensate_pump_work.fix()
    m.fs.steam_turbine_auxiliary.fix()

    @m.fs.Constraint(m.fs.time)
    def gross_power_constraint_mod(fs, t):
        return (fs.gross_power[t] == fs.stack_power_AC[t]
                + fs.air_turbine.control_volume.work[t]*-1e-3
                )

    @m.fs.Constraint(m.fs.time)
    def steam_cycle_heat_constraint_mod(fs, t):
        return (0 <=
                fs.HRSG_heat_duty[t] -
                fs.ASU_HP_steam_heat[t] -
                fs.hrsg_heater.heat_duty[0]*1e-3
                )

    @m.fs.Constraint(m.fs.time)
    def net_power_constraint_mod(fs, t):
        return (fs.net_power[t] ==
                fs.gross_power[t] -
                fs.auxiliary_load[t] -
                fs.transformer_losses[t])


def main():
    # create charge model
    # m = ConcreteModel()
    # m.fs = FlowsheetBlock(default={"dynamic": False})

    # build SOFC model
    m = SOFC.get_model()

    # build storage model
    m = build_caes_discharge(m)

    # integrate storage with SOFC flowsheet
    integrate_storage(m)

    # add costing
    build_costing(m)

    return m


if __name__ == "__main__":

    optarg = {
        "max_iter": 300,
        "halt_on_ampl_error": "yes",
        'bound_push': 1e-23,
        "mu_init": 1e-6,
    }

    solver = get_solver("ipopt", optarg)

    m = main()

    m.fs.storage_tank.control_volume.\
        properties_out[0].flow_mol.unfix()  # kmol/s
    m.fs.flow_valve.control_volume.\
        properties_out[0].flow_mol.fix(6)  # kmol/s

    solver.solve(m, tee=True, symbolic_solver_labels=True)

    print("              ")
    print("Turbine Power", value(
        m.fs.air_turbine.control_volume.work[0])*-1e-3)
    print("Heater Duty", value(
        m.fs.hrsg_heater.heat_duty[0])*1e-3)
    print("Net Power", value(m.fs.net_power[0]))
