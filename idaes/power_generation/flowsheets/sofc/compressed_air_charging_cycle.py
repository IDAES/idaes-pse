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
This is an example flowsheet for the charge cycle of a natural gas fuel cell
(NGFC) power plant with carbon capture integrated with compressed air energy
storage. The model uses a reduced order model created by PNNL to calculate
the performance of the solid oxide fuel cell.
The power plant model builds a simple steam cycle using the the efficiency
and energy balance.
For modeling compressed air energy storage, a single stage compressor and a
single cooler for interstage cooling is assumed. A compressed gas tank model
with fixed volume is used to model gas storage.
"""

# Import Pyomo libraries
from pyomo.environ import (
    TransformationFactory,
    ConcreteModel,
    value,
    Var,
    Constraint,
    Block)
from pyomo.network import Arc

# IDAES Imports
import idaes.logger as idaeslog
from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.core.util import get_solver
from idaes.core.util.initialization import propagate_state
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock,
)
from idaes.generic_models.unit_models import Heater, PressureChanger
from idaes.generic_models.unit_models.pressure_changer import \
    ThermodynamicAssumption
import idaes.core.util.unit_costing as icost
from idaes.power_generation.flowsheets.sofc import sofc as SOFC
from idaes.power_generation.flowsheets.\
    sofc.properties.natural_gas_PR_scaled_units import (
    get_NG_properties,
)

# Import tank model
from compressed_gas_tank import CompressedGasTank


def build_model(m):

    # create property packages for air
    air_config = get_NG_properties(components=["H2O", "CO2", "N2", "O2", "Ar"])
    m.fs.air_props = GenericParameterBlock(default=air_config)

    # create unti model instances
    m.fs.air_compressor = PressureChanger(
        default={
            "property_package": m.fs.air_props,
            "compressor": True,
            "material_balance_type": MaterialBalanceType.componentTotal,
            "thermodynamic_assumption": ThermodynamicAssumption.isentropic,
        }
    )

    m.fs.interstage_cooler = Heater(
        default={
            "dynamic": False,
            "property_package": m.fs.air_props,
            "has_pressure_change": False,
        }
    )

    m.fs.storage_tank = CompressedGasTank(
        default={"property_package": m.fs.air_props, "dynamic": False}
    )

    # create arcs

    # compressor to cooler
    m.fs.compressor_to_cooler = Arc(
        source=m.fs.air_compressor.outlet,
        destination=m.fs.interstage_cooler.inlet
    )
    # cooler to tank
    m.fs.cooler_to_tank = Arc(
        source=m.fs.interstage_cooler.outlet,
        destination=m.fs.storage_tank.inlet
    )
    TransformationFactory("network.expand_arcs").apply_to(m)


def set_inputs(m):

    # *********** COMPRESSOR INPUTS ***********
    # Air inlet to compressor
    m.fs.air_compressor.inlet.pressure[0].fix(101.325)  # kPa
    m.fs.air_compressor.inlet.temperature[0].fix(303)  # K
    m.fs.air_compressor.inlet.flow_mol.fix(5)  # kmol/s
    m.fs.air_compressor.inlet.mole_frac_comp[0, "H2O"].fix(0.0104)
    m.fs.air_compressor.inlet.mole_frac_comp[0, "CO2"].fix(0.0003)
    m.fs.air_compressor.inlet.mole_frac_comp[0, "N2"].fix(0.7726)
    m.fs.air_compressor.inlet.mole_frac_comp[0, "O2"].fix(0.2077)
    m.fs.air_compressor.inlet.mole_frac_comp[0, "Ar"].fix(0.0094)

    m.fs.air_compressor.efficiency_isentropic.fix(0.9)
    m.fs.air_compressor.outlet.pressure[0].fix(7e3)  # kPa

    # *********** COOLER INPUTS ***********
    m.fs.interstage_cooler.outlet.temperature[0].fix(323)  # K

    # *********** TANK INPUTS ***********
    # Fix tank geometry
    m.fs.storage_tank.tank_diameter.fix(8)
    m.fs.storage_tank.tank_length.fix(20)

    # Fix initial state of tank
    m.fs.storage_tank.previous_state[0].temperature.fix(323)  # K
    m.fs.storage_tank.previous_state[0].pressure.fix(4e3)  # kPa
    m.fs.storage_tank.previous_state[0].mole_frac_comp["H2O"].fix(0.0104)
    m.fs.storage_tank.previous_state[0].mole_frac_comp["N2"].fix(0.7722)
    m.fs.storage_tank.previous_state[0].mole_frac_comp["O2"].fix(0.2077)
    m.fs.storage_tank.previous_state[0].mole_frac_comp["Ar"].fix(0.0094)

    # Fix Duration of Operation (Time Step =  1hr)
    m.fs.storage_tank.dt[0].fix(3600)

    # Fix the outlet flow to zero: during charge
    m.fs.storage_tank.control_volume.properties_out[0].flow_mol.fix(0)


def add_bounds(m):

    # Add required bounds to the variables

    m.fs.air_compressor.control_volume.properties_out[0].\
        pressure.setub(1e15)
    m.fs.air_compressor.properties_isentropic[0].\
        pressure.setub(1e15)

    m.fs.interstage_cooler.control_volume.properties_in[0].\
        pressure.setub(1e15)
    m.fs.interstage_cooler.control_volume.properties_out[0].\
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
    }
    solver = get_solver(solver, optarg)

    iscale.calculate_scaling_factors(m)

    # initialize compressor
    m.fs.air_compressor.\
        initialize(outlvl=outlvl, optarg=solver.options)

    # initialize cooler
    propagate_state(m.fs.compressor_to_cooler)
    m.fs.interstage_cooler.\
        initialize(outlvl=outlvl, optarg=solver.options)

    # initialize tank
    propagate_state(m.fs.cooler_to_tank)

    m.fs.storage_tank.previous_state[0].\
        mole_frac_comp["CO2"].value = 0.0003
    m.fs.storage_tank.previous_state[0].\
        mole_frac_phase_comp[
        "Vap", "H2O"
    ].value = 0.0104
    m.fs.storage_tank.previous_state[0].\
        mole_frac_phase_comp[
        "Vap", "CO2"
    ].value = 0.0003
    m.fs.storage_tank.previous_state[0].\
        mole_frac_phase_comp["Vap", "N2"].value = 0.7722
    m.fs.storage_tank.previous_state[0].\
        mole_frac_phase_comp["Vap", "O2"].value = 0.2077
    m.fs.storage_tank.previous_state[0].\
        mole_frac_phase_comp["Vap", "Ar"].value = 0.0094
    m.fs.storage_tank.previous_state[0].\
        flow_mol_phase["Vap"].value = 0

    m.fs.storage_tank.\
        initialize(outlvl=outlvl, optarg=solver.options)

    print("*************  Storage Unit Models Initialized   **************")


def build_costing(m):

    # Chemical engineering cost index for 2019
    m.CE_index = 607.5  
    # Number of years for annulaizing capital cost
    m.number_of_years = 15

    m.fs.charge_capital_cost = Var(
        initialize=1000000,
        doc="Annualized capital cost")

    m.fs.air_compressor.get_costing()
    m.fs.air_compressor.costing.CE_index = m.CE_index
    icost.initialize(m.fs.air_compressor.costing)

    # creating a dummy block to add capital costs associated with gas storage
    # this is a placeholder and should be updated with actual cost data
    m.fs.storage_tank.costing = Block()
    m.fs.storage_tank.purchase_cost = Var(
    initialize=100000,
    doc="Capital costs associated gas storage")
    m.fs.storage_tank.purchase_cost.fix()

    # No costing added for the interstage_cooler block as this is a dummy unit

    # Computing annualized capital cost
    def charge_capital_cost_rule(b):
        return m.fs.charge_capital_cost == ((
            m.fs.air_compressor.costing.purchase_cost +
            m.fs.storage_tank.purchase_cost) / m.number_of_years
        )
    m.fs.charge_capital_cost_eq = Constraint(rule=charge_capital_cost_rule)


def build_caes_charge(m):

    # build the storage model, set model inputs, and add bounds
    build_model(m)
    set_inputs(m)
    add_bounds(m)

    # check if the storage model is complete by asserting DOF = 0
    assert degrees_of_freedom(m) == 0

    # initialize the charge storage model
    initialize(m)

    # using storage cavern volume as the variable instead of L and D
    m.fs.storage_tank.volume_cons.deactivate()
    m.fs.storage_tank.control_volume.volume[:].fix(300000)

    return m


def integrate_storage(m):

    # using storage cavern volume as the variable instead of L and D
    m.fs.storage_tank.volume_cons.deactivate()
    m.fs.storage_tank.control_volume.volume[:].fix(300000)

    # Deactivating and writing new constraints to integrate storage
    m.fs.steam_cycle_heat_constraint.deactivate()
    m.fs.net_power_constraint.deactivate()

    @m.fs.Constraint(m.fs.time)
    def steam_cycle_heat_constraint_mod(fs, t):
        return (
            fs.steam_cycle_heat[t]
            == fs.HRSG_heat_duty[t]
            - fs.ASU_HP_steam_heat[t]
            + fs.interstage_cooler.heat_duty[0] * 1e-3
        )

    @m.fs.Constraint(m.fs.time)
    def net_power_constraint_mod(fs, t):
        return (
            fs.net_power[t]
            == fs.gross_power[t]
            - fs.auxiliary_load[t]
            - m.fs.air_compressor.control_volume.work[t] * 1e-3
            - fs.transformer_losses[t]
        )


def main():
    # # create charge model
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # build SOFC model
    m = SOFC.get_model()

    # build storage model
    m = build_caes_charge(m)

    # initialize unit models in storage
    integrate_storage(m)

    # add costing
    build_costing(m)

    return m


if __name__ == "__main__":

    optarg = {
        "max_iter": 300,
        "halt_on_ampl_error": "yes",
        "bound_push": 1e-23,
        "mu_init": 1e-6,
    }

    solver = get_solver("ipopt", optarg)

    m = main()

    solver.solve(m, tee=True, symbolic_solver_labels=True)

    print("              ")
    print("Compressor Duty", value(
        m.fs.air_compressor.control_volume.work[0]) * 1e-3)
    print("Cooler Duty", value(
        m.fs.interstage_cooler.heat_duty[0]) * 1e-3)
    print("Net Power", value(m.fs.net_power[0]))
