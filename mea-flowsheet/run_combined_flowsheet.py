# Python imports
import copy
import math
import logging
import sys

# Pyomo imports
import pyomo.environ as pyo
from pyomo.network import Arc, Port
from pyomo.common.config import ConfigValue, In, Bool, ListOf
from pyomo.util.check_units import assert_units_consistent

# IDAES imports
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlockData, declare_process_block_class

from idaes.models.unit_models import Mixer, MomentumMixingType, Translator

from idaes.models_extra.column_models.solvent_reboiler import SolventReboiler
from idaes.models_extra.column_models.solvent_condenser import SolventCondenser
from idaes.models_extra.column_models.MEAsolvent_column import MEAColumn

from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale
import idaes.core.util.model_serializer as ms

from idaes.core.surrogate.alamopy import AlamoSurrogate

# External IDAES imports (workspace)
sys.path.insert(0, 'flowsheets')
from mea_properties import (
    MEALiquidParameterBlock,
    MEAVaporParameterBlock,
    scale_mea_liquid_params,
    scale_mea_vapor_params,
    switch_liquid_to_parmest_params
)
from combined_flowsheet import MEACombinedFlowsheet, MEACombinedFlowsheetData
from mea_absorber_reformulated import MEAAbsorberFlowsheet, MEAAbsorberFlowsheetData
from mea_stripper_reformulated import MEAStripperFlowsheet, MEAStripperFlowsheetData
from run_absorber import set_absorber_linking_conditions
from run_stripper import set_stripper_linking_conditions

def set_stripper_linking_conditions(fs):
    # Fix conditions in flowsheet linking stream
    fs.reflux_mixer.rich_solvent.flow_mol.fix(24000)  # mol/sec
    fs.reflux_mixer.rich_solvent.temperature.fix(378) #380.15) # K
    fs.reflux_mixer.rich_solvent.pressure.fix(183700)  # Pa
    fs.reflux_mixer.rich_solvent.mole_frac_comp[0, "CO2"].fix(0.04)
    fs.reflux_mixer.rich_solvent.mole_frac_comp[0, "H2O"].fix(0.84)
    fs.reflux_mixer.rich_solvent.mole_frac_comp[0, "MEA"].fix(0.12)


def unfix_stripper_linking_conditions(fs):
    fs.reflux_mixer.rich_solvent.flow_mol.unfix()
    fs.reflux_mixer.rich_solvent.temperature.unfix()
    fs.reflux_mixer.rich_solvent.pressure.unfix()
    fs.reflux_mixer.rich_solvent.mole_frac_comp[...].unfix()



if __name__ == "__main__":
    logging.getLogger('pyomo.repn.plugins.nl_writer').setLevel(logging.ERROR)
    _log = idaeslog.getLogger("mea_flowsheet")
    # Create finite element set for stripper
    nfe = 40
    grid = [i / nfe for i in range(nfe + 1)]

    stripper_surrogate = AlamoSurrogate.load_from_file("alamo_surrogate_stripper.json")
    absorber_surrogate = AlamoSurrogate.load_from_file("alamo_surrogate_absorber.json")

    # Build flowsheet
    m = pyo.ConcreteModel()
    m.fs = MEACombinedFlowsheet(
        time_set=[0],
        absorber_finite_element_set=grid,
        stripper_finite_element_set=grid,
        # stripper_surrogate_enhancement_factor_model=stripper_surrogate,
        # absorber_surrogate_enhancement_factor_model=stripper_surrogate
    )
    # strip_statevar_bounds_for_stripper_initialization(m.stripper_section)
    switch_liquid_to_parmest_params(m.fs.stripper_section.liquid_properties, ions=True)
    switch_liquid_to_parmest_params(m.fs.stripper_section.liquid_properties_no_ions, ions=False)
    switch_liquid_to_parmest_params(m.fs.absorber_section.liquid_properties, ions=True)
    switch_liquid_to_parmest_params(m.fs.absorber_section.liquid_properties_no_ions, ions=False)

    # Set Initial guesses for inlets to sub-flowsheets
    m.fs.absorber_section.makeup_mixer.h2o_makeup.flow_mol.fix(1300)  # mol/sec
    iscale.calculate_scaling_factors(m.fs)
    m.fs.initialize_build(
        outlvl=idaeslog.DEBUG,
        optarg={
            # 'bound_push' : 1e-22,
            'nlp_scaling_method': 'user-scaling',
            'linear_solver': 'ma57',
            'OF_ma57_automatic_scaling': 'yes',
            'max_iter': 300,
            'tol': 1e-8,
            'halt_on_ampl_error': 'no',
            # 'mu_strategy': 'monotone',
        }
    )

    # Solve flowsheet
    solver = get_solver()
    results = solver.solve(m, tee=True)

    print("\nProcess model solved, starting to build costing equations...\n")
    MEACombinedFlowsheetData.add_costing(
        m.fs,
        export_economic_results=True,
        overwrite_economic_results=True
        )  # call costing module
    assert_units_consistent(m)

    # Print absorber, stripper results
    print()
    print("Absorber Column Results:")
    MEAAbsorberFlowsheetData.print_column_design_parameters(m.fs.absorber_section)
    print()
    print("Stripper Column Results:")
    MEAStripperFlowsheetData.print_column_design_parameters(m.fs.stripper_section)
    print()

    print("\nClean Gas to Stack")
    m.fs.absorber_section.clean_gas.display()
    print("\nCO2 Stream")
    m.fs.stripper_section.distillate.display()
    print("\nLean Solvent leaving Stripper Section")
    m.fs.stripper_section.lean_solvent.display()
    print("\nRich Solvent entering Stripper Section")
    m.fs.stripper_section.rich_solvent.display()
    print("\nH2O Makeup")
    m.fs.absorber_section.makeup.display()
    print("\nCondenser Duty")
    m.fs.stripper_section.condenser.heat_duty.display()
    print("\nReboiler Duty")
    m.fs.stripper_section.reboiler.heat_duty.display()
    print("\nRich Solvent Pump Pressure Change")
    m.fs.absorber_section.rich_solvent_pump.deltaP.display()
    print("\nLean Solvent Pump Pressure Change")
    m.fs.absorber_section.lean_solvent_pump.deltaP.display()
    print("\nAbsorber Liquid Inlet Pressure")
    m.fs.absorber_section.absorber.liquid_inlet.pressure.display()
    print("\nLean Solvent Loading")
    m.fs.lean_loading.display()
    print("\nRich Solvent Loading")
    m.fs.rich_loading.display()

    print("\nNGCC With Cap Costing Report")
    m.fs.ngcccap.report()
    print("\nNGCC No Cap Costing Report")
    m.fs.ngccnocap.report()
    print("\nLCOE Results")
    print()
    
    LCOE_results_list = ["capital_lcoe", "capital_lcoe_nocap", "fixed_lcoe",
                         "fixed_lcoe_nocap", "nonfuel_variable_lcoe",
                         "nonfuel_variable_lcoe_nocap", "fuel_lcoe",
                         "fuel_lcoe_nocap", "variable_lcoe", "variable_lcoe_nocap",
                         "transport_lcoe", "transport_lcoe_nocap", "LCOE",
                         "LCOE_nocap", "Cost_of_capture", "Cost_CO2_avoided",
                         "RemovalSystem_TPC", "RemovalSystem_Equip_Adjust",
                         "CC5_4", "FG_Cleanup_TPC", "NGCC_TPC", "plant_net_power",
                         "absorber_column_cost", "absorber_packing_cost",
                         "stripper_column_cost", "stripper_packing_cost",
                         "stripper_reboiler_cost", "stripper_condenser_cost",
                         "lean_rich_hex_cost", "FG_blower_cost", "DCC_column_cost",
                         "DCC_packing_cost", "DCC_pump_cost", "DCC_cooler_cost",
                         # "absorber_ic1_cost", "absorber_ic2_cost", "absorber_ic1pump_cost",
                         #"absorber_ic2pump_cost",
                         "ww_pump_cost", "ww_cooler_cost",
                         "rich_solvent_pump_cost", "stripper_reflux_drum_cost",
                         "lean_solvent_pump_cost",
                         # "lean_solvent_cooler_cost",
                         "solvent_stripper_reclaimer_cost", "solvent_filtration_cost",
                         "solvent_storage_tank_cost"]
    for result in LCOE_results_list:
        obj = getattr(m.fs.costing, result)
        print(result, ": ", pyo.value(obj), " ", pyo.units.get_units(obj))

    def print_units(variable):
        print(pyo.units.get_units(getattr(m.fs.costing, variable)))
        return
