# Python imports
import copy
import math
import logging
import sys

# Pyomo imports
import pyomo.environ as pyo
from pyomo.network import Arc, Port
from pyomo.common.config import ConfigValue, In, Bool, ListOf

# IDAES imports
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import FlowsheetBlockData, declare_process_block_class

from idaes.models.unit_models import Mixer, MomentumMixingType

from idaes.models_extra.column_models.solvent_reboiler import SolventReboiler
from idaes.models_extra.column_models.solvent_condenser import SolventCondenser
from idaes.models_extra.column_models.MEAsolvent_column import MEAColumn

from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale

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
from mea_stripper_reformulated import MEAStripperFlowsheet

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
    # Create finite element set for stripper
    nfe = 40
    grid = [i / nfe for i in range(nfe + 1)]

    surrogate = AlamoSurrogate.load_from_file("alamo_surrogate_stripper.json")

    # Build flowsheet
    m = pyo.ConcreteModel()
    m.fs = MEAStripperFlowsheet(
        time_set=[0],
        column_finite_element_set=grid,
        # surrogate_enhancement_factor_model=surrogate,
    )
    # strip_statevar_bounds_for_stripper_initialization(m.stripper_section)
    switch_liquid_to_parmest_params(m.fs.liquid_properties, ions=True)
    switch_liquid_to_parmest_params(m.fs.liquid_properties_no_ions, ions=False)

    # Set some inputs for test case
    set_stripper_linking_conditions(m.fs)

    # Apply flowsheet scaling
    iscale.calculate_scaling_factors(m.fs)

    # Initialize flowsheet
    # results = initialize_stripper_flowsheet(m.stripper_section,
    #                                         outlvl=idaeslog.DEBUG,  tee=True)
    m.fs.strip_statevar_bounds_for_stripper_initialization()
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

    # check_scaling(m.stripper_section)

    # m.stripper_section.stripper.vapor_inlet.display()
    # m.stripper_section.reboiler.vapor_reboil.display()

    # m.stripper_section.condenser.vapor_outlet.display()
    # m.stripper_section.reboiler.bottoms.display()
    # m.stripper_section.reboiler.heat_duty.display()

    # define_column_design_parameters(m.stripper_section)
    print("\n-------- stripper Simulation Results --------")
    m.fs.print_column_design_parameters()

    print("\nStripper - Condenser connection")
    print("------------------------------------------------------------------")
    m.fs.stripper.vapor_outlet.display()
    m.fs.condenser.inlet.display()
    print("------------------------------------------------------------------")
    print("\nCondenser - Reflux Mixer connection")
    print("------------------------------------------------------------------")
    m.fs.condenser.reflux.display()
    m.fs.reflux_mixer.reflux.display()
    print("\n Rich solvent stream (reflux mixer inlet):")
    m.fs.reflux_mixer.rich_solvent.display()
    print("------------------------------------------------------------------")
    print("\nReflux Mixer - Stripper connection")
    print("------------------------------------------------------------------")
    m.fs.reflux_mixer.outlet.display()
    m.fs.stripper.liquid_inlet.display()
    print("------------------------------------------------------------------")
    print("\nStripper - Reboiler LIQUID connection")
    print("------------------------------------------------------------------")
    m.fs.stripper.liquid_outlet.display()
    m.fs.reboiler.inlet.display()
    print("------------------------------------------------------------------")
    print("\nReboiler - Stripper VAPOR connection")
    print("------------------------------------------------------------------")
    m.fs.reboiler.vapor_reboil.display()
    m.fs.stripper.vapor_inlet.display()
    print("\nHot lean solvent stream (reboiler outlet):")
    m.fs.reboiler.bottoms.display()

    # m.fs.stripper.liquid_phase.properties[:,:].mole_frac_comp.pprint()
    # m.fs.stripper.vapor_phase.properties[:,:].mole_frac_comp.pprint()
    # m.fs.stripper.vapor_phase.properties[:,:].pressure.pprint()
    # m.fs.stripper.liquid_phase.properties[:,:].temperature.pprint()
    # m.fs.stripper.mass_transfer_coeff_liq.pprint()
    # m.fs.stripper.mass_transfer_coeff_vap.pprint()
    print("ok, boomer")