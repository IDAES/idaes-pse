# Python imports
import copy
import math
import logging
import sys

# Pyomo imports
import pyomo.environ as pyo
from pyomo.network import Arc, Port

# IDAES imports
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.models.unit_models import Pump, HeatExchangerNTU, Mixer
from idaes.models.unit_models.mixer import MomentumMixingType
from idaes.models.properties.modular_properties import GenericParameterBlock

from idaes.models_extra.column_models.properties.MEA_vapor import (
    flue_gas as vapor_config,
)
from idaes.models_extra.column_models.properties.MEA_solvent import (
    configuration as liquid_config,
)

from idaes.core.solvers.get_solver import get_solver
import idaes.logger as idaeslog
from idaes.core.util.initialization import propagate_state
import idaes.core.util.scaling as iscale
from idaes.core.surrogate.alamopy import AlamoSurrogate

sys.path.insert(0,'flowsheets')
from mea_properties import (
    MEALiquidParameterBlock,
    MEAVaporParameterBlock,
    scale_mea_liquid_params,
    scale_mea_vapor_params,
    switch_liquid_to_parmest_params,
)
from mea_absorber_reformulated import MEAAbsorberFlowsheet

def set_absorber_linking_conditions(fs):
    # Fix conditions in flowsheet linking stream
    # Lean solvent feed to sub-flowsheet
    fs.lean_solvent_pump.inlet.flow_mol.fix(23000)  # mol/sec
    fs.lean_solvent_pump.inlet.temperature.fix(394)  # K
    fs.lean_solvent_pump.inlet.pressure.fix(183700)  # Pa
    fs.lean_solvent_pump.inlet.mole_frac_comp[0, "CO2"].fix(0.02735)
    fs.lean_solvent_pump.inlet.mole_frac_comp[0, "H2O"].fix(0.85097)
    fs.lean_solvent_pump.inlet.mole_frac_comp[0, "MEA"].fix(0.12180)


def unfix_absorber_linking_conditions(fs):
    fs.lean_solvent_pump.inlet.flow_mol.unfix()
    fs.lean_solvent_pump.inlet.temperature.unfix()
    fs.lean_solvent_pump.inlet.pressure.unfix()
    fs.lean_solvent_pump.inlet.mole_frac_comp[...].unfix()



if __name__ == "__main__":
    logging.getLogger('pyomo.repn.plugins.nl_writer').setLevel(logging.ERROR)
    # Create finite element set for absorber
    nfe = 40
    grid = [i / nfe for i in range(nfe + 1)]

    # surrogate = AlamoSurrogate.load_from_file("alamo_surrogate_absorber.json")

    # Build flowsheet
    m = pyo.ConcreteModel()
    m.fs = MEAAbsorberFlowsheet(
        time_set=[0],
        column_finite_element_set=grid,
        # surrogate_enhancement_factor_model=surrogate,
    )

    switch_liquid_to_parmest_params(m.fs.liquid_properties, ions=True)
    switch_liquid_to_parmest_params(m.fs.liquid_properties_no_ions, ions=False)

    # Set some inputs for test case
    set_absorber_linking_conditions(m.fs)

    # Apply flowsheet scaling
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

    m.fs.absorber.vapor_outlet.display()
    m.fs.absorber.liquid_outlet.display()
    m.fs.lean_rich_heat_exchanger.report()
    m.fs.makeup_mixer.report()

    print("\n-------- Absorber Simulation Results --------")
    m.fs.print_column_design_parameters()

    # m.fs.absorber.liquid_phase.properties[:,:].mole_frac_comp.pprint()
    # m.fs.absorber.vapor_phase.properties[:,:].mole_frac_comp.pprint()
    # m.fs.absorber.vapor_phase.properties[:,:].pressure.pprint()
    # m.fs.absorber.liquid_phase.properties[:,:].temperature.pprint()
    # m.fs.absorber.mass_transfer_coeff_liq.pprint()
    # m.fs.absorber.mass_transfer_coeff_vap.pprint()
