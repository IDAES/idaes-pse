##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

import pyomo.environ as pe
import pyomo.network as network
from pyomo.network.plugins import expand_arcs
from idaes.core import FlowsheetBlock
from idaes.unit_models import Feed, Heater, Product, Flash, StoichiometricReactor, Mixer, Separator
from idaes.unit_models.mixer import MomentumMixingType
from methanol_reaction import MethanolReactionParameterBlock
from compressor import IdealGasIsentropicCompressor
from expander import Expander
from methanol_param_VLE import PhysicalParameterBlock
from pyomo.common.config import ConfigBlock
from idaes.core.util.misc import add_object_reference
import math
import sys


def create_model():
    m = pe.ConcreteModel()
    m.fs = fs = FlowsheetBlock(default={"dynamic": False})
    fs.vapor_props = vapor_props = PhysicalParameterBlock(default={"valid_phase": 'Vap'})
    fs.properties = props = PhysicalParameterBlock(default={"valid_phase": ('Vap', 'Liq')})
    fs.reaction_params = MethanolReactionParameterBlock(default={'property_package': vapor_props})

    fs.feed = feed = Feed(default={"property_package": vapor_props})
    fs.compressor1 = IdealGasIsentropicCompressor(default={"property_package": vapor_props,
                                                           "has_phase_equilibrium": False})
    fs.cooler1 = Heater(default={"property_package": vapor_props, "has_phase_equilibrium": False})
    fs.compressor2 = IdealGasIsentropicCompressor(default={"property_package": vapor_props,
                                                           "has_phase_equilibrium": False})
    fs.equal_electric = pe.Constraint(expr=fs.compressor1.work[0.0] == fs.compressor2.work[0.0])
    fs.mixer = Mixer(default={'property_package': vapor_props,
                              'inlet_list': ['feed', 'recycle'],
                              'momentum_mixing_type': MomentumMixingType.equality})

    # Reactor
    fs.reactor = StoichiometricReactor(default={'property_package': vapor_props,
                                                'reaction_package': fs.reaction_params,
                                                'has_heat_of_reaction': True,
                                                'has_pressure_change': False})
    fs.reactor.conversion_eq = pe.Var()
    fs.reactor.t_inv = pe.Var()
    fs.reactor.p_sq_inv = pe.Var()
    fs.reactor.conversion = pe.Var()
    fs.reactor.consumption_rate = pe.Var()
    fs.reactor.t_inv_con = pe.Constraint(expr=fs.reactor.t_inv * fs.reactor.outlet.temperature[0] == 1)
    fs.reactor.p_sq_inv_con = pe.Constraint(expr=fs.reactor.p_sq_inv * fs.reactor.inlet.pressure[0] ** 2 == 1)
    fs.reactor.conversion_eq_con = pe.Constraint(expr=(fs.reactor.conversion_eq == 0.415 * (1 - 26.25 * pe.exp(-18 * fs.reactor.t_inv) * fs.reactor.p_sq_inv)))
    fs.reactor.conversion_con = pe.Constraint(expr=(fs.reactor.conversion == fs.reactor.conversion_eq * (1 - pe.exp(-5)) * (fs.reactor.inlet.mole_frac[0, "H2"] + fs.reactor.inlet.mole_frac[0, "CO"] + fs.reactor.inlet.mole_frac[0, "CH3OH"])))
    fs.reactor.consumption_rate_con = pe.Constraint(expr=(fs.reactor.consumption_rate == fs.reactor.conversion * fs.reactor.inlet.mole_frac[0, "H2"] * fs.reactor.inlet.flow_mol[0]))
    fs.reactor.h2_consumption_con = pe.Constraint(expr=(fs.reactor.outlet.flow_mol[0] * fs.reactor.outlet.mole_frac[0, "H2"] == fs.reactor.inlet.flow_mol[0] * fs.reactor.inlet.mole_frac[0, "H2"] - fs.reactor.consumption_rate))

    fs.expander = Expander(default={'property_package': vapor_props,
                                    'has_phase_equilibrium': False})
    fs.cooler2 = Heater(default={"property_package": vapor_props, "has_phase_equilibrium": False})
    fs.flash = Flash(default={"property_package": props})
    fs.purge_splitter = Separator(default={'property_package': vapor_props,
                                           'outlet_list': ['purge', 'recycle'],
                                           'ideal_separation': False})
    fs.compressor3 = IdealGasIsentropicCompressor(default={"property_package": vapor_props,
                                                           "has_phase_equilibrium": False})


    ###########################
    # Set scaling factors
    ###########################
    fs.compressor1.control_volume.scaling_factor_energy.value = 1
    fs.compressor2.control_volume.scaling_factor_energy.value = 1
    fs.cooler1.control_volume.scaling_factor_energy.value = 1
    fs.flash.control_volume.scaling_factor_energy.value = 1
    fs.reactor.control_volume.scaling_factor_energy.value = 1
    fs.cooler2.control_volume.scaling_factor_energy.value = 1
    fs.compressor3.control_volume.scaling_factor_energy.value = 1
    fs.mixer.scaling_factor_energy.value = 1

    fs.cooler1.control_volume.scaling_factor_pressure.value = 1
    fs.flash.control_volume.scaling_factor_pressure.value = 1
    fs.reactor.control_volume.scaling_factor_pressure.value = 1
    fs.cooler2.control_volume.scaling_factor_pressure.value = 1

    ###########################
    # Objective
    ###########################
    m.objective = pe.Objective(expr=(-fs.flash.liq_outlet.flow_mol[0.0]))

    ###########################
    # Connect Units
    ###########################
    fs.stream1 = network.Arc(source=feed.outlet, destination=fs.compressor1.inlet)
    fs.stream2 = network.Arc(source=fs.compressor1.outlet, destination=fs.cooler1.inlet)
    fs.stream3 = network.Arc(source=fs.cooler1.outlet, destination=fs.compressor2.inlet)
    fs.stream4 = network.Arc(source=fs.compressor2.outlet, destination=fs.mixer.feed)
    fs.stream5 = network.Arc(source=fs.mixer.outlet, destination=fs.reactor.inlet)
    fs.stream6 = network.Arc(source=fs.reactor.outlet, destination=fs.expander.inlet)
    fs.stream7 = network.Arc(source=fs.expander.outlet, destination=fs.cooler2.inlet)
    fs.stream8 = network.Arc(source=fs.cooler2.outlet, destination=fs.flash.inlet)
    fs.stream9 = network.Arc(source=fs.flash.vap_outlet, destination=fs.purge_splitter.inlet)
    fs.stream10 = network.Arc(source=fs.purge_splitter.recycle, destination=fs.compressor3.inlet)
    fs.stream11 = network.Arc(source=fs.compressor3.outlet, destination=fs.mixer.recycle)
    pe.TransformationFactory("network.expand_arcs").apply_to(m)

    ###########################
    # Set problem specs
    ###########################
    feed.flow_mol.fix(3.40898)
    feed.pressure.fix(1)
    feed.temperature.fix(3)
    feed.mole_frac[0.0, "CH4"].fix(0.05)
    feed.mole_frac[0.0, "CO"].fix(0.3)
    feed.mole_frac[0.0, "H2"].fix(0.6)
    feed.mole_frac[0.0, "CH3OH"].fix(0.05)

    fs.cooler1.heat_duty[0.0].setub(0)  # it is a cooler
    fs.cooler1.outlet.temperature[0.0].setlb(3)

    fs.flash.heat_duty.fix(0)
    fs.flash.deltaP.fix(0)

    fs.cooler2.heat_duty[0.0].setub(0)  # it is a cooler
    fs.cooler2.outlet.temperature[0.0].setlb(3)

    fs.compressor2.outlet.pressure[0.0].setub(10)
    fs.flash.liq_outlet.mole_frac[0.0, 'CH3OH'].expr.setlb(0.9)

    fs.purge_splitter.split_fraction[0.0, 'purge'].fix(0.05)

    ###########################
    # Prepare for initialization
    ###########################
    fs.compressor1.outlet.pressure.fix(5)
    fs.compressor2.outlet.pressure.fix(10)
    fs.cooler1.outlet.temperature.fix(3)
    fs.expander.outlet.pressure[0.0].fix(5)
    fs.cooler2.outlet.temperature[0.0].fix(3)
    fs.flash.liq_outlet.mole_frac[0.0, 'CH3OH'].expr.setlb(None)

    # Setup decomposition process
    seq = network.SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 5

    # Determine tear stream and calculation order
    G = seq.create_graph(m)
    heu_result = seq.tear_set_arcs(G, method="heuristic")
    order = seq.calculation_order(G)

    # Display tear stream and calculation order
    print("Tear")
    for o in heu_result:
        print(o.name)
    print("Order")
    for o in order:
        for oo in o:
            print(oo.name)

    # Set guesses for tear stream
    tear_guesses = {
        "flow_mol": {0: 10},
        "mole_frac": {(0, 'CH3OH'): 0.06,
                      (0, 'CH4'): 0.21,
                      (0, 'CO'): 0.24,
                      (0, 'H2'): 0.50},
        "temperature": {0: 4.4},
        "pressure": {0: 15}}
    seq.set_guesses_for(m.fs.reactor.inlet, tear_guesses)

    # Define method for initialising each block
    def function(unit):
        unit.initialize(outlvl=1)

    # Run sequential initialisation
    seq.run(m, function)

    ###########################
    # Unfix vars that were fixed for initialization
    ###########################
    m.fs.compressor1.outlet.pressure.unfix()
    m.fs.compressor2.outlet.pressure.unfix()
    m.fs.cooler1.outlet.temperature.unfix()
    m.fs.expander.outlet.pressure.unfix()
    m.fs.cooler2.outlet.temperature.unfix()
    m.fs.flash.liq_outlet.mole_frac[0.0, 'CH3OH'].expr.setlb(0.9)

    return m


m = create_model()

opt = pe.SolverFactory('ipopt')
opt.options['linear_solver'] = 'mumps'
res = opt.solve(m, tee=True)

print('******************************')
print('Compressor 1:')
m.fs.compressor1.outlet.display()
print('\n\n')
print('******************************')
print('Cooler 1:')
m.fs.cooler1.outlet.display()
print('\n\n')
print('******************************')
print('Compressor 2:')
m.fs.compressor2.outlet.display()
print('\n\n')
print('******************************')
print('Mixer:')
m.fs.mixer.outlet.display()
print('\n\n')
print('******************************')
print('Reactor:')
m.fs.reactor.outlet.display()
print('\n\n')
print('******************************')
print('Expander:')
m.fs.expander.outlet.display()
print('\n\n')
print('******************************')
print('Cooler 2:')
m.fs.cooler2.outlet.display()
print('\n\n')
print('******************************')
print('Flash Liquid:')
m.fs.flash.liq_outlet.display()
print('\n\n')
print('******************************')
print('Flash Vapor:')
m.fs.flash.vap_outlet.display()
print('\n\n')
print('******************************')
print('Purge:')
m.fs.purge_splitter.purge.display()
print('\n\n')
print('******************************')
print('Recycle:')
m.fs.purge_splitter.recycle.display()
print('\n\n')
print('******************************')
print('Compressor 3:')
m.fs.compressor3.outlet.display()
