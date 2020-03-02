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
"""
Make sure the supercritical steam cycle example solves.
"""

__author__ = "Miguel Zamarripa"

import pytest
import pyomo.environ as pyo
from pyomo.network import Arc
import idaes.power_generation.flowsheets.supercritical_steam_cycle as steam_cycle
import idaes.power_generation.flowsheets.supercritical_power_plant.boiler_subflowsheet_build as blr
from idaes.power_generation.flowsheets.supercritical_power_plant.SCPC_full_plant import import_steam_cycle
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              activated_equalities_generator)
from idaes.generic_models.properties import iapws95
import argparse

solver_available = pyo.SolverFactory('ipopt').available()
prop_available = iapws95.iapws95_available()


@pytest.mark.slow
@pytest.mark.solver
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(not solver_available, reason="Solver not available")
def test_init():
    m, solver = blr.main()
    
    # initialize each unit at the time
    blr.initialize(m)
    blr.unfix_inlets(m)
    # check that the model solved properly and has 0 degrees of freedom
    assert(degrees_of_freedom(m)==0)


@pytest.mark.slow
@pytest.mark.solver
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(not solver_available, reason="Solver not available")
def test_boiler():
    m, solver = blr.main()
        
    # initialize each unit at the time
    blr.initialize(m)
   # unfix inlets to build arcs at the flowsheet level
    blr.unfix_inlets(m)
    m.fs.ATMP1.outlet.enth_mol[0].fix(62710.01)
    m.fs.ATMP1.SprayWater.flow_mol[0].unfix()
    result = solver.solve(m, tee=False)
    assert result.solver.termination_condition == \
            pyo.TerminationCondition.optimal
    assert result.solver.status == pyo.SolverStatus.ok
        
    assert pyo.value(m.fs.ECON.side_1.properties_out[0].temperature) == \
        pytest.approx(521.009,1)
#    assert gross_power_mw(m) == pytest.approx(620.8100259113626, abs=1e-3)


@pytest.mark.slow
@pytest.mark.solver
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(not solver_available, reason="Solver not available")
def test_power_plan():
    # import steam cycle and build concrete model
    m, solver = steam_cycle.main()
    print(degrees_of_freedom(m))
    #at this point we have a flowsheet with "steam cycle" that solves 
    # correctly, with 0 degrees of freedom.
    
    # next step is to import and build the boiler heat exchanger network
    # importing the boiler heat exchanger network from (boiler_subflowsheet_build.py)
    # will basically append all the unit models into our model ("m") 
    # model "m" has been created a few lines above
    
        # import the models (ECON, WW, PrSH, PlSH, FSH, Spliter, Mixer, Reheater)
        # see boiler_subflowhseet_build.py for a beter description
    blr.build_boiler(m.fs)
    #initialize boiler network models (one by one)
    blr.initialize(m)
    # at this point we have both flowsheets (steam cycle + boiler network)
    # in the same model/concrete object ("m")
    # however they are disconnected. Here we want to solve them at the same time
    # this is a square problem (i.e. degrees of freedom = 0)
#    print('solving square problem disconnected')
    results = solver.solve(m, tee=True)
    
    # at this point we want to connect the units in both flowsheets
    # Economizer inlet = Feed water heater 8 outlet (water)
    # HP inlet = Attemperator outlet (steam)
    # Reheater inlet (steam) = HP split 7 outlet (last stage of HP turbine)
    # IP inlet = Reheater outlet steam7
    blr.unfix_inlets(m)
    
    # deactivate constraints linking the FWH8 to HP turbine
    m.fs.boiler_pressure_drop.deactivate()
    m.fs.close_flow.deactivate()
    m.fs.turb.constraint_reheat_flow.deactivate()
    m.fs.turb.constraint_reheat_press.deactivate()
    m.fs.turb.constraint_reheat_temp.deactivate()
    m.fs.turb.inlet_split.inlet.enth_mol.unfix()
    m.fs.turb.inlet_split.inlet.pressure.unfix()
    
    m.fs.FHWtoECON = Arc(source = m.fs.fwh8.desuperheat.outlet_2,
                      destination = m.fs.ECON.side_1_inlet)
    
    m.fs.Att2HP = Arc(source = m.fs.ATMP1.outlet,
                   destination = m.fs.turb.inlet_split.inlet)
    
    m.fs.HPout2RH = Arc(source = m.fs.turb.hp_split[7].outlet_1,
                     destination = m.fs.RH.side_1_inlet)
    
    m.fs.RHtoIP = Arc(source = m.fs.RH.side_1_outlet,
                   destination =m.fs.turb.ip_stages[1].inlet)
    
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    
    #unfix boiler connections
    m.fs.ECON.side_1_inlet.flow_mol.unfix()
    m.fs.ECON.side_1_inlet.enth_mol[0].unfix()
    m.fs.ECON.side_1_inlet.pressure[0].unfix()
    m.fs.RH.side_1_inlet.flow_mol.unfix()
    m.fs.RH.side_1_inlet.enth_mol[0].unfix()
    m.fs.RH.side_1_inlet.pressure[0].unfix()
    m.fs.hotwell.makeup.flow_mol[:].setlb(-1.0)
    
#    if user has trouble with infeasible solutions, an easy test 
#    is to deactivate the link to HP turbine (m.fs.Att2HP_expanded "enth_mol and pressure" equalities) 
#    and fix inlet pressure and enth_mol to turbine (m.fs.turb.inlet_split.inlet)
    m.fs.turb.inlet_split.inlet.pressure.fix(2.423e7)
#   finally, since we want to maintain High Pressure (HP) inlet temperature constant (~866 K)
#   we need to fix Attemperator enthalpy outlet and unfix heat duty to Platen superheater
#   note fixing enthalpy to control temperature is only valid because pressure is also fixed
    m.fs.ATMP1.outlet.enth_mol[0].fix(62710.01)
    m.fs.PlSH.heat_duty[:].unfix()

    
#    print(degrees_of_freedom(m))
    solver.options = {
        "tol": 1e-6,
        "linear_solver": "ma27",
        "max_iter": 40,
    }
    #square problems tend to work better without bounds
    strip_bounds = pyo.TransformationFactory('contrib.strip_var_bounds')
    strip_bounds.apply_to(m, reversible=True)
    # this is the final solve with both flowsheets connected
    results = solver.solve(m, tee=True)
    strip_bounds.revert(m)