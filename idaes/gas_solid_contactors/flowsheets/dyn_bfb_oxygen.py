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
An example of how a user might want to set up and solve (simulate) a
BFB CLC oxidation flowsheet

Author: Robert Parker and Chinedu Okoli
"""
import time
import sys
import os

from pyomo.environ import (ConcreteModel, SolverFactory, value,
        TransformationFactory, Var)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.dae import DerivativeVar
from pyomo.kernel import ComponentSet

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

from idaes.core import FlowsheetBlock, EnergyBalanceType
from idaes.core.util.model_statistics import degrees_of_freedom
from unit_models.bubbling_fluidized_bed import BubblingFluidizedBed
from property_packages.oxygen_iron_OC_oxidation.gas_phase_thermo \
        import Gas_Phase_Thermo_ParameterBlock
from property_packages.oxygen_iron_OC_oxidation.solid_phase_thermo \
        import Solid_Phase_Thermo_ParameterBlock
from property_packages.oxygen_iron_OC_oxidation.hetero_reactions \
        import HeteroReactionParameterBlock

from idaes.core.util.dyn_utils import (copy_values_at_time, 
        copy_non_time_indexed_values)
from idaes.core.util.initialization import initialize_by_time_element
import idaes.logger as idaeslog

from util.dyn_utils import (fix_initial_conditions,
        write_violated_equalities, remove_bounds_from,
        get_differential_equations)

def main():

    # Create solver
    solver_options = {'tol': 1e-5,
                      'mu_init': 1e-8,
                      'linear_solver': 'ma27',
                      'bound_push': 1e-8,
                      'halt_on_ampl_error': 'yes'}
    solver = SolverFactory('ipopt')
    solver.options = solver_options

    # Create spacial grid
    fe_set = [0]
    nxfe = 5
    nxcp = 3
    fe_a = 0.25
    fe_b = 0.2
    for i in range(1, nxfe+1):
        if i < nxfe*fe_a:
            fe_set.append(i*fe_b/(nxfe*fe_a))
        elif i == nxfe:
            fe_set.append(1)
        else:
            fe_set.append(fe_b + (i-nxfe*fe_a)*(1-fe_b)/(nxfe*(1-fe_a)))

    # Create time grid (5 second sampling interval)
    horizon = 60
    time_set = [0, horizon]
    ntfe = 12
    ntcp = 1

    # Create top-level and flowsheet models
    m = ConcreteModel(name='Full-horizon dynamic BFB-CLC reducer flowsheet model')
    m.fs = FlowsheetBlock(default={"dynamic": True,
                                   "time_set": time_set})

    # Set up thermo props and reaction props
    m.fs.gas_properties = Gas_Phase_Thermo_ParameterBlock()
    m.fs.solid_properties = Solid_Phase_Thermo_ParameterBlock()
    m.fs.hetero_reactions = HeteroReactionParameterBlock(
            default={"solid_property_package": m.fs.solid_properties,
                     "gas_property_package": m.fs.gas_properties})

    # Create MovingBed unit model
    m.fs.BFB = BubblingFluidizedBed(default={"flow_type": "co_current",
                     "finite_elements": nxfe,
                     "has_holdup": True,
                     "length_domain_set": fe_set,
                     "transformation_method": "dae.collocation",
                     "collocation_points": nxcp,
                     "transformation_scheme": "LAGRANGE-RADAU",
                     "bubble_config":
                     {"property_package": m.fs.gas_properties},
                     "gas_emulsion_config":
                     {"property_package": m.fs.gas_properties,
                      "has_pressure_change": True},
                     "solid_emulsion_config":
                     {"property_package": m.fs.solid_properties,
                      "reaction_package": m.fs.hetero_reactions
                      }})

    # Time-discretization
    discretizer = TransformationFactory('dae.collocation')
    discretizer.apply_to(m, wrt=m.fs.time, nfe=ntfe, ncp=ntcp,
            scheme='LAGRANGE-RADAU')
    #discretizer = TransformationFactory('dae.finite_difference')
    #discretizer.apply_to(m, wrt=m.fs.time, nfe=ntfe, scheme='BACKWARD')

    # Fix geometry variables
    m.fs.BFB.number_orifice.fix(2500) # [-]
    m.fs.BFB.bed_diameter.fix(6.5)    # [m]
    m.fs.BFB.bed_height.fix(5)        # [m]

    # Fix inlet variables for all t
    for t in m.fs.time:
        m.fs.BFB.gas_inlet.flow_mol[t].fix(1567.79)
        m.fs.BFB.solid_inlet.flow_mass[t].fix(1230.865)
        m.fs.BFB.gas_inlet.pressure[t].fix(1.86)

        m.fs.BFB.gas_inlet.temperature[t].fix(373)
        m.fs.BFB.gas_inlet.mole_frac[t, "O2"].fix(0.2095)
        m.fs.BFB.gas_inlet.mole_frac[t, "N2"].fix(0.7808)
        m.fs.BFB.gas_inlet.mole_frac[t, "CO2"].fix(0.0004)
        m.fs.BFB.gas_inlet.mole_frac[t, "H2O"].fix(0.0093)

        m.fs.BFB.solid_inlet.temperature[t].fix(1173.9)
        m.fs.BFB.solid_inlet.mass_frac[t, "Fe2O3"].fix(0.244)
        m.fs.BFB.solid_inlet.mass_frac[t, "Fe3O4"].fix(0.202)
        m.fs.BFB.solid_inlet.mass_frac[t, "Al2O3"].fix(0.554)
        '''
        For some reason, applying non-trivial inputs at t = 0 will
        cause this model to fail in a solve for consistent initial
        conditions.
        Needs further investigation.
        '''
        if t != m.fs.time.first():
            m.fs.BFB.gas_inlet.mole_frac[t, "O2"].fix(0.0595)
            m.fs.BFB.gas_inlet.mole_frac[t, "N2"].fix(0.8808)
            m.fs.BFB.gas_inlet.mole_frac[t, "CO2"].fix(0.0004)
            m.fs.BFB.gas_inlet.mole_frac[t, "H2O"].fix(0.0593)
#            m.fs.BFB.gas_in_flow[t].fix(372.81)

    for t in m.fs.time:
        for x in m.fs.BFB.length_domain:
            m.fs.BFB.solid_emulsion_region.properties[t, x]. \
                    _params.particle_dia = 1.5e-3

    t_start = time.time()

    # # #
    # Fix initial conditions
    # Options:
    # a) Fix num_diff_var variables to some desired values,
    #    then solve for rest of variables.
    #    These num_diff_var variables must be be sufficient to
    #    specify the differential variables (holdups).
    #    E.g. initial holdups could be fixed, or initial x/y/Tg/Ts
    #    could be fixed (as long as initial x/y are self-consistent).
    #    This includes fixing time derivatives to zero.
    #
    # b) Fix all variables to some consistent initial conditions.
    #    (Then deactivate model blocks/constraints at t = 0.)
    #
    # c) Specify a steady state to use as initial conditions.
    #
    # Here I use option c
    m_ss = ConcreteModel(name='Steady state BFB CLC reducer flowsheet model')
    m_ss.fs = FlowsheetBlock(default={"dynamic": False})
    m_ss.fs.gas_properties = Gas_Phase_Thermo_ParameterBlock()
    m_ss.fs.solid_properties = Solid_Phase_Thermo_ParameterBlock()
    m_ss.fs.hetero_reactions = HeteroReactionParameterBlock(
            default={"solid_property_package": m_ss.fs.solid_properties,
                     "gas_property_package": m_ss.fs.gas_properties})
    m_ss.fs.BFB = BubblingFluidizedBed(default={"finite_elements": nxfe,
                     "has_holdup": True,
                     "length_domain_set": fe_set,
                     "transformation_method": "dae.collocation",
                     "collocation_points": nxcp,
                     "transformation_scheme": "LAGRANGE-RADAU",
                     "bubble_config": {"property_package": m_ss.fs.gas_properties},
                     "gas_emulsion_config":
                     {"property_package": m_ss.fs.gas_properties,
                      "has_pressure_change": True},
                      "solid_emulsion_config":
                     {"property_package": m_ss.fs.solid_properties,
                      "reaction_package": m_ss.fs.hetero_reactions
                      }})

    # Fix inputs to steady state
    m_ss.fs.BFB.number_orifice.fix(2500) # [-]
    m_ss.fs.BFB.bed_diameter.fix(6.5)      # [m]
    m_ss.fs.BFB.bed_height.fix(5)     # [m]

    m_ss.fs.BFB.gas_inlet.flow_mol.fix(1567.79)
    m_ss.fs.BFB.gas_inlet.temperature.fix(373)
    m_ss.fs.BFB.gas_inlet.pressure.fix(1.86)
    m_ss.fs.BFB.gas_inlet.mole_frac[0, "O2"].fix(0.2095)
    m_ss.fs.BFB.gas_inlet.mole_frac[0, "N2"].fix(0.7808)
    m_ss.fs.BFB.gas_inlet.mole_frac[0, "CO2"].fix(0.0004)
    m_ss.fs.BFB.gas_inlet.mole_frac[0, "H2O"].fix(0.0093)

    m_ss.fs.BFB.solid_inlet.flow_mass.fix(1230.865)
    m_ss.fs.BFB.solid_inlet.temperature.fix(1173.9)
    m_ss.fs.BFB.solid_inlet.mass_frac[0, "Fe2O3"].fix(0.244)
    m_ss.fs.BFB.solid_inlet.mass_frac[0, "Fe3O4"].fix(0.202)
    m_ss.fs.BFB.solid_inlet.mass_frac[0, "Al2O3"].fix(0.554)

    # Update particle diameter:
    for t in m_ss.fs.time:
        for x in m_ss.fs.BFB.length_domain:
            m_ss.fs.BFB.solid_emulsion_region.properties[t, x]. \
                    _params.particle_dia = 1.5e-3


    # State arguments for initializing property state blocks
    blk = m_ss.fs.BFB
    gas_phase_state_args = {'flow_mol': blk.gas_inlet.flow_mol[0].value,
                            'temperature': blk.solid_inlet.temperature[0].value,
                            'pressure': blk.gas_inlet.pressure[0].value,
                            'mole_frac': {
                                'O2': blk.gas_inlet.mole_frac[0, 'O2'].value,
                                'N2': blk.gas_inlet.mole_frac[0, 'N2'].value,
                                'CO2': blk.gas_inlet.mole_frac[0, 'CO2'].value,
                                'H2O': blk.gas_inlet.mole_frac[0, 'H2O'].value}
                            }
    solid_phase_state_args = {'flow_mass': blk.solid_inlet.flow_mass[0].value,
                              'temperature': blk.solid_inlet.temperature[0].value,
                              'mass_frac': {
                              'Fe2O3': blk.solid_inlet.mass_frac
                                          [0, 'Fe2O3'].value,
                              'Fe3O4': blk.solid_inlet.mass_frac
                                          [0, 'Fe3O4'].value,
                              'Al2O3': blk.solid_inlet.mass_frac
                                          [0, 'Al2O3'].value}
                              }

    m_ss.fs.BFB.initialize(outlvl=idaeslog.WARNING,
                           gas_phase_state_args=gas_phase_state_args,
                           solid_phase_state_args=solid_phase_state_args)
    results = solver.solve(m_ss)

    if debug:
        write_violated_equalities(m_ss)

    # Initialize dynamic model (except fixed inputs) to this steady state
    for t in m.fs.time:
        copy_values_at_time(m.fs, m_ss.fs, t, 0.0, copy_fixed=False,
                outlvl=idaeslog.ERROR)

    copy_non_time_indexed_values(m.fs, m_ss.fs)

    # Fix initial conditions to current values (from above steady state)
    fix_initial_conditions(m.fs, m.fs.time)

    # Correct for input overlap, if necessary due to choice of inputs:
    for j in m.fs.gas_properties.component_list:
        m.fs.BFB.bubble_region.material_holdup[0, 0, 'Vap', j].unfix()
        m.fs.BFB.gas_emulsion_region.material_holdup[0, 0, 'Vap', j].unfix()
    for j in m.fs.solid_properties.component_list:
        m.fs.BFB.solid_emulsion_region.material_holdup[0, 0, 'Sol', j].unfix()

    m.fs.BFB.bubble_region.energy_holdup[0, 0, 'Vap'].unfix()
    m.fs.BFB.gas_emulsion_region.energy_holdup[0, 0, 'Vap'].unfix()
    m.fs.BFB.solid_emulsion_region.energy_holdup[0, 0, 'Sol'].unfix()

    # # #

    assert degrees_of_freedom(m.fs) == 0

    remove_bounds_from(m.fs)

    for var in m.fs.component_objects(Var):
        if isinstance(var, DerivativeVar):
            if m.fs.time in ComponentSet(
                    var.get_continuousset_list()):
                for index in var:
                    var[index].set_value(0)

    if debug:
        write_violated_equalities(m)

    initialize_by_time_element(m.fs, m.fs.time, solver=solver,
            outlvl=idaeslog.INFO)
    assert degrees_of_freedom(m.fs) == 0

    results = solver.solve(m.fs, tee=True)
    if debug:
        with open('dyn_sol.txt', 'w') as f:
            m.display(ostream=f)

    print('Outlet conversions:')
    print('      Solid phase')
    for t in m.fs.time:
        print(f't = {t}:',
                m.fs.BFB.solid_emulsion_region.reactions[t, 1].OC_conv.value)

#%%        # ---------------------------------------------------------------------

if __name__ == "__main__":
    debug = False
    main()

