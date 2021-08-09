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
moving bed CLC flowsheet

Author: Robert Parker and Chinedu Okoli
"""
import time
import sys
import os

from pyomo.environ import (ConcreteModel, SolverFactory, value,
        TransformationFactory, Constraint, Var)
from pyomo.dae import DerivativeVar
from pyomo.kernel import ComponentSet
from pyomo.util.calc_var_value import calculate_variable_from_constraint

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

from idaes.core import FlowsheetBlock, EnergyBalanceType
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog

from unit_models.moving_bed import MovingBed
from property_packages.methane_iron_OC_reduction.gas_phase_thermo \
        import Gas_Phase_Thermo_ParameterBlock
from property_packages.methane_iron_OC_reduction.solid_phase_thermo \
        import Solid_Phase_Thermo_ParameterBlock
from property_packages.methane_iron_OC_reduction.hetero_reactions \
        import HeteroReactionParameterBlock

from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util.dyn_utils import (copy_values_at_time, 
        copy_non_time_indexed_values)

from util.dyn_utils import (fix_initial_conditions,
        write_violated_equalities, remove_bounds_from,
        get_differential_equations)

def main():

    # Create solver
    solver_options = {'tol': 1e-5,
                      'mu_init': 1e-8,
                      'linear_solver': 'ma27',
                      'halt_on_ampl_error': 'yes',
                      'bound_push': 1e-8}
    solver = SolverFactory('ipopt')
    solver.options = solver_options

    # Create spacial grid
    fe_set = [0]
    nxfe = 10
    nxcp = 3

    fe_set = [0, 0.004]
    fe_a = 1/4.0
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
    m = ConcreteModel(name='Full-horizon dynamic MBCLC flowsheet model')
    m.fs = FlowsheetBlock(default={"dynamic": True,
                                   "time_set": time_set})

    # Set up thermo props and reaction props
    m.fs.gas_properties = Gas_Phase_Thermo_ParameterBlock()
    m.fs.solid_properties = Solid_Phase_Thermo_ParameterBlock()
    m.fs.hetero_reactions = HeteroReactionParameterBlock(
#            default={"property_package": m.fs.solid_properties})
            default={"solid_property_package": m.fs.solid_properties,
                     "gas_property_package": m.fs.gas_properties})

    # Create MovingBed unit model
    m.fs.MB = MovingBed(default={"finite_elements": nxfe,
                                 "has_holdup": True,
                                 "length_domain_set": fe_set,
                                 "transformation_method": "dae.collocation",
                                 "collocation_points": nxcp,
                                 "transformation_scheme": "LAGRANGE-RADAU",
                                 "gas_phase_config":
                                    {"property_package": m.fs.gas_properties,
                                     "has_pressure_change": True,
                                     "pressure_drop_type": "ergun_correlation"},
                                 "solid_phase_config":
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
    m.fs.MB.bed_diameter.fix(6.5)
    m.fs.MB.bed_height.fix(5)

    # Fix inlet variables for all t
    for t in m.fs.time:
        m.fs.MB.gas_inlet.flow_mol[t].fix(128.20513)
        m.fs.MB.gas_inlet.pressure[t].fix(2.00)

        m.fs.MB.solid_inlet.flow_mass[t].fix(591.4)

        if t != m.fs.time.first():
        # "Differential variables" are already specified by the initial 
        # conditions at time.first(). Alternatively, initial condition
        # specification could be omitted at the inlets.
            m.fs.MB.gas_inlet.temperature[t].fix(298.15)
            m.fs.MB.gas_inlet.mole_frac[t, "CO2"].fix(0.5)
            m.fs.MB.gas_inlet.mole_frac[t, "H2O"].fix(0.0)
            m.fs.MB.gas_inlet.mole_frac[t, "CH4"].fix(0.5)

            m.fs.MB.solid_inlet.temperature[t].fix(1183.15)
            m.fs.MB.solid_inlet.mass_frac[t, "Fe2O3"].fix(0.45)
            m.fs.MB.solid_inlet.mass_frac[t, "Fe3O4"].fix(1e-9)
            m.fs.MB.solid_inlet.mass_frac[t, "Al2O3"].fix(0.55)

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
    # b) Fix all variables to some consistent initial conditions.
    #    In this case the model equations at time.first() must be
    #    deactivated so the model is not overspecified.
    # c) Specify a steady state to use as initial conditions.
    #
    # Here I use option c
    m_ss = ConcreteModel(name='Full-horizon dynamic MBCLC flowsheet model')
    m_ss.fs = FlowsheetBlock(default={"dynamic": False})
    m_ss.fs.gas_properties = Gas_Phase_Thermo_ParameterBlock()
    m_ss.fs.solid_properties = Solid_Phase_Thermo_ParameterBlock()
    m_ss.fs.hetero_reactions = HeteroReactionParameterBlock(
            default={"solid_property_package": m_ss.fs.solid_properties,
                     "gas_property_package": m_ss.fs.gas_properties})
    m_ss.fs.MB = MovingBed(default={"finite_elements": nxfe,
                                 "has_holdup": True,
                                 "length_domain_set": fe_set,
                                 "transformation_method": "dae.collocation",
                                 "collocation_points": nxcp,
                                 "transformation_scheme": "LAGRANGE-RADAU",
                                 "gas_phase_config":
                                    {"property_package": m_ss.fs.gas_properties,
                                     "has_pressure_change": True,
                                     "pressure_drop_type": "ergun_correlation"},
                                 "solid_phase_config":
                                    {"property_package": m_ss.fs.solid_properties,
                                     "reaction_package": m_ss.fs.hetero_reactions
                                     }})

    # Fix inputs to desired initial steady state conditions
    m_ss.fs.MB.bed_diameter.fix(6.5)
    m_ss.fs.MB.bed_height.fix(5)

    m_ss.fs.MB.gas_inlet.flow_mol[0].fix(128.20513)
    m_ss.fs.MB.gas_inlet.temperature[0].fix(298.15)
    m_ss.fs.MB.gas_inlet.pressure[0].fix(2.00)
    m_ss.fs.MB.gas_inlet.mole_frac[0, "CO2"].fix(0.02499)
    m_ss.fs.MB.gas_inlet.mole_frac[0, "H2O"].fix(0.00001)
    m_ss.fs.MB.gas_inlet.mole_frac[0, "CH4"].fix(0.97500)

    m_ss.fs.MB.solid_inlet.flow_mass[0].fix(591.4)
    m_ss.fs.MB.solid_inlet.temperature[0].fix(1183.15)
    m_ss.fs.MB.solid_inlet.mass_frac[0, "Fe2O3"].fix(0.45)
    m_ss.fs.MB.solid_inlet.mass_frac[0, "Fe3O4"].fix(1e-9)
    m_ss.fs.MB.solid_inlet.mass_frac[0, "Al2O3"].fix(0.55)

    blk = m_ss.fs.MB
    gas_phase_state_args = {
            'flow_mol': blk.gas_inlet.flow_mol[0].value,
            'temperature': blk.solid_inlet.temperature[0].value,
            'pressure': blk.gas_inlet.pressure[0].value,
            'mole_frac': {
                'CH4': blk.gas_inlet.mole_frac[0, 'CH4'].value,
                'CO2': blk.gas_inlet.mole_frac[0, 'CO2'].value,
                'H2O': blk.gas_inlet.mole_frac[0, 'H2O'].value}}
    solid_phase_state_args = {
            'flow_mass': blk.solid_inlet.flow_mass[0].value,
            'temperature': blk.solid_inlet.temperature[0].value,
            'mass_frac': {
                'Fe2O3': blk.solid_inlet.mass_frac[0, 'Fe2O3'].value,
                'Fe3O4': blk.solid_inlet.mass_frac[0, 'Fe3O4'].value,
                'Al2O3': blk.solid_inlet.mass_frac[0, 'Al2O3'].value}}

    # Initialize and solve steady state problem
    m_ss.fs.MB.initialize(outlvl=idaeslog.WARNING,
                          gas_phase_state_args=gas_phase_state_args,
                          solid_phase_state_args=solid_phase_state_args,
                          optarg={
                                  'tol': 1e-6,
                                  'mu_init': 1e-8,
                                  'linear_solver': 'ma27',
                                  'halt_on_ampl_error': 'yes',
                                  'bound_push': 1e-8})
    results = solver.solve(m_ss)

    # Initialize dynamic model (except fixed inputs) to this steady state
    # Initial conditions have not yet been fixed, so they will be copied
    # over (as desired).
    for t in m.fs.time:
        copy_values_at_time(m.fs, m_ss.fs, t, 0.0, copy_fixed=False,
                outlvl=idaeslog.ERROR)

    # Non-time indexed variables are things like design variables, uncertain
    # parameters, variables explicitly tied to some point in time (initial
    # or terminal, usually), or variables that store an average over time.
    # 
    # For simulation, these variables should usually be fixed or inactive,
    # but I copy their values anyway to try and ensure that all variables
    # are initialized.
    copy_non_time_indexed_values(m.fs, m_ss.fs)

    # Fix initial conditions to current values (from above steady state)
    # This function is pretty crude. It just fixes variables that have time 
    # derivatives, regardless of whether those variables are ~actually~ 
    # differential.
    fix_initial_conditions(m.fs, m.fs.time)

    assert degrees_of_freedom(m.fs) == 0

    # I don't think we want to be performing a dynamic simulation with bounds.
    # This does more than Pyomo's strip_bounds transformation. It also changes
    # domains to Reals.
    remove_bounds_from(m.fs)

    # Derivatives are not present in the steady state model, so they will
    # not necessarily have values here. I could just assign them values
    # of zero, but instead I calculate them from their differential equations.
    for var in m.fs.component_objects(Var):
        if isinstance(var, DerivativeVar):
            if m.fs.time in ComponentSet(
                    var.get_continuousset_list()):
                for index in var:
                    var[index].set_value(0)

    # For sanity, make sure the only equations that are violated
    # are those I expect to be violated. (Those due to perturbed
    # inputs and possibly initial discretization equations.)
    if debug: 
        write_violated_equalities(m)

    # Initialize by time element
    initialize_by_time_element(m.fs, m.fs.time, solver=solver,
            outlvl=idaeslog.INFO)
    print('dof after integrating:', degrees_of_freedom(m.fs))

    # And solve. This will be trivial here for the square problem.
    results = solver.solve(m.fs, tee=True)
    if debug:
        with open('dyn_sol.txt', 'w') as f:
            m.display(ostream=f)

    # Display some very basic information
    print('Outlet conversions:')
    print('      Solid phase')
    for t in m.fs.time:
        print(f't = {t}:', m.fs.MB.solid_phase.reactions[t, 0].OC_conv.value)

#%%        # ---------------------------------------------------------------------

if __name__ == "__main__":
    debug = False
    main()
