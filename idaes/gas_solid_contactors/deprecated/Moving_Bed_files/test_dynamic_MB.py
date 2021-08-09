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
Tests for ControlVolumeBlockData, and for initializing the moving bed module

Author: Chinedu Okoli and Robert Parker
"""
import os

from pyomo.environ import (ConcreteModel, SolverFactory,
                           TransformationFactory)

from idaes.core import FlowsheetBlock, EnergyBalanceType
from Models.moving_bed import MovingBed
from Models.gas_phase_thermo import Gas_Phase_Thermo_ParameterBlock
from Models.solid_phase_thermo import Solid_Phase_Thermo_ParameterBlock
from Models.hetero_reactions import HeteroReactionParameterBlock

from idaes.core.util.model_statistics import (
        degrees_of_freedom,
        activated_equalities_generator,
        unfixed_variables_in_activated_equalities_set)
from dyn_utils import (fix_initial_conditions)


# -----------------------------------------------------------------------------
# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None


# -----------------------------------------------------------------------------
def test_build(debug=False):
    # space grid
    nxfe = 10
    nxcp = 3
    xfe_set = [0]

    # time grid
    horizon = 10
    time_set = [0, horizon]
    ntfe = 2
    ntcp = 2

    # for each holdup(8), have a diff var for each point in space, except inlet
    n_diff_var = nxfe*nxcp*8

    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": True,
                                   "time_set": time_set})

    # Set up thermo props and reaction props
    m.fs.gas_properties = Gas_Phase_Thermo_ParameterBlock()
    m.fs.solid_properties = Solid_Phase_Thermo_ParameterBlock()
    m.fs.hetero_reactions = HeteroReactionParameterBlock(
            default={"property_package": m.fs.solid_properties})

    m.fs.MB = MovingBed(
            default={
                    "finite_elements": nxfe,
                    "collocation_points": nxcp,
                    "has_holdup": True,
                    "length_domain_set": xfe_set,
                    "transformation_method": "dae.collocation",
                    "transformation_scheme": "LAGRANGE-RADAU",
                    "energy_balance_type": EnergyBalanceType.enthalpyTotal,
                    "gas_phase_config":
                    {"property_package": m.fs.gas_properties,
                     "has_pressure_change": True,
                     "pressure_drop_type": "ergun_correlation"},
                    "solid_phase_config":
                    {"property_package": m.fs.solid_properties,
                     "reaction_package": m.fs.hetero_reactions
                     }})

    # Time discretization
    discretizer = TransformationFactory('dae.collocation')
    discretizer.apply_to(m, wrt=m.fs.time, nfe=ntfe, ncp=ntcp,
                         scheme='LAGRANGE-RADAU')

    nt = len(m.fs.time)
    nx = len(m.fs.MB.length_domain)

    assert m.fs.MB.config.flow_type == "counter_current"
    assert m.fs.MB.config.transformation_scheme == "LAGRANGE-RADAU"
    assert m.fs.MB.config.collocation_points == nxcp
    assert m.fs.MB.config.finite_elements == nxfe
    assert nx == nxfe*nxcp + 1
    assert nt == ntfe*ntcp + 1

    if debug:
        if not os.path.isdir('debug'):
            os.mkdir('debug')
        with open('debug/eqn_postdisc', 'w') as f:
            for con in activated_equalities_generator(m):
                f.write(con.name + '\n')
        with open('debug/var_postdisc', 'w') as f:
            for var in unfixed_variables_in_activated_equalities_set(m):
                f.write(var.name + '\n')

    # There should be (n_diff_var + nt*11 + 2) DOFs in this model:
    # Geometry - 2 (bed length and height)
    # Gas feed - 6*nt (inlet flow, temperature, pressure and mole fractions(3))
    # Solid feed - 5*nt (inlet flow, temperature, and mass fractions (3))
    # Initial conditions - n_diff_var
    assert degrees_of_freedom(m) == n_diff_var + nt*11 + 2

    fix_initial_conditions(m, m.fs.time)
    # Currently fixes initial conditions at inlets
    # => Initially, inlets for 8 'differential variables'
    #    (temp/comp) are no longer degrees of freedom.
    #    This may or may not be desired behavior.

    assert degrees_of_freedom(m) == nt*11 + 2 - 8


if __name__ == '__main__':
    test_build(debug=False)
