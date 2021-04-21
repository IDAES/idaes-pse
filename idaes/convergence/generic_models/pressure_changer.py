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
import idaes.core.util.convergence.convergence_base as cb
import pyomo.environ as pe

from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models.pressure_changer import (
    PressureChanger,
    ThermodynamicAssumption,
)

# Import property package for testing
from idaes.generic_models.properties import iapws95 as pp
'''
This module contains the code for convergence testing of the
PressureChanger model
'''

@cb.register_convergence_class("PressureChanger")
class PressureChangerConvergenceEvaluation(cb.ConvergenceEvaluation):
    def get_specification(self):
        """
        Returns the convergence evaluation specification for the
        PressureChanger unit model

        Returns
        -------
           ConvergenceEvaluationSpecification
        """
        s = cb.ConvergenceEvaluationSpecification()

        s.add_sampled_input(
                name='Inlet_Flowrate',
                pyomo_path='fs.pc.control_volume.properties_in[0].flow_mol',
                lower=1, upper=1e6, mean=5e3, std=5e3)

        s.add_sampled_input(
                name='Inlet_Pressure',
                pyomo_path='fs.pc.control_volume.properties_in[0].pressure',
                lower=1e5, upper=1e7, mean=1e6, std=3e6)
        return s

    def get_initialized_model(self):
        """
        Returns an initialized model for the PressureChanger unit model
        convergence evaluation

        Returns
        -------
           Pyomo model : returns a pyomo model of the PressureChanger unit
        """
        m = pe.ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.props = pp.Iapws95ParameterBlock()

        m.fs.pc = PressureChanger(default={
                "property_package": m.fs.props,
                "thermodynamic_assumption": ThermodynamicAssumption.isothermal})

        m.fs.pc.deltaP.fix(-1e3)
        m.fs.pc.inlet[:].flow_mol.fix(27.5e3)
        m.fs.pc.inlet[:].enth_mol.fix(4000)
        m.fs.pc.inlet[:].pressure.fix(2e6)

        init_state = {
                "flow_mol": 27.5e3,
                "pressure": 2e6,
                "enth_mol": 4000
                }

        m.fs.pc.initialize(state_args=init_state)

        # Create a solver for initialization
        opt = self.get_solver()
        opt.solve(m)

        # return the initialized model
        return m

    def get_solver(self):
        """
        Returns an ipopt solver object with the desired options for
        convergence evaluation (and initialization)
        Returns
        -------
           Pyomo solver
        """
        opt = pe.SolverFactory('ipopt')
        opt.options = {'tol': 1e-6,
                       'mu_init': 1e-8,
                       'bound_push': 1e-8}
        return opt
