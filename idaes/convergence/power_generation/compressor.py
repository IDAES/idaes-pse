
import pyomo.environ as pyo
import idaes.core as idaes_core
from idaes.power_generation.unit_models.helm import HelmIsentropicCompressor
import idaes.core.util.convergence.convergence_base as cb
from idaes.generic_models.properties import iapws95
from idaes.core.util import get_solver

from pyomo.environ import units as pyunits
import idaes

def create_isentropic_compressor(f=1000, T_in=500, p_in=1e6, ratioP=1.5):
    m = pyo.ConcreteModel()
    m.fs = idaes_core.FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = HelmIsentropicCompressor(default={"property_package": m.fs.properties})
    hin = iapws95.htpx(T_in*pyunits.K, p_in*pyunits.Pa) # J/mol
    m.fs.unit.inlet.flow_mol[0].fix(f)
    m.fs.unit.inlet.enth_mol[0].fix(hin)
    m.fs.unit.inlet.pressure[0].fix(p_in)
    m.fs.unit.ratioP[0].fix(ratioP)
    m.fs.unit.efficiency_isentropic.fix(0.9)
    m.fs.unit.initialize()
    return m

@cb.register_convergence_class("HelmIsentropicCompressor")
class HelmIsentropicCompressorConvergenceEvaluation(cb.ConvergenceEvaluation):
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
                pyomo_path='fs.unit.control_volume.properties_in[0].flow_mol',
                lower=1, upper=1e6, distribution="uniform")

        s.add_sampled_input(
                name='Inlet_Enthalpy',
                pyomo_path='fs.unit.control_volume.properties_in[0].enth_mol',
                lower=5000, upper=60000, distribution="uniform")

        s.add_sampled_input(
                name='Inlet_Pressure',
                pyomo_path='fs.unit.control_volume.properties_in[0].pressure',
                lower=1e5, upper=5e6, distribution="uniform")

        s.add_sampled_input(
                name='Pressure_Ratio',
                pyomo_path='fs.unit.ratioP[0]',
                lower=1, upper=4, distribution="uniform")

        s.add_sampled_input(
                name='Efficiency',
                pyomo_path='fs.unit.efficiency_isentropic[0]',
                lower=0.2, upper=1, distribution="uniform")

        return s

    def get_initialized_model(self):
        """
        Returns an initialized model for the PressureChanger unit model
        convergence evaluation

        Returns
        -------
           Pyomo model : returns a pyomo model of the PressureChanger unit
        """
        return create_isentropic_compressor()

    def get_solver(self):
        """
        Returns a solver object with the desired options for
        convergence evaluation (and initialization)

        Returns
        -------
           Pyomo solver
        """
        opt = get_solver(options={"max_iter": 25})
        return opt
