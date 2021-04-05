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
"""
Methods for eNRTL activity coefficient method.

Only applicable to liquid/electrolyte phases
"""
from pyomo.environ import Expression

from .eos_base import EoSBase
from .enrtl_submethods import ConstantAlpha, ConstantTau
import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


class ENRTL(EoSBase):
    # Add attribute indicating support for electrolyte systems
    electrolyte_support = True

    def build_parameters(b):

        # Check options for alpha rule
        if (b.config.equation_of_state_options is not None and
                "alpha_rule" in b.config.equation_of_state_options):
            b.config.equation_of_state_options[
                "alpha_rule"].build_parameters(b)
        else:
            ConstantAlpha.build_parameters(b)

        # Check options for tau rule
        if (b.config.equation_of_state_options is not None and
                "tau_rule" in b.config.equation_of_state_options):
            b.config.equation_of_state_options[
                "tau_rule"].build_parameters(b)
        else:
            ConstantTau.build_parameters(b)

    def common(b, pobj):
        def rule_X(b, j):
            if j in b.params.cation_set or j in b.params.anion_set:
                cobj = b.params.get_component(j)
                return (b.mole_frac_phase_comp_true[pobj.local_name, j] *
                        abs(cobj.config.charge))
            else:
                return b.mole_frac_phase_comp_true[pobj.local_name, j]
        b._X = Expression(b.params.true_species_set,
                          rule=rule_X,
                          doc="Charge x mole fraction term")

        def rule_Y(b, j):
            cobj = b.params.get_component(j)
            if cobj.config.charge < 0:
                # Anion
                dom = b.params.anion_set
            else:
                dom = b.params.cation_set

            return b._X[j]/sum(b._X[i] for i in dom)

        b._Y = Expression(b.params.ion_set,
                          rule=rule_Y,
                          doc="Charge composition")

    @staticmethod
    def calculate_scaling_factors(b, pobj):
        pass

    @staticmethod
    def dens_mol_phase(b, p):
        return 55e3

    @staticmethod
    def enth_mol_phase(b, p):
        return 1e2*b.temperature

    @staticmethod
    def enth_mol_phase_comp(b, p, j):
        return 1e2*b.temperature
