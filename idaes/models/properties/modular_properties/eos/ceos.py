#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Methods for cubic equations of state.

Currently only supports liquid and vapor phases
"""
from enum import Enum
from copy import deepcopy

from pyomo.environ import (
    exp,
    Expression,
    log,
    Reals,
    sqrt,
    Var,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In
from idaes.models.properties.modular_properties.base.utility import (
    get_method,
    get_component_object as cobj,
)
from idaes.core.util.math import safe_log
from .eos_base import EoSBase
import idaes.logger as idaeslog
from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    PropertyNotSupportedError,
)
from idaes.models.properties.modular_properties.eos.ceos_common import (
    EoS_param,
    cubic_roots_available,
    CubicThermoExpressions,
    CubicType,
)

# pylint: disable=invalid-name

# Set up logger
_log = idaeslog.getLogger(__name__)


"""
References:
[1]. Poling, B.E., Prausnitz, J.M. and O’connell, J.P., 2001.
     Properties of gases and liquids. McGraw-Hill Education.

[2]. Trujillo, M.F., O'Rourke, P. and Torres, D., 2002.
     Generalizing the Thermodynamics State Relationships in KIVA-3V
     https://www.osti.gov/servlets/purl/809947 (Last accessed: 08/13/2021)
"""


class MixingRuleA(Enum):
    # Rule to calculate am for cubic equations of state
    default = 0


class MixingRuleB(Enum):
    # Rule to calculate bm for cubic equations of state
    default = 0


# Value for smoothing epsilon for SafeLog when used
eps_SL = 1e-8

CubicConfig = ConfigBlock()
CubicConfig.declare(
    "type",
    ConfigValue(
        domain=In(CubicType),
        description="Equation of state to use",
        doc="Enum indicating type of cubic equation of state to use.",
    ),
)
CubicConfig.declare(
    "mixing_rule_a",
    ConfigValue(
        default=MixingRuleA.default,
        domain=In(MixingRuleA),
        description="Mixing rule for parameter am for cubic EoS",
        doc="Enum indicating type of mixing rule for parameter am to use.",
    ),
)
CubicConfig.declare(
    "mixing_rule_b",
    ConfigValue(
        default=MixingRuleB.default,
        domain=In(MixingRuleB),
        description="Mixing rule for parameter bm for cubic EoS",
        doc="Enum indicating type of mixing rule for parameter bm to use.",
    ),
)


class Cubic(EoSBase):
    @staticmethod
    def common(b, pobj):
        # TODO: determine if Henry's Law applies to Cubic EoS systems
        # For now, raise an exception if found
        # Follow on questions:
        # If Henry's law is used for a component, how does that effect
        # calculating A, B and phi?
        for j in b.component_list:
            cobj = b.params.get_component(j)
            if (
                cobj.config.henry_component is not None
                and pobj.local_name in cobj.config.henry_component
            ):
                raise PropertyNotSupportedError(
                    "{} Cubic equations of state do not support Henry's "
                    "components [{}, {}].".format(b.name, pobj.local_name, j)
                )

        ctype = pobj._cubic_type
        cname = pobj.config.equation_of_state_options["type"].name
        mixing_rule_a = pobj._mixing_rule_a
        mixing_rule_b = pobj._mixing_rule_b

        if hasattr(b, cname + "_fw"):
            # Common components already constructed by previous phase
            return

        # Create expressions for coefficients
        def rule_fw(m, j):
            func_fw = getattr(m.params, cname + "_func_fw")
            cobj = m.params.get_component(j)
            return func_fw(cobj)

        b.add_component(
            cname + "_fw",
            Expression(b.component_list, rule=rule_fw, doc="EoS S factor"),
        )

        def rule_a_crit(m, j):
            cobj = m.params.get_component(j)
            return EoS_param[ctype]["omegaA"] * (
                (Cubic.gas_constant(b) * cobj.temperature_crit) ** 2
                / cobj.pressure_crit
            )

        b.add_component(
            cname + "_a_crit",
            Expression(
                b.component_list,
                rule=rule_a_crit,
                doc="Component a coefficient at T_crit",
            ),
        )

        def rule_a(m, j):
            cobj = m.params.get_component(j)
            fw = getattr(m, cname + "_fw")[j]
            ac = getattr(m, cname + "_a_crit")[j]
            func_alpha = getattr(m.params, cname + "_func_alpha")

            return ac * func_alpha(m.temperature, fw, cobj)

        b.add_component(
            cname + "_a",
            Expression(b.component_list, rule=rule_a, doc="Component a coefficient"),
        )

        def rule_da_dT(m, j):
            cobj = m.params.get_component(j)
            fw = getattr(m, cname + "_fw")[j]
            ac = getattr(m, cname + "_a_crit")[j]
            func_dalpha_dT = getattr(m.params, cname + "_func_dalpha_dT")

            return ac * func_dalpha_dT(m.temperature, fw, cobj)

        b.add_component(
            cname + "_da_dT",
            Expression(
                b.component_list,
                rule=rule_da_dT,
                doc="Temperature derivative of component a",
            ),
        )

        def rule_d2a_dT2(m, j):
            cobj = m.params.get_component(j)
            fw = getattr(m, cname + "_fw")[j]
            ac = getattr(m, cname + "_a_crit")[j]
            func_d2alpha_dT2 = getattr(m.params, cname + "_func_d2alpha_dT2")

            return ac * func_d2alpha_dT2(m.temperature, fw, cobj)

        b.add_component(
            cname + "_d2a_dT2",
            Expression(
                b.component_list,
                rule=rule_d2a_dT2,
                doc="Second temperature derivative" "of component a",
            ),
        )

        def func_b(m, j):
            cobj = m.params.get_component(j)
            return (
                EoS_param[ctype]["coeff_b"]
                * Cubic.gas_constant(b)
                * cobj.temperature_crit
                / cobj.pressure_crit
            )

        b.add_component(
            cname + "_b",
            Expression(b.component_list, rule=func_b, doc="Component b coefficient"),
        )

        if mixing_rule_a == MixingRuleA.default:

            def rule_am(m, p):
                a = getattr(m, cname + "_a")
                return rule_am_default(m, cname, a, p)

            b.add_component(cname + "_am", Expression(b.phase_list, rule=rule_am))

            def rule_daij_dT(m, i, j):
                a = getattr(m, cname + "_a")
                da_dT = getattr(m, cname + "_da_dT")
                k = getattr(m.params, cname + "_kappa")

                # Include temperature derivative of k for future extension
                dk_ij_dT = 0

                return sqrt(a[i] * a[j]) * (
                    -dk_ij_dT + (1 - k[i, j]) / 2 * (da_dT[i] / a[i] + da_dT[j] / a[j])
                )

            b.add_component(
                cname + "_daij_dT",
                Expression(b.component_list, b.component_list, rule=rule_daij_dT),
            )

            def rule_dam_dT(m, p):
                daij_dT = getattr(m, cname + "_daij_dT")
                return sum(
                    sum(
                        m.mole_frac_phase_comp[p, i]
                        * m.mole_frac_phase_comp[p, j]
                        * daij_dT[i, j]
                        for j in m.components_in_phase(p)
                    )
                    for i in m.components_in_phase(p)
                )

            b.add_component(
                cname + "_dam_dT", Expression(b.phase_list, rule=rule_dam_dT)
            )

            def rule_d2am_dT2(m, p):
                k = getattr(m.params, cname + "_kappa")
                a = getattr(m, cname + "_a")
                da_dT = getattr(m, cname + "_da_dT")
                d2a_dT2 = getattr(m, cname + "_d2a_dT2")
                # Placeholders for if temperature dependent k is needed
                dk_dT = 0
                d2k_dT2 = 0

                # Initialize loop variable
                d2am_dT2 = 0

                for i in m.components_in_phase(p):
                    for j in m.components_in_phase(p):
                        d2aij_dT2 = sqrt(a[i] * a[j]) * (
                            -d2k_dT2
                            - dk_dT * (da_dT[i] / a[i] + da_dT[j] / a[j])
                            + (1 - k[i, j])
                            / 2
                            * (
                                d2a_dT2[i] / a[i]
                                + d2a_dT2[j] / a[j]
                                - 1 / 2 * (da_dT[i] / a[i] - da_dT[j] / a[j]) ** 2
                            )
                        )
                        d2am_dT2 += (
                            m.mole_frac_phase_comp[p, i]
                            * m.mole_frac_phase_comp[p, j]
                            * d2aij_dT2
                        )
                return d2am_dT2

            b.add_component(
                cname + "_d2am_dT2", Expression(b.phase_list, rule=rule_d2am_dT2)
            )

            def rule_delta(m, p, i):
                # See pg. 145 in Properties of Gases and Liquids
                a = getattr(m, cname + "_a")
                am = getattr(m, cname + "_am")
                kappa = getattr(m.params, cname + "_kappa")
                return (
                    2
                    * sqrt(a[i])
                    / am[p]
                    * sum(
                        m.mole_frac_phase_comp[p, j] * sqrt(a[j]) * (1 - kappa[i, j])
                        for j in b.components_in_phase(p)
                    )
                )

            b.add_component(
                cname + "_delta", Expression(b.phase_component_set, rule=rule_delta)
            )

        else:
            raise ConfigurationError(
                "{} Unrecognized option for Equation of State "
                "mixing_rule_a: {}. Must be an instance of MixingRuleA "
                "Enum.".format(b.name, mixing_rule_a)
            )

        if mixing_rule_b == MixingRuleB.default:

            def rule_bm(m, p):
                b = getattr(m, cname + "_b")
                return rule_bm_default(m, b, p)

            b.add_component(cname + "_bm", Expression(b.phase_list, rule=rule_bm))
        else:
            raise ConfigurationError(
                "{} Unrecognized option for Equation of State "
                "mixing_rule_a: {}. Must be an instance of MixingRuleB "
                "Enum.".format(b.name, mixing_rule_b)
            )

        def rule_A(m, p):
            am = getattr(m, cname + "_am")
            return am[p] * m.pressure / (Cubic.gas_constant(b) * m.temperature) ** 2

        b.add_component(cname + "_A", Expression(b.phase_list, rule=rule_A))

        def rule_B(m, p):
            bm = getattr(m, cname + "_bm")
            return bm[p] * m.pressure / (Cubic.gas_constant(b) * m.temperature)

        b.add_component(cname + "_B", Expression(b.phase_list, rule=rule_B))

        # Add components at equilibrium state if required
        if b.params.config.phases_in_equilibrium is not None and (
            not b.config.defined_state or b.always_flash
        ):

            def func_a_eq(m, p1, p2, j):
                cobj = m.params.get_component(j)
                fw = getattr(m, cname + "_fw")[j]
                ac = getattr(m, cname + "_a_crit")[j]
                func_alpha = getattr(m.params, cname + "_func_alpha")

                return ac * func_alpha(m._teq[p1, p2], fw, cobj)

            b.add_component(
                "_" + cname + "_a_eq",
                Expression(
                    b.params._pe_pairs,
                    b.component_list,
                    rule=func_a_eq,
                    doc="Component a coefficient at Teq",
                ),
            )

            def rule_am_eq(m, p1, p2, p3):
                try:
                    rule = m.params.get_phase(p3).config.equation_of_state_options[
                        "mixing_rule_a"
                    ]
                except (KeyError, TypeError):
                    rule = MixingRuleA.default

                a = getattr(m, "_" + cname + "_a_eq")
                if rule == MixingRuleA.default:
                    return rule_am_default(m, cname, a, p3, (p1, p2))
                else:
                    raise ConfigurationError(
                        "{} Unrecognized option for Equation of State "
                        "mixing_rule_a: {}. Must be an instance of MixingRuleA "
                        "Enum.".format(m.name, rule)
                    )

            b.add_component(
                "_" + cname + "_am_eq",
                Expression(b.params._pe_pairs, b.phase_list, rule=rule_am_eq),
            )

            def rule_A_eq(m, p1, p2, p3):
                am_eq = getattr(m, "_" + cname + "_am_eq")
                return (
                    am_eq[p1, p2, p3]
                    * m.pressure
                    / (Cubic.gas_constant(b) * m._teq[p1, p2]) ** 2
                )

            b.add_component(
                "_" + cname + "_A_eq",
                Expression(b.params._pe_pairs, b.phase_list, rule=rule_A_eq),
            )

            def rule_B_eq(m, p1, p2, p3):
                bm = getattr(m, cname + "_bm")
                return bm[p3] * m.pressure / (Cubic.gas_constant(b) * m._teq[p1, p2])

            b.add_component(
                "_" + cname + "_B_eq",
                Expression(b.params._pe_pairs, b.phase_list, rule=rule_B_eq),
            )

            def rule_delta_eq(m, p1, p2, p3, i):
                # See pg. 145 in Properties of Gases and Liquids
                a = getattr(m, "_" + cname + "_a_eq")
                am = getattr(m, "_" + cname + "_am_eq")
                kappa = getattr(m.params, cname + "_kappa")
                return (
                    2
                    * sqrt(a[p1, p2, i])
                    / am[p1, p2, p3]
                    * sum(
                        m.mole_frac_phase_comp[p3, j]
                        * sqrt(a[p1, p2, j])
                        * (1 - kappa[i, j])
                        for j in m.components_in_phase(p3)
                    )
                )

            b.add_component(
                "_" + cname + "_delta_eq",
                Expression(
                    b.params._pe_pairs, b.phase_component_set, rule=rule_delta_eq
                ),
            )

    @staticmethod
    def calculate_scaling_factors(b, pobj):
        pass

    @staticmethod
    def build_parameters(b):
        param_block = b.parent_block()
        if not (b.is_vapor_phase() or b.is_liquid_phase()):
            raise PropertyNotSupportedError(
                "{} received unrecognized phase "
                "name {}. Cubic equation of state supports only Vap and Liq "
                "phases.".format(param_block.name, b)
            )

        if b.config.equation_of_state_options["type"] not in set(
            item for item in CubicType
        ):
            raise ConfigurationError(
                "{} Unrecognized option for equation of "
                "state type: {}. Must be an instance of CubicType "
                "Enum.".format(b.name, b.config.equation_of_state_options["type"])
            )

        ctype = b.config.equation_of_state_options["type"]
        b._cubic_type = ctype
        cname = ctype.name

        # Check to see if ConfigBlock was created by previous phase
        if hasattr(param_block, cname + "_eos_options"):
            ConfigBlock = getattr(param_block, cname + "_eos_options")
            for key, value in b.config.equation_of_state_options.items():
                if ConfigBlock[key] != value:
                    raise ConfigurationError(
                        "In {}, different {} equation of "
                        "state options for {} are set in different phases, which is "
                        "not supported.".format(b.name, cname, key)
                    )
            # Once the options have been validated, we don't have anything
            # left to do
            mixing_rule_a = ConfigBlock["mixing_rule_a"]
            mixing_rule_b = ConfigBlock["mixing_rule_b"]
            b._mixing_rule_a = mixing_rule_a
            b._mixing_rule_b = mixing_rule_b
            return

        setattr(param_block, cname + "_eos_options", deepcopy(CubicConfig))
        ConfigBlock = getattr(param_block, cname + "_eos_options")
        ConfigBlock.set_value(b.config.equation_of_state_options)

        mixing_rule_a = ConfigBlock["mixing_rule_a"]
        mixing_rule_b = ConfigBlock["mixing_rule_b"]
        b._mixing_rule_a = mixing_rule_a
        b._mixing_rule_b = mixing_rule_b

        kappa_data = param_block.config.parameter_data[cname + "_kappa"]
        param_block.add_component(
            cname + "_kappa",
            Var(
                param_block.component_list,
                param_block.component_list,
                within=Reals,
                initialize=kappa_data,
                doc=cname + " binary interaction parameters",
                units=None,
            ),
        )

        if b._cubic_type == CubicType.PR:
            func_fw = func_fw_PR
        elif b._cubic_type == CubicType.SRK:
            func_fw = func_fw_SRK
        else:
            raise BurntToast(
                "{} received unrecognized cubic type. This should "
                "never happen, so please contact the IDAES developers "
                "with this bug.".format(b.name)
            )
        setattr(param_block, cname + "_func_fw", func_fw)
        setattr(param_block, cname + "_func_alpha", func_alpha_soave)
        setattr(param_block, cname + "_func_dalpha_dT", func_dalpha_dT_soave)
        setattr(param_block, cname + "_func_d2alpha_dT2", func_d2alpha_dT2_soave)

    @staticmethod
    def compress_fact_phase(b, p):
        pobj = b.params.get_phase(p)
        cname = pobj._cubic_type.name
        A = getattr(b, cname + "_A")
        B = getattr(b, cname + "_B")

        expr_write = CubicThermoExpressions(b)
        if pobj.is_vapor_phase():
            return expr_write.z_vap(eos=pobj._cubic_type, A=A[p], B=B[p])
        elif pobj.is_liquid_phase():
            return expr_write.z_liq(eos=pobj._cubic_type, A=A[p], B=B[p])
        raise BurntToast(
            "{} non-vapor or liquid phase called for cubic "
            "EoS compressability factor. This should never "
            "happen, so please contact the IDAES developers "
            "with this bug.".format(b.name)
        )

    @staticmethod
    def cp_mol_phase(blk, p):
        pobj = blk.params.get_phase(p)
        cname = pobj._cubic_type.name

        am = getattr(blk, cname + "_am")[p]
        bm = getattr(blk, cname + "_bm")[p]
        B = getattr(blk, cname + "_B")[p]

        dam_dT = getattr(blk, cname + "_dam_dT")[p]
        d2am_dT2 = getattr(blk, cname + "_d2am_dT2")[p]

        T = blk.temperature
        R = Cubic.gas_constant(blk)
        Z = blk.compress_fact_phase[p]
        dZdT = _dZ_dT(blk, p)

        EoS_u = EoS_param[pobj._cubic_type]["u"]
        EoS_w = EoS_param[pobj._cubic_type]["w"]
        EoS_p = sqrt(EoS_u**2 - 4 * EoS_w)

        expression1 = 2 * Z + (EoS_u + EoS_p) * B
        expression2 = 2 * Z + (EoS_u - EoS_p) * B
        expression3 = B * (dZdT + Z / T) / (Z**2 + Z * EoS_u * B + EoS_w * B**2)

        cp_ideal_gas = sum(
            blk.mole_frac_phase_comp[p, j]
            * get_method(blk, "cp_mol_ig_comp", j)(blk, cobj(blk, j), T)
            for j in blk.components_in_phase(p)
        )

        # Derived from the relations in Chapter 6 of [1]
        cp_departure = (
            R * (T * dZdT + Z - 1)
            + (T * d2am_dT2 / (EoS_p * bm))
            * safe_log(expression1 / expression2, eps=eps_SL)
            + ((am - T * dam_dT) * expression3 / bm)
        )

        return cp_ideal_gas + cp_departure

    @staticmethod
    def cv_mol_phase(blk, p):
        pobj = blk.params.get_phase(p)
        cname = pobj._cubic_type.name
        am = getattr(blk, cname + "_am")[p]
        bm = getattr(blk, cname + "_bm")[p]
        cp = blk.cp_mol_phase[p]
        V = 1 / blk.dens_mol_phase[p]
        dam_dT = getattr(blk, cname + "_dam_dT")[p]

        EoS_u = EoS_param[pobj._cubic_type]["u"]
        EoS_w = EoS_param[pobj._cubic_type]["w"]

        dPdV = -((Cubic.gas_constant(blk) * blk.temperature) / (V - bm) ** 2) + (
            am * (2 * V + EoS_u * bm) / (V**2 + EoS_u * bm * V + EoS_w * bm**2) ** 2
        )

        dPdT = (Cubic.gas_constant(blk) / (V - bm)) - (
            1 / (V**2 + EoS_u * bm * V + EoS_w * bm**2)
        ) * dam_dT

        # See Chapter 6 in [1]
        return cp + blk.temperature * dPdT**2 / dPdV

    @staticmethod
    def dens_mass_phase(b, p):
        return b.dens_mol_phase[p] * b.mw_phase[p]

    @staticmethod
    def dens_mol_phase(b, p):
        pobj = b.params.get_phase(p)
        return b.pressure / (
            Cubic.gas_constant(b) * b.temperature * b.compress_fact_phase[p]
        )

    @staticmethod
    def energy_internal_mol_phase(blk, p):
        pobj = blk.params.get_phase(p)

        cname = pobj._cubic_type.name
        am = getattr(blk, cname + "_am")[p]
        bm = getattr(blk, cname + "_bm")[p]
        B = getattr(blk, cname + "_B")[p]
        dam_dT = getattr(blk, cname + "_dam_dT")[p]
        Z = blk.compress_fact_phase[p]

        EoS_u = EoS_param[pobj._cubic_type]["u"]
        EoS_w = EoS_param[pobj._cubic_type]["w"]
        EoS_p = sqrt(EoS_u**2 - 4 * EoS_w)

        # Derived from equation on pg. 120 in Properties of Gases and Liquids
        # Departure function for U is similar to H minus the RT(Z-1) term
        return (
            (blk.temperature * dam_dT - am)
            * safe_log(
                (2 * Z + B * (EoS_u + EoS_p)) / (2 * Z + B * (EoS_u - EoS_p)),
                eps=eps_SL,
            )
        ) / (bm * EoS_p) + sum(
            blk.mole_frac_phase_comp[p, j]
            * EoSBase.energy_internal_mol_ig_comp_pure(blk, j)
            for j in blk.components_in_phase(p)
        )

    @staticmethod
    def energy_internal_mol_phase_comp(blk, p, j):
        pobj = blk.params.get_phase(p)

        return (
            blk.enth_mol_phase_comp[p, j] - blk.pressure * blk.vol_mol_phase_comp[p, j]
        )

    @staticmethod
    def enth_mol_phase(blk, p):
        pobj = blk.params.get_phase(p)

        cname = pobj._cubic_type.name
        am = getattr(blk, cname + "_am")[p]
        bm = getattr(blk, cname + "_bm")[p]
        B = getattr(blk, cname + "_B")[p]
        dam_dT = getattr(blk, cname + "_dam_dT")[p]
        Z = blk.compress_fact_phase[p]
        R = Cubic.gas_constant(blk)
        T = blk.temperature

        EoS_u = EoS_param[pobj._cubic_type]["u"]
        EoS_w = EoS_param[pobj._cubic_type]["w"]
        EoS_p = sqrt(EoS_u**2 - 4 * EoS_w)

        enth_ideal = sum(
            blk.mole_frac_phase_comp[p, j]
            * get_method(blk, "enth_mol_ig_comp", j)(blk, cobj(blk, j), blk.temperature)
            for j in blk.components_in_phase(p)
        )

        # Derived from equation on pg. 120 in Properties of Gases and Liquids
        enth_departure = R * T * (Z - 1) + (T * dam_dT - am) / (bm * EoS_p) * safe_log(
            (2 * Z + B * (EoS_u + EoS_p)) / (2 * Z + B * (EoS_u - EoS_p)), eps=eps_SL
        )
        return enth_ideal + enth_departure

    @staticmethod
    def enth_mol_phase_comp(blk, p, j):
        pobj = blk.params.get_phase(p)

        dlogphi_j_dT = _d_log_fug_coeff_dT_phase_comp(blk, p, j)

        enth_ideal_gas = get_method(blk, "enth_mol_ig_comp", j)(
            blk, cobj(blk, j), blk.temperature
        )

        enth_departure = -Cubic.gas_constant(blk) * blk.temperature**2 * dlogphi_j_dT

        return enth_ideal_gas + enth_departure

    @staticmethod
    def entr_mol_phase(blk, p):
        pobj = blk.params.get_phase(p)

        cname = pobj._cubic_type.name
        bm = getattr(blk, cname + "_bm")[p]
        B = getattr(blk, cname + "_B")[p]
        dam_dT = getattr(blk, cname + "_dam_dT")[p]
        Z = blk.compress_fact_phase[p]

        EoS_u = EoS_param[pobj._cubic_type]["u"]
        EoS_w = EoS_param[pobj._cubic_type]["w"]
        EoS_p = sqrt(EoS_u**2 - 4 * EoS_w)

        R = Cubic.gas_constant(blk)

        entr_ideal_gas = -R * safe_log(
            blk.pressure / blk.params.pressure_ref, eps=eps_SL
        )
        for j in blk.components_in_phase(p):
            entr_j = get_method(blk, "entr_mol_ig_comp", j)(
                blk, cobj(blk, j), blk.temperature
            )
            xj = blk.mole_frac_phase_comp[p, j]
            log_xj = blk.log_mole_frac_phase_comp[p, j]

            entr_ideal_gas += xj * (entr_j - R * log_xj)

        # See pg. 102 in Properties of Gases and Liquids
        # or pg. 208 of Sandler, 4th Ed.
        entr_departure = R * safe_log((Z - B), eps=eps_SL) + dam_dT / (
            bm * EoS_p
        ) * safe_log(
            (2 * Z + B * (EoS_u + EoS_p)) / (2 * Z + B * (EoS_u - EoS_p)), eps=eps_SL
        )

        return entr_ideal_gas + entr_departure

    @staticmethod
    def entr_mol_phase_comp(blk, p, j):
        pobj = blk.params.get_phase(p)

        logphi_j = _log_fug_coeff_phase_comp(blk, p, j)
        dlogphi_j_dT = _d_log_fug_coeff_dT_phase_comp(blk, p, j)

        R = Cubic.gas_constant(blk)

        entr_ideal_gas = get_method(blk, "entr_mol_ig_comp", j)(
            blk, cobj(blk, j), blk.temperature
        ) - R * (
            safe_log(blk.pressure / blk.params.pressure_ref, eps=eps_SL)
            + blk.log_mole_frac_phase_comp[p, j]
        )

        entr_departure = -R * logphi_j - R * blk.temperature * dlogphi_j_dT

        return entr_ideal_gas + entr_departure

    @staticmethod
    def fug_phase_comp(b, p, j):
        return b.mole_frac_phase_comp[p, j] * b.pressure * b.fug_coeff_phase_comp[p, j]

    @staticmethod
    def fug_phase_comp_eq(b, p, j, pp):
        return (
            b.mole_frac_phase_comp[p, j]
            * b.pressure
            * exp(_log_fug_coeff_phase_comp_eq(b, p, j, pp))
        )

    @staticmethod
    def log_fug_phase_comp_eq(b, p, j, pp):
        return (
            b.log_mole_frac_phase_comp[p, j]
            + log(b.pressure / b.params.pressure_ref)
            + _log_fug_coeff_phase_comp_eq(b, p, j, pp)
        )

    @staticmethod
    def fug_coeff_phase_comp(blk, p, j):
        pobj = blk.params.get_phase(p)
        ctype = pobj._cubic_type

        cname = pobj._cubic_type.name
        b = getattr(blk, cname + "_b")[j]
        bm = getattr(blk, cname + "_bm")[p]
        A = getattr(blk, cname + "_A")[p]
        B = getattr(blk, cname + "_B")[p]
        delta = getattr(blk, cname + "_delta")[p, j]
        Z = blk.compress_fact_phase[p]

        return exp(_log_fug_coeff_method(A, b, bm, B, delta, Z, ctype))

    @staticmethod
    def fug_coeff_phase_comp_eq(blk, p, j, pp):
        return exp(_log_fug_coeff_phase_comp_eq(blk, p, j, pp))

    @staticmethod
    def log_fug_phase_comp_Tbub(blk, p, j, pp):
        return _bubble_dew_log_fug_coeff_method(blk, p, j, pp, blk.temperature_bubble)

    @staticmethod
    def log_fug_phase_comp_Tdew(blk, p, j, pp):
        return _bubble_dew_log_fug_coeff_method(blk, p, j, pp, blk.temperature_dew)

    @staticmethod
    def log_fug_phase_comp_Pbub(blk, p, j, pp):
        return _bubble_dew_log_fug_coeff_method(blk, p, j, pp, blk.pressure_bubble)

    @staticmethod
    def log_fug_phase_comp_Pdew(blk, p, j, pp):
        return _bubble_dew_log_fug_coeff_method(blk, p, j, pp, blk.pressure_dew)

    @staticmethod
    def gibbs_mol_phase(b, p):
        return b.enth_mol_phase[p] - b.entr_mol_phase[p] * b.temperature

    @staticmethod
    def gibbs_mol_phase_comp(blk, p, j):
        # Calling the enthalpy and entropy directly adds a lot of overhead
        # because expressions involving the derivative of the fugacity coefficient
        # are generated. Those terms cancel mathematically, but I suspect Pyomo
        # leaves them in, causing trouble.
        R = Cubic.gas_constant(blk)
        T = blk.temperature
        logphi_j = _log_fug_coeff_phase_comp(blk, p, j)

        enth_ideal_gas = get_method(blk, "enth_mol_ig_comp", j)(blk, cobj(blk, j), T)

        entr_ideal_gas = get_method(blk, "entr_mol_ig_comp", j)(
            blk, cobj(blk, j), T
        ) - R * (
            safe_log(blk.pressure / blk.params.pressure_ref, eps=eps_SL)
            + blk.log_mole_frac_phase_comp[p, j]
        )
        gibbs_ideal_gas = enth_ideal_gas - T * entr_ideal_gas
        gibbs_departure = R * T * logphi_j

        return gibbs_ideal_gas + gibbs_departure

    @staticmethod
    def isentropic_speed_sound_phase(blk, p):
        # See Reference [2]
        return (
            sqrt(blk.heat_capacity_ratio_phase[p]) * blk.isothermal_speed_sound_phase[p]
        )

    @staticmethod
    def isothermal_speed_sound_phase(blk, p):
        pobj = blk.params.get_phase(p)
        cname = pobj._cubic_type.name
        am = getattr(blk, cname + "_am")[p]
        bm = getattr(blk, cname + "_bm")[p]
        V = 1 / blk.dens_mol_phase[p]
        mw = blk.mw
        rho = blk.dens_mass_phase[p]

        EoS_u = EoS_param[pobj._cubic_type]["u"]
        EoS_w = EoS_param[pobj._cubic_type]["w"]

        dPdV = -((Cubic.gas_constant(blk) * blk.temperature) / (V - bm) ** 2) + (
            am * (2 * V + EoS_u * bm) / (V**2 + EoS_u * bm * V + EoS_w * bm**2) ** 2
        )

        # see reference [2]
        return sqrt(-dPdV * mw / rho**2)

    @staticmethod
    def vol_mol_phase(b, p):
        return (
            Cubic.gas_constant(b)
            * b.temperature
            * b.compress_fact_phase[p]
            / b.pressure
        )

    @staticmethod
    def vol_mol_phase_comp(b, p, j):
        return (
            Cubic.gas_constant(b)
            * b.temperature
            / b.pressure
            * (b.compress_fact_phase[p] + _N_dZ_dNj(b, p, j))
        )


def _dZ_dT(blk, p):
    pobj = blk.params.get_phase(p)
    cname = pobj._cubic_type.name
    am = getattr(blk, cname + "_am")[p]
    A = getattr(blk, cname + "_A")[p]
    B = getattr(blk, cname + "_B")[p]
    dam_dT = getattr(blk, cname + "_dam_dT")[p]
    Z = blk.compress_fact_phase[p]
    T = blk.temperature

    EoS_u = EoS_param[pobj._cubic_type]["u"]
    EoS_w = EoS_param[pobj._cubic_type]["w"]

    dBdT = -B / T
    dAdT = (A / am) * dam_dT - (2 * A / T)

    K2 = (EoS_u - 1) * B - 1
    K3 = A - EoS_u * B - (EoS_u - EoS_w) * B**2
    K4 = -(A * B + EoS_w * B**2 + EoS_w * B**3)

    dK2dT = (EoS_u - 1) * dBdT
    dK3dT = dAdT - EoS_u * dBdT - 2 * (EoS_u - EoS_w) * B * dBdT
    dK4dT = -(A * dBdT + B * dAdT + 2 * EoS_w * B * dBdT + 3 * EoS_w * B**2 * dBdT)

    return -(Z**2 * dK2dT + Z * dK3dT + dK4dT) / (3 * Z**2 + 2 * K2 * Z + K3)


def _N_dZ_dNj(blk, p, j):
    pobj = blk.params.get_phase(p)
    cname = pobj._cubic_type.name

    if not (
        pobj._mixing_rule_a == MixingRuleA.default
        and pobj._mixing_rule_b == MixingRuleB.default
    ):
        # Any user adding more mixing rules will need to explicitly add
        # support for these functions
        raise NotImplementedError(
            "Block {} called for a property "
            "that is not supported by this choice "
            "of mixing rules.".format(blk.name)
        )

    a = getattr(blk, cname + "_a")
    b = getattr(blk, cname + "_b")
    k = getattr(blk.params, cname + "_kappa")
    am = getattr(blk, cname + "_am")[p]
    bm = getattr(blk, cname + "_bm")[p]
    A = getattr(blk, cname + "_A")[p]
    B = getattr(blk, cname + "_B")[p]

    if pobj._mixing_rule_a == MixingRuleA.default:
        N_dam_dNj = 2 * (
            -am
            + sum(
                blk.mole_frac_phase_comp[p, i] * (1 - k[i, j]) * sqrt(a[i] * a[j])
                for i in blk.components_in_phase(p)
            )
        )

    if pobj._mixing_rule_b == MixingRuleB.default:
        N_dbm_dNj = b[j] - bm

    Z = blk.compress_fact_phase[p]
    R = Cubic.gas_constant(blk)
    T = blk.temperature
    P = blk.pressure

    EoS_u = EoS_param[pobj._cubic_type]["u"]
    EoS_w = EoS_param[pobj._cubic_type]["w"]

    N_dA_dNj = P / (R * T) ** 2 * N_dam_dNj
    N_dB_dNj = P / (R * T) * N_dbm_dNj

    K2 = (EoS_u - 1) * B - 1
    K3 = A - EoS_u * B - (EoS_u - EoS_w) * B**2
    K4 = -(A * B + EoS_w * B**2 + EoS_w * B**3)

    N_dK2_dNj = (EoS_u - 1) * N_dB_dNj
    N_dK3_dNj = N_dA_dNj - EoS_u * N_dB_dNj - 2 * (EoS_u - EoS_w) * B * N_dB_dNj
    N_dK4_dNj = -(
        A * N_dB_dNj
        + B * N_dA_dNj
        + 2 * EoS_w * B * N_dB_dNj
        + 3 * EoS_w * B**2 * N_dB_dNj
    )

    return -(Z**2 * N_dK2_dNj + Z * N_dK3_dNj + N_dK4_dNj) / (
        3 * Z**2 + 2 * K2 * Z + K3
    )


def _log_fug_coeff_phase_comp_eq(blk, p, j, pp):
    pobj = blk.params.get_phase(p)

    cname = pobj._cubic_type.name
    b = getattr(blk, cname + "_b")
    bm = getattr(blk, cname + "_bm")
    Aeq = getattr(blk, "_" + cname + "_A_eq")
    Beq = getattr(blk, "_" + cname + "_B_eq")
    delta_eq = getattr(blk, "_" + cname + "_delta_eq")

    expr_write = CubicThermoExpressions(blk)
    if pobj.is_vapor_phase():

        def Zeq(p):
            return expr_write.z_vap(eos=pobj._cubic_type, A=Aeq[pp, p], B=Beq[pp, p])

    elif pobj.is_liquid_phase():

        def Zeq(p):
            return expr_write.z_liq(eos=pobj._cubic_type, A=Aeq[pp, p], B=Beq[pp, p])

    return _log_fug_coeff_method(
        Aeq[pp, p],
        b[j],
        bm[p],
        Beq[pp, p],
        delta_eq[pp, p, j],
        Zeq(p),
        pobj._cubic_type,
    )


def _log_fug_coeff_phase_comp(blk, p, j):
    pobj = blk.params.get_phase(p)

    if not (
        pobj._mixing_rule_a == MixingRuleA.default
        and pobj._mixing_rule_b == MixingRuleB.default
    ):
        # Any
        raise NotImplementedError(
            "Block {} called for a property "
            "that is not supported by this choice "
            "of mixing rules.".format(blk.name)
        )

    cname = pobj._cubic_type.name
    b = getattr(blk, cname + "_b")
    bm = getattr(blk, cname + "_bm")
    A = getattr(blk, cname + "_A")
    B = getattr(blk, cname + "_B")
    delta = getattr(blk, cname + "_delta")

    expr_write = CubicThermoExpressions(blk)
    if pobj.is_vapor_phase():

        def Z(p):
            return expr_write.z_vap(eos=pobj._cubic_type, A=A[p], B=B[p])

    elif pobj.is_liquid_phase():

        def Z(p):
            return expr_write.z_liq(eos=pobj._cubic_type, A=A[p], B=B[p])

    return _log_fug_coeff_method(
        A[p], b[j], bm[p], B[p], delta[p, j], Z(p), pobj._cubic_type
    )


def _log_fug_coeff_method(A, b, bm, B, delta, Z, cubic_type):
    u = EoS_param[cubic_type]["u"]
    w = EoS_param[cubic_type]["w"]
    p = sqrt(u**2 - 4 * w)

    return (
        b / bm * (Z - 1) * (B * p)
        - safe_log(Z - B, eps=eps_SL) * (B * p)
        + A
        * (b / bm - delta)
        * safe_log((2 * Z + B * (u + p)) / (2 * Z + B * (u - p)), eps=eps_SL)
    ) / (B * p)


def _d_log_fug_coeff_dT_phase_comp(blk, p, j):
    pobj = blk.params.get_phase(p)

    if not (
        pobj._mixing_rule_a == MixingRuleA.default
        and pobj._mixing_rule_b == MixingRuleB.default
    ):
        # Any
        raise NotImplementedError(
            "Block {} called for a property "
            "that is not supported by this choice "
            "of mixing rules.".format(blk.name)
        )

    cname = pobj._cubic_type.name
    am = getattr(blk, cname + "_am")[p]
    daij_dT = getattr(blk, cname + "_daij_dT")
    dam_dT = getattr(blk, cname + "_dam_dT")[p]
    b = getattr(blk, cname + "_b")[j]
    bm = getattr(blk, cname + "_bm")[p]
    A = getattr(blk, cname + "_A")[p]
    B = getattr(blk, cname + "_B")[p]
    delta = getattr(blk, cname + "_delta")[p, j]

    Z = blk.compress_fact_phase[p]
    dZ_dT = _dZ_dT(blk, p)
    T = blk.temperature

    u = EoS_param[pobj._cubic_type]["u"]
    w = EoS_param[pobj._cubic_type]["w"]
    EoS_p = sqrt(u**2 - 4 * w)

    expr = (
        A
        / (EoS_p * B)
        * (
            b / bm * (1 / am * dam_dT - 1 / T)
            + delta / T
            - 2
            / am
            * sum(
                blk.mole_frac_phase_comp[p, i] * daij_dT[i, j]
                for i in blk.component_list
            )
        )
    )

    return (
        b / bm * dZ_dT
        - (dZ_dT + B / T) / (Z - B)
        - A * (b / bm - delta) * (Z / T + dZ_dT) / (Z**2 + u * B * Z + w * B**2)
        + log((2 * Z + B * (u + EoS_p)) / (2 * Z + B * (u - EoS_p))) * expr
    )


def _bubble_dew_log_fug_coeff_method(blk, p, j, pp, pt_var):
    pobj = blk.params.get_phase(p)
    ctype = pobj._cubic_type
    cname = pobj.config.equation_of_state_options["type"].name

    # Ditch the m.fs.unit.control_volume...
    short_name = pt_var.name.split(".")[-1]
    # import pdb; pdb.set_trace()
    if short_name.startswith("temperature"):
        abbrv = "t"
        T = pt_var[pp]
        P = blk.pressure
    elif short_name.startswith("pressure"):
        abbrv = "p"
        P = pt_var[pp]
        T = blk.temperature
    else:
        raise BurntToast(
            "Users shouldn't be calling this function. "
            "If you're a dev, you know what you did."
        )

    if short_name.endswith("bubble"):
        abbrv += "bub"
        if pobj.is_liquid_phase():
            x = blk.mole_frac_comp
            log_mole_frac = blk.log_mole_frac_comp
            xidx = ()
        elif pobj.is_vapor_phase():
            x = getattr(blk, "_mole_frac_" + abbrv)
            log_mole_frac = getattr(blk, "log_mole_frac_" + abbrv)
            xidx = pp
        else:
            raise BurntToast(
                "Users shouldn't be calling this function. "
                "If you're a dev, you know what you did."
            )

    elif short_name.endswith("dew"):
        abbrv += "dew"
        if pobj.is_vapor_phase():
            x = blk.mole_frac_comp
            log_mole_frac = blk.log_mole_frac_comp
            xidx = ()
        elif pobj.is_liquid_phase():
            x = getattr(blk, "_mole_frac_" + abbrv)
            log_mole_frac = getattr(blk, "log_mole_frac_" + abbrv)
            xidx = pp
        else:
            raise BurntToast(
                "Users shouldn't be calling this function. "
                "If you're a dev, you know what you did."
            )
    else:
        raise BurntToast(
            "Users shouldn't be calling this function. "
            "If you're a dev, you know what you did."
        )

    def a(k):
        cobj = blk.params.get_component(k)
        fw = getattr(blk, cname + "_fw")[k]
        ac = getattr(blk, cname + "_a_crit")[k]
        func_alpha = getattr(blk.params, cname + "_func_alpha")

        return ac * func_alpha(T, fw, cobj)

    kappa = getattr(blk.params, cname + "_kappa")
    am = sum(
        sum(
            x[xidx, i] * x[xidx, j] * sqrt(a(i) * a(j)) * (1 - kappa[i, j])
            for j in blk.component_list
        )
        for i in blk.component_list
    )

    b = getattr(blk, cname + "_b")
    bm = sum(x[xidx, i] * b[i] for i in blk.component_list)
    R = Cubic.gas_constant(blk)

    A = am * P / (R * T) ** 2
    B = bm * P / (R * T)

    delta = (
        2
        * sqrt(a(j))
        / am
        * sum(x[xidx, i] * sqrt(a(i)) * (1 - kappa[j, i]) for i in blk.component_list)
    )

    expr_write = CubicThermoExpressions(blk)
    if pobj.is_vapor_phase():
        Z = expr_write.z_vap(eos=pobj._cubic_type, A=A, B=B)
    elif pobj.is_liquid_phase():
        Z = expr_write.z_liq(eos=pobj._cubic_type, A=A, B=B)

    return (
        _log_fug_coeff_method(A, b[j], bm, B, delta, Z, ctype)
        + log_mole_frac[xidx, j]
        + log(P / blk.params.pressure_ref)
    )


# -----------------------------------------------------------------------------
# Default rules for cubic expressions
def func_fw_PR(cobj):
    return 0.37464 + 1.54226 * cobj.omega - 0.26992 * cobj.omega**2


def func_fw_SRK(cobj):
    return 0.48 + 1.574 * cobj.omega - 0.176 * cobj.omega**2


def func_alpha_soave(T, fw, cobj):
    Tc = cobj.temperature_crit
    Tr = T / Tc
    return (1 + fw * (1 - sqrt(Tr))) ** 2


def func_dalpha_dT_soave(T, fw, cobj):
    Tc = cobj.temperature_crit
    Tr = T / Tc
    return 1 / Tc * (-fw / sqrt(Tr)) * (1 + fw * (1 - sqrt(Tr)))


def func_d2alpha_dT2_soave(T, fw, cobj):
    Tc = cobj.temperature_crit
    Tr = T / Tc
    return 1 / Tc**2 * ((fw**2 + fw) / (2 * Tr * sqrt(Tr)))


# -----------------------------------------------------------------------------
# Mixing rules
def rule_am_default(m, cname, a, p, pp=()):
    k = getattr(m.params, cname + "_kappa")
    return sum(
        sum(
            m.mole_frac_phase_comp[p, i]
            * m.mole_frac_phase_comp[p, j]
            * sqrt(a[pp, i] * a[pp, j])
            * (1 - k[i, j])
            for j in m.components_in_phase(p)
        )
        for i in m.components_in_phase(p)
    )


def rule_bm_default(m, b, p):
    return sum(m.mole_frac_phase_comp[p, i] * b[i] for i in m.components_in_phase(p))
