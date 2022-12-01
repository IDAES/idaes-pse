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
Tests for eNRTL methods

Reference:

[1] Song, Y. and Chen, C.-C., Symmetric Electrolyte Nonrandom Two-Liquid Activity
Coefficient Model, Ind. Eng. Chem. Res., 2009, Vol. 48, pgs. 7788–7797

Figures digitized using WebPlotDigitizer, https://apps.automeris.io/wpd/,
May 2021

Author: Andrew Lee
"""
import pytest
from math import log

from pyomo.environ import ConcreteModel, units as pyunits, value

from idaes.core import AqueousPhase, Solvent, Apparent, Anion, Cation
from idaes.models.properties.modular_properties.eos.enrtl import ENRTL
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
    StateIndex,
)
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.pure.electrolyte import (
    relative_permittivity_constant,
)


def rho_H2O(b, *args, **kwargs):
    return 1000 / 18e-3 * pyunits.mol / pyunits.m**3


def rho_MeOH(b, *args, **kwargs):
    return 792 / 32e-3 * pyunits.mol / pyunits.m**3


def rho_EtOH(b, *args, **kwargs):
    return 789.45 / 46e-3 * pyunits.mol / pyunits.m**3


configuration = {
    "components": {
        "H2O": {
            "type": Solvent,
            "dens_mol_liq_comp": rho_H2O,
            "relative_permittivity_liq_comp": relative_permittivity_constant,
            "parameter_data": {
                "mw": (18e-3, pyunits.kg / pyunits.mol),
                "relative_permittivity_liq_comp": 78.54,
            },
        },
        "MeOH": {
            "type": Solvent,
            "dens_mol_liq_comp": rho_MeOH,
            "relative_permittivity_liq_comp": relative_permittivity_constant,
            "parameter_data": {
                "mw": (32e-3, pyunits.kg / pyunits.mol),
                "relative_permittivity_liq_comp": 32.6146,
            },
        },
        "EtOH": {
            "type": Solvent,
            "dens_mol_liq_comp": rho_EtOH,
            "relative_permittivity_liq_comp": relative_permittivity_constant,
            "parameter_data": {
                "mw": (46e-3, pyunits.kg / pyunits.mol),
                "relative_permittivity_liq_comp": 24.113,
            },
        },
        "NaBr": {"type": Apparent, "dissociation_species": {"Na+": 1, "Br-": 1}},
        "Na+": {"type": Cation, "charge": +1},
        "Br-": {"type": Anion, "charge": -1},
    },
    "phases": {"Liq": {"type": AqueousPhase, "equation_of_state": ENRTL}},
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "state_definition": FTPx,
    "state_components": StateIndex.true,
    "pressure_ref": 1e5,
    "temperature_ref": 300,
    "parameter_data": {
        "Liq_alpha": {
            ("H2O", "EtOH"): 0.3031,
            ("H2O", "MeOH"): 0.2994,
            ("MeOH", "EtOH"): 0.3356,
            ("H2O", "Na+, Br-"): 0.2,
            ("MeOH", "Na+, Br-"): 0.2,
            ("EtOH", "Na+, Br-"): 0.1,
        },
        "Liq_tau": {
            ("H2O", "MeOH"): 1.4265,
            ("MeOH", "H2O"): -0.42864,
            ("H2O", "EtOH"): 2.2485,
            ("EtOH", "H2O"): -0.18514,
            ("MeOH", "EtOH"): -0.04394,
            ("EtOH", "MeOH"): 0.02147,
            ("H2O", "Na+, Br-"): 9.527,
            ("Na+, Br-", "H2O"): -4.790,
            ("MeOH", "Na+, Br-"): 5.910,
            ("Na+, Br-", "MeOH"): -3.863,
            ("EtOH", "Na+, Br-"): 6.118,
            ("Na+, Br-", "EtOH"): -4.450,
        },
    },
}


class Test_H2O_MeOH_EtOH(object):
    # Test case for having parameters for a second salt with 0 concentration
    # Results should be the same as for the single salt case
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.params = GenericParameterBlock(**configuration)

        m.state = m.params.build_state_block([1])

        # Need to set a value of T for checking expressions later
        m.state[1].temperature.set_value(298.15)

        return m

    @pytest.mark.unit
    def test_H2O_EtOH(self, model):
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 5 [1]
        # Form %H2O in solvent by mass, %NaBr in mix by mass
        data = {
            1.90023164606316: 4.07279483333809,
            2.52553840508245: 4.60341051665865,
            3.07268181922432: 5.10150067253983,
            3.6198252333662: 5.70745038179749,
            4.3232953372629: 6.59686802660166,
            5.18309213091442: 7.50180397873668,
            6.35554230407559: 8.65548586720136,
            10.420036237701: 11.7202417243546,
            12.2959565147588: 12.9721116679189,
            14.0155501020619: 14.2110461502489,
            15.9696337239972: 15.5348240281863,
            18.0018806908098: 16.9196471596141,
            20.1122910024999: 18.3011439277624,
            21.9658788953071: 19.4438459565312,
            23.4174838716019: 20.2934710794969,
            25.0365817297768: 20.9436527661739,
            26.9906653517121: 21.8041782373378,
            29.1792390082796: 23.0611962099913,
            30.7425059058279: 23.7038456814246,
            32.8529162175179: 24.7168170838022,
            35.1196532189629: 25.7943190088345,
            37.230063530653: 26.6777599200545,
            39.4968005320979: 27.3055144261466,
            40.9037407398913: 27.8418581569274,
            42.8578243618266: 28.7814740254471,
            45.1245613632715: 29.6986205793095,
            47.2349716749616: 30.2248062217835,
            53.01905919589: 32.5290639405499,
            55.2076328524575: 33.3573486018869,
            57.3180431641476: 33.9680033246278,
            62.8520079814683: 36.2484067677242,
            65.0562143070113: 36.6504595612375,
            66.9929282967516: 37.5059452940427,
            73.0288754845073: 39.062432020278,
            75.1392857961974: 39.8922616615831,
            76.780716038623: 40.8064682302478,
            83.7372537327126: 42.7720496832428,
            85.4568473200156: 43.2448153854284,
            87.1764409073187: 43.3361172522717,
            88.1144010458476: 43.9206762302395,
            92.6478750487374: 44.9451625065632,
            94.3674686360405: 45.3523294496401,
            96.0870622233435: 45.6278737343614,
            97.8066558106466: 45.9958029050737,
        }

        for x, g in data.items():
            n_salt = g / 102.9e-3  # MW of NaBr
            n_H2O = (100 - g) * (x / 100) / 18e-3
            n_EtOH = (100 - g) * ((100 - x) / 100) / 46e-3
            n_total = n_H2O + n_EtOH + 2 * n_salt

            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(n_H2O / n_total)
            model.state[1].mole_frac_phase_comp["Liq", "EtOH"].set_value(
                n_EtOH / n_total
            )
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(
                n_salt / n_total
            )
            model.state[1].mole_frac_phase_comp["Liq", "Br-"].set_value(
                n_salt / n_total
            )

            # Check NaBr solubility product
            # Note the Ksp given in [1] is actually ln(Ksp)
            ln_Ksp = value(
                model.state[1].Liq_log_gamma["Na+"]
                + log(value(model.state[1].mole_frac_phase_comp["Liq", "Na+"]))
                + model.state[1].Liq_log_gamma["Br-"]
                + log(value(model.state[1].mole_frac_phase_comp["Liq", "Br-"]))
            )
            assert pytest.approx(-7.157, rel=2e-2) == ln_Ksp

    @pytest.mark.unit
    def test_H2O_MeOH(self, model):
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 5 [1]
        # Form %H2O in solvent by mass, %NaBr in mix by mass
        data = {
            6.82452237334007: 19.1520757612024,
            8.59622485722806: 19.9360728637851,
            10.2637095479462: 20.7022143982515,
            11.9833031352492: 21.984026817603,
            13.7028967225522: 22.8450059962131,
            15.4224903098553: 23.6114376640682,
            17.1420838971583: 24.3195929462338,
            18.8616774844614: 24.9062851989822,
            20.5812710717644: 25.4301958788082,
            22.3008646590675: 26.3135810224415,
            24.0204582463705: 26.9298015887213,
            25.7400518336735: 27.5604530219357,
            27.4596454209766: 28.261495638838,
            29.1792390082796: 28.7765060560814,
            30.8988325955827: 29.4697504108841,
            32.6184261828857: 30.1301893048218,
            34.3380197701888: 30.7447999546189,
            36.0576133574918: 31.3102035876584,
            37.7772069447949: 31.91742945333,
            39.4968005320979: 32.732522595915,
            41.2163941194009: 32.803533475634,
            42.935987706704: 33.4957645464418,
            44.655581294007: 33.9441655334417,
            46.3751748813101: 34.4833831285068,
            48.0947684686131: 34.9794442095881,
            49.8143620559162: 36.0640561473066,
            51.5339556432192: 35.8727337255362,
            53.2535492305222: 36.3672811244236,
            54.9731428178253: 36.8359756714104,
            56.6927364051283: 37.3474811618958,
            58.4123299924314: 37.6614971193476,
            60.1319235797344: 38.6004592649369,
            61.8515171670375: 38.5624551461233,
            63.5711107543405: 38.9441861952893,
            65.2907043416436: 39.2909739234897,
            67.0102979289466: 39.6603783446311,
            68.7298915162496: 40.0924508332419,
            70.4494851035527: 40.890132583986,
            72.1690786908557: 40.6091546532659,
            73.6541822435265: 40.9707684873752,
            81.8613334556547: 42.0336203709758,
            83.5809270429577: 42.4706249979091,
            85.3005206302608: 42.7855500518895,
            87.0201142175638: 43.2516691324203,
            88.7397078048669: 43.78176630887,
            92.178894979473: 44.5269393858065,
            93.898488566776: 44.7467835048108,
            95.618082154079: 45.0120297029858,
            97.3376757413821: 45.4128499806034,
            98.9009426389303: 45.5991226308649,
        }

        for x, g in data.items():
            n_salt = g / 102.9e-3  # MW of NaBr
            n_H2O = (100 - g) * (x / 100) / 18e-3
            n_MeOH = (100 - g) * ((100 - x) / 100) / 32e-3
            n_total = n_H2O + n_MeOH + 2 * n_salt

            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(n_H2O / n_total)
            model.state[1].mole_frac_phase_comp["Liq", "MeOH"].set_value(
                n_MeOH / n_total
            )
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(
                n_salt / n_total
            )
            model.state[1].mole_frac_phase_comp["Liq", "Br-"].set_value(
                n_salt / n_total
            )

            # Check NaBr solubility product
            # Note the Ksp given in [1] is actually ln(Ksp)
            ln_Ksp = value(
                model.state[1].Liq_log_gamma["Na+"]
                + log(value(model.state[1].mole_frac_phase_comp["Liq", "Na+"]))
                + model.state[1].Liq_log_gamma["Br-"]
                + log(value(model.state[1].mole_frac_phase_comp["Liq", "Br-"]))
            )
            assert pytest.approx(-7.157, rel=2.5e-2) == ln_Ksp

    @pytest.mark.unit
    def test_MeOH_EtOH(self, model):
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 5 [1]
        # Form %MeOH in solvent by mass, %NaBr in mix by mass
        data = {
            2.71536367121331: 2.73554829623076,
            4.28421366482421: 2.86945960942034,
            6.08755369306734: 3.02744863943992,
            7.84064585674641: 3.18919593486594,
            13.1036444118254: 3.73345296717809,
            15.2214988515991: 3.93266625747611,
            17.3765739317906: 4.14197504888717,
            23.0824981078416: 4.76850518869678,
            25.2264069959077: 4.97587135803405,
            26.9125020068347: 5.17220110227576,
            28.1631155248733: 5.35033381027391,
            33.2883977104064: 5.88583802085873,
            35.6146877365198: 6.19031092340952,
            37.7325421762935: 6.44369463450957,
            43.293305854715: 7.09823866325609,
            45.5488766640346: 7.40118664900166,
            47.4381963716429: 7.60676829681603,
            53.2088844620209: 8.35096922885007,
            55.3974581185884: 8.66033648961119,
            57.7423584649107: 8.94647646211022,
            64.8217242723791: 9.94790728412982,
            67.4792779982111: 10.3476934868586,
            73.0623740608833: 11.102759854227,
            74.7484690718103: 11.3459914299202,
            76.6578879252442: 11.5676617376356,
            83.0672822051919: 12.4888362818425,
            84.7533772161189: 12.7303712917151,
            86.9084522963104: 13.0824620221641,
            93.3178465762581: 13.9936939870009,
            95.5734173855777: 14.3668192656925,
            97.3376757413821: 14.6361763117594,
        }

        for x, g in data.items():
            n_salt = g / 102.9e-3  # MW of NaBr
            n_MeOH = (100 - g) * (x / 100) / 32e-3
            n_EtOH = (100 - g) * ((100 - x) / 100) / 46e-3
            n_total = n_MeOH + n_EtOH + 2 * n_salt

            model.state[1].mole_frac_phase_comp["Liq", "MeOH"].set_value(
                n_MeOH / n_total
            )
            model.state[1].mole_frac_phase_comp["Liq", "EtOH"].set_value(
                n_EtOH / n_total
            )
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(
                n_salt / n_total
            )
            model.state[1].mole_frac_phase_comp["Liq", "Br-"].set_value(
                n_salt / n_total
            )

            # Check NaBr solubility product
            # Note the Ksp given in [1] is actually ln(Ksp)
            ln_Ksp = value(
                model.state[1].Liq_log_gamma["Na+"]
                + log(value(model.state[1].mole_frac_phase_comp["Liq", "Na+"]))
                + model.state[1].Liq_log_gamma["Br-"]
                + log(value(model.state[1].mole_frac_phase_comp["Liq", "Br-"]))
            )
            assert pytest.approx(-7.157, rel=2.3e-2) == ln_Ksp
