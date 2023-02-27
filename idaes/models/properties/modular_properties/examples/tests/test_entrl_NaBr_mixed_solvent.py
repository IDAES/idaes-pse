#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for eNRTL example of NaBr in mixed solvents

Reference:

[1] Song, Y. and Chen, C.-C., Symmetric Electrolyte Nonrandom Two-Liquid Activity
Coefficient Model, Ind. Eng. Chem. Res., 2009, Vol. 48, pgs. 7788–7797

Figures digitized using WebPlotDigitizer, https://apps.automeris.io/wpd/,
May 2021

Author: Andrew Lee
"""
import pytest
from math import log

from pyomo.environ import ConcreteModel, value

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.examples.enrtl_NaBr_mixed_solvent import (
    configuration,
)


@pytest.mark.unit
def model():
    m = ConcreteModel()
    m.params = GenericParameterBlock(**configuration)

    m.state = m.params.build_state_block([1])

    # Need to set a value of T for checking expressions later
    m.state[1].temperature.set_value(298.15)

    # Using 0 results in division by zero errors
    for k in m.state[1].mole_frac_phase_comp:
        m.state[1].mole_frac_phase_comp[k].set_value(1e-12)

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

        m.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(n_H2O / n_total)
        m.state[1].mole_frac_phase_comp["Liq", "EtOH"].set_value(n_EtOH / n_total)
        m.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(n_salt / n_total)
        m.state[1].mole_frac_phase_comp["Liq", "Br-"].set_value(n_salt / n_total)

        # Check NaBr solubility product
        # Note the Ksp given in [1] is actually ln(Ksp)
        ln_Ksp = value(
            m.state[1].Liq_log_gamma["Na+"]
            + log(value(m.state[1].mole_frac_phase_comp["Liq", "Na+"]))
            + m.state[1].Liq_log_gamma["Br-"]
            + log(value(m.state[1].mole_frac_phase_comp["Liq", "Br-"]))
        )
        assert pytest.approx(-7.157, rel=2e-2) == ln_Ksp
