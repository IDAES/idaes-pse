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
Tests for eNRTL example with aqueous NaCl and KCL

Reference:

[1] Local Composition Model for Excess Gibbs Energy of Electrolyte Systems, Pt 1.
Chen, C.-C., Britt, H.I., Boston, J.F., Evans, L.B.,
AIChE Journal, 1982, Vol. 28(4), pgs. 588-596

[2] Song, Y. and Chen, C.-C., Symmetric Electrolyte Nonrandom Two-Liquid Activity
Coefficient Model, Ind. Eng. Chem. Res., 2009, Vol. 48, pgs. 7788–7797

[3] New Data on Activity Coefficients of Potassium, Nitrate, and Chloride Ions
in Aqueous Solutions of KNO3 and KCl by Ion Selective Electrodes
Dash, D., Kumar, S., Mallika, C., Kamachi Mudali, U.,
ISRN Chemical Engineering, 2012, doi:10.5402/2012/730154

Figures digitized using WebPlotDigitizer, https://apps.automeris.io/wpd/,
May 2021

Author: Andrew Lee
"""
import pytest

from pyomo.environ import ConcreteModel, value

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.examples.enrtl_H2O_NaCl_KCl import (
    configuration,
)


class TestSymmetric_0KCl:
    # Test case for having parameters for a second salt with 0 concentration
    # Results should be the same as for the single salt case
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.params = GenericParameterBlock(**configuration)

        m.state = m.params.build_state_block([1])

        # Need to set a value of T for checking expressions later
        m.state[1].temperature.set_value(298.15)

        # Set parameters to those used in 1982 paper
        m.params.Liq.tau["H2O", "Na+, Cl-"].set_value(8.885)
        m.params.Liq.tau["Na+, Cl-", "H2O"].set_value(-4.549)

        return m

    @pytest.mark.unit
    def test_parameters(self, model):
        assert model.params.Liq.tau["H2O", "Na+, Cl-"].value == 8.885
        assert model.params.Liq.tau["Na+, Cl-", "H2O"].value == -4.549

    @pytest.mark.unit
    def test_log_gamma(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 6 [1]
        data = {
            0.01620941361168171: 1.4776225827856462,
            0.03630284789262042: 1.4722881224536513,
            0.05640439199500977: 1.4714951621340302,
            0.07650593609739911: 1.4707022018144094,
            0.09660453117380643: 1.468257786944834,
            0.11669575368525856: 1.4616847357003733,
            0.13678697619671068: 1.4551116844559127,
            0.1568759869386763: 1.4473000422989863,
            0.17697310750209252: 1.4440299001544339,
            0.19706727903952678: 1.4391083034599275,
            0.21714891721653726: 1.4271680249281156,
            0.2372312926500432: 1.415640610033792,
            0.2573195661355132: 1.4074161042393771,
            0.27739235723457745: 1.3905214620577024,
            0.29747473266808344: 1.378994047163379,
            0.3175489982801387: 1.3629251322566815,
            0.3376203148662119: 1.3452047628000297,
            0.3576997412737358: 1.332025893355752,
            0.3777658970643405: 1.3114154784366803,
            0.3978268920594765: 1.2879150180551888,
            0.41788714979811703: 1.2640016940362084,
            0.43794077222829786: 1.2363725972798314,
            0.4579907083760012: 1.2066791823360115,
            0.47803400921524486: 1.1732699946547944,
            0.49807952182397497: 1.1410993978860429,
            0.5181139755852725: 1.102735846554963,
            0.5381454803205878: 1.0627208406739292,
            0.5581718242604347: 1.0198157893304751,
            0.5781907956353264: 0.9727821016121361,
            0.5982009199322718: 0.9207940502439338,
            0.618197773612298: 0.8613744534009378,
            0.6381975763183062: 0.8036063111078957,
            0.6581877946898728: 0.7404709415275024,
            0.6781794875744304: 0.678161299222086,
            0.6981645451505283: 0.6121358841792732,
            0.718139281135689: 0.5403303782116202,
            0.7381140171208497: 0.46852487224396766,
            0.7580821177975509: 0.39300359353891734,
            0.7780472694482699: 0.3158308602839135,
            0.7980079975600158: 0.2361809452039778,
            0.8209839409036424: 0.13934585427154245,
            0.8379095478581289: 0.06573379683191538,
            0.8569581176584067: -0.013013060958402267,
            0.8772461460653336: -0.10618700754926813,
            0.8913177190714701: -0.15982445783107968,
            0.9086308771904323: -0.24383529059164832,
            0.928588656276196: -0.3251366602215384,
            0.9486120511900609: -0.3696931661149465,
            0.9686917233497501: -0.3827344143467273,
            0.9855967963455463: -0.3541317749389892,
            # Error gets too large at this point
            # 0.9946128753251389: -0.19481723605693535
        }

        # Need to correct for different reference state
        OFFSET = 1.462

        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(x)
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value((1 - x) / 2)
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value((1 - x) / 2)

            assert pytest.approx(g - OFFSET, rel=4e-2, abs=2e-2) == value(
                model.state[1].Liq_log_gamma["Na+"]
            )
            assert pytest.approx(g - OFFSET, rel=4e-2, abs=2e-2) == value(
                model.state[1].Liq_log_gamma["Cl-"]
            )
