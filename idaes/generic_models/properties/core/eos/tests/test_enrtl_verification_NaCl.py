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
Tests for eNRTL methods

References:

[1] Local Composition Model for Excess Gibbs Energy of Electrolyte Systems, Pt 1.
Chen, C.-C., Britt, H.I., Boston, J.F., Evans, L.B.,
AIChE Journal, 1982, Vol. 28(4), pgs. 588-596

Author: Andrew Lee
"""
import pytest
from math import exp

from pyomo.environ import (ConcreteModel,
                           units as pyunits,
                           value)

from idaes.core import (AqueousPhase,
                        Solvent,
                        Apparent,
                        Anion,
                        Cation)
from idaes.generic_models.properties.core.eos.enrtl import ENRTL
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock, StateIndex)
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.pure.electrolyte import \
    relative_permittivity_constant


def dummy_method(b, *args, **kwargs):
    return 1000/18e-3*pyunits.mol/pyunits.m**3


configuration = {
    "components": {
        "H2O": {"type": Solvent,
                "dens_mol_liq_comp": dummy_method,
                "relative_permittivity_liq_comp":
                    relative_permittivity_constant,
                "parameter_data": {
                    "mw": (18E-3, pyunits.kg/pyunits.mol),
                    "relative_permittivity_liq_comp": 73.41964622627609}},
        "NaCl": {"type": Apparent,
                 "dissociation_species": {"Na+": 1, "Cl-": 1}},
        "Na+": {"type": Cation,
                "charge": +1},
        "Cl-": {"type": Anion,
                "charge": -1}},
    "phases": {
        "Liq": {"type": AqueousPhase,
                "equation_of_state": ENRTL}},
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "state_definition": FTPx,
    "state_components": StateIndex.true,
    "pressure_ref": 1e5,
    "temperature_ref": 300,
    "parameter_data": {
        "Liq_tau": {
            ("H2O", "Na+, Cl-"): 8.885,  # Table 1, [1]
            ("Na+, Cl-", "H2O"): -4.549}}}  # Table 1, [1]


class TestStateBlockSymmetric(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.params = GenericParameterBlock(default=configuration)

        m.state = m.params.build_state_block([1])

        # Need to set a value of T for checking expressions later
        m.state[1].temperature.set_value(313)

        return m

    @pytest.mark.unit
    def test_parameters(self, model):
        assert model.params.Liq.tau["H2O", "Na+, Cl-"].value == 8.885
        assert model.params.Liq.tau["Na+, Cl-", "H2O"].value == -4.549

    @pytest.mark.unit
    def test_log_gamma_h2o(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 6 [1]
        data = {0.0166772823326004: -1.32412227463648,
                0.0368235348413794: -1.29987852737794,
                0.0569697873501583: -1.2756347801194,
                0.0771212007755797: -1.24850091954107,
                0.0972755632962254: -1.21971556563714,
                0.117432874912095: -1.18927871840762,
                0.137590186527965: -1.1588418711781,
                0.157748235417641: -1.12799215061718,
                0.177912182497766: -1.09383944340507,
                0.198075392304084: -1.06009960952436,
                0.218237127562791: -1.02718552230644,
                0.238402549190528: -0.992207068431534,
                0.255813423719564: -0.964639854693858,
                0.315405637317823: -0.854261500929733,
                0.335560737112275: -0.825063273694407,
                0.355727633287624: -0.789259073156702,
                0.375898215832003: -0.751390505962001,
                0.396062900185934: -0.716824925418492,
                0.416232008182701: -0.67978210488659,
                0.436393743441407: -0.646868017668675,
                0.45655990234295: -0.611476690462368,
                0.476725323970687: -0.576498236587459,
                0.496889271050811: -0.542345529375349,
                0.517057641773772: -0.505715582174845,
                0.537221588853897: -0.471562874962733,
                0.557387747755439: -0.436171547756426,
                0.57755243210937: -0.401605967212917,
                0.597715641915689: -0.367866133332205,
                0.617878114448201: -0.334539172782892,
                0.638036900611684: -0.303276578890572,
                0.658196424048972: -0.271601111666854,
                0.678351523843423: -0.242402884431529,
                0.698502199995039: -0.215681897184597,
                0.718656562515684: -0.186896543280671,
                0.7388072386673: -0.160175556033739,
                0.758950542080854: -0.137583302100794,
                0.779093845494409: -0.11499104816785,
                0.799235674360352: -0.0932245408977037,
                0.820938311873743: -0.0746252682763595,
                0.83951122208037: -0.0542331330027968,
                0.887081111079702: -0.0241043023855885,
                0.907200084457657: -0.0151368683888049,
                0.927323481478448: -0.00369219440362834,
                0.947432870296923: -0.0000921137150293738,
                0.967534886377337: -0.0006207663404183,
                0.982152912433448: -0.00191356230614259,
                0.999998: 0}

        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(x)
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(
                (1-x)/2)
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(
                (1-x)/2)

            assert pytest.approx(exp(g), rel=2e-2) == exp(value(
                model.state[1].Liq_log_gamma["H2O"]))

    @pytest.mark.unit
    def test_log_gamma_pdh(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 6 [1]
        data = {0.0102202097143692: -0.934567107147103,
                0.030381425381713: -0.931468127697603,
                0.0505354137471605: -0.932910302939649,
                0.0706887450851628: -0.934765310426382,
                0.0908453615603908: -0.934556156689684,
                0.110996721816058: -0.937649660910475,
                0.131153995318731: -0.93702767492909,
                0.151309297739069: -0.937644185681764,
                0.171467228269187: -0.936609367455696,
                0.191622530689525: -0.937225878208369,
                0.211781775274533: -0.935365395492927,
                0.231939705804651: -0.934330577266858,
                0.252104863636666: -0.92875460434924,
                0.321728090914091: -0.933511483051937,
                0.34189390577355: -0.927522677889634,
                0.36205840657812: -0.922359537216703,
                0.382220936300354: -0.91843489327783,
                0.402379523857917: -0.916987242807075,
                0.422542053580151: -0.913062598868203,
                0.442695384918153: -0.914917606354935,
                0.462852001393382: -0.914708452618238,
                0.48301124597839: -0.912847969902795,
                0.503172461645734: -0.909748990453295,
                0.523331049203297: -0.908301339982539,
                0.543494892980421: -0.903551031554295,
                0.563651509455649: -0.903341877817597,
                0.583805497821097: -0.904784053059643,
                0.603966713488441: -0.901685073610144,
                0.624131871320455: -0.896109100692527,
                0.644295058070134: -0.891771624508968,
                0.664454302655143: -0.889911141793526,
                0.684614204267596: -0.887637826833398,
                0.704786589401507: -0.877520699224235,
                0.724946491013961: -0.875247384264106,
                0.74511099181853: -0.870084243591175,
                0.765274835595655: -0.865333935162931,
                0.785453133976572: -0.851501317351594,
                0.805626833165373: -0.840558525253058,
                0.827641523255411: -0.824211134895509,
                0.84781883609516: -0.8109977654512,
                0.868006989887754: -0.790972663969571,
                0.888198428817574: -0.768883401264512,
                0.90840235126885: -0.738950325910417,
                0.928616786159249: -0.702411934641345,
                0.948852902955337: -0.652250079297633,
                0.968231171867699: -0.565426968552611,
                0.981271802390295: -0.43102358298104,
                0.994543298706282: -0.295478442814278}

        # Need to correct for different reference state
        OFFSET = 0.971

        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(x)
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(
                (1-x)/2)
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(
                (1-x)/2)

            assert pytest.approx(OFFSET, abs=0.041) == value(
                model.state[1].Liq_log_gamma_pdh["Na+"] - g)
            assert pytest.approx(OFFSET, abs=0.041) == value(
                model.state[1].Liq_log_gamma_pdh["Cl-"] - g)

    @pytest.mark.unit
    def test_log_gamma_lc(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 6 [1]
        data = {0.00642910940949107: 2.44391469691802,
                0.0261925045768978: 2.42178852215264,
                0.046347839144311: 2.42099343992118,
                0.0665026808921592: 2.41978548440975,
                0.086651608805228: 2.41362304953867,
                0.106801522357427: 2.40828636122754,
                0.126949464631366: 2.40129817979652,
                0.147097406905304: 2.39430999836549,
                0.167250277374893: 2.39145054973418,
                0.187395262731442: 2.38198512862334,
                0.207538276809731: 2.3708682143926,
                0.227685726264105: 2.36346715968161,
                0.24782331932718: 2.3478086393712,
                0.267966826225035: 2.33710459842044,
                0.288102940829415: 2.32020745827011,
                0.308242505170749: 2.30620043107958,
                0.328377141316434: 2.28806467108934,
                0.348510299003425: 2.26869029125919,
                0.36864394950998: 2.24972878470901,
                0.388771193362191: 2.2253999255192,
                0.408894987477447: 2.1981809533696,
                0.429018781592703: 2.17096198122,
                0.449137154692745: 2.13920140299072,
                0.469255527792787: 2.10744082476143,
                0.489368972697179: 2.07155151373244,
                0.509481924782006: 2.03524932942348,
                0.52958649893423: 1.99192829935501,
                0.549690580266888: 1.94819439600656,
                0.569786776486507: 1.89785452017858,
                0.589878537330041: 1.84379878483086,
                0.609963398699666: 1.78396282356355,
                0.630044810332336: 1.72123674933644,
                0.650120800949792: 1.65396906902965,
                0.670193834649858: 1.58422414904303,
                0.690261940154275: 1.5103504962567,
                0.710324624643477: 1.4319352373907,
                0.730382873756595: 1.34980411900495,
                0.750435701854498: 1.26313139453952,
                0.770486065854576: 1.17439430367423,
                0.790531501659005: 1.08152848000924,
                0.807210557418041: 1.01934951372455,
                0.832434963627757: 0.879209621277547,
                0.851556165620075: 0.779608815922458,
                0.869765897827013: 0.683966446828523,
                0.88797129322178: 0.584690792870842,
                0.906174520210461: 0.483598496481289,
                0.9243842524174: 0.387956127387353,
                0.943499737702764: 0.283565991984599,
                0.962618672725083: 0.182065969541643,
                0.98081688527469: 0.0767726875283845,
                0.99266736005891: 0.0264587394852946}

        # Need to correct for different reference state
        OFFSET = -2.4439

        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(x)
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(
                (1-x)/2)
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(
                (1-x)/2)

            assert pytest.approx(OFFSET, rel=3e-2) == value(
                model.state[1].Liq_log_gamma_lc["Na+"] - g)
            assert pytest.approx(OFFSET, rel=3e-2) == value(
                model.state[1].Liq_log_gamma_lc["Cl-"] - g)

    @pytest.mark.unit
    def test_pure_water(self, model):
        # Start by setting all mole fractions to small number
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Test pure water
        model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(1)

        # Check mixing expressions
        assert value(model.state[1].Liq_X["H2O"]) == pytest.approx(1, rel=1e-8)
        assert value(model.state[1].Liq_X["Na+"]) == pytest.approx(
            1e-12, rel=1e-8)
        assert value(model.state[1].Liq_X["Cl-"]) == pytest.approx(
            1e-12, rel=1e-8)

        for v in model.state[1].Liq_Y.values():
            assert value(v) == pytest.approx(1, rel=1e-8)

        for k, v in model.state[1].Liq_alpha.items():
            if k == ("H2O", "H2O"):
                assert value(v) == 0.3
            elif k in [("Na+", "Cl-"), ("Cl-", "Na+")]:
                assert value(v) == 0
            else:
                assert value(v) == 0.2

        assert value(model.state[1].Liq_G["H2O", "H2O"]) == 1
        assert value(model.state[1].Liq_G["H2O", "Na+"]) == exp(-0.2*8.885)
        assert value(model.state[1].Liq_G["Na+", "H2O"]) == exp(-0.2*-4.549)
        assert value(model.state[1].Liq_G["H2O", "Cl-"]) == exp(-0.2*8.885)
        assert value(model.state[1].Liq_G["Cl-", "H2O"]) == exp(-0.2*-4.549)
        assert value(model.state[1].Liq_G["Na+", "Cl-"]) == 1
        assert value(model.state[1].Liq_G["Cl-", "Na+"]) == 1

        assert value(model.state[1].Liq_tau["H2O", "H2O"]) == 0
        assert value(model.state[1].Liq_tau["H2O", "Na+"]) == pytest.approx(
            8.885, rel=1e-8)
        assert value(model.state[1].Liq_tau["Na+", "H2O"]) == pytest.approx(
            -4.549, rel=1e-8)
        assert value(model.state[1].Liq_tau["H2O", "Cl-"]) == pytest.approx(
            8.885, rel=1e-8)
        assert value(model.state[1].Liq_tau["Cl-", "H2O"]) == pytest.approx(
            -4.549, rel=1e-8)
        assert value(model.state[1].Liq_tau["Na+", "Cl-"]) == 0
        assert value(model.state[1].Liq_tau["Cl-", "Na+"]) == 0

        # Check activity coefficient contributions
        assert (value(model.state[1].Liq_log_gamma_pdh["H2O"]) ==
                pytest.approx(0, abs=1e-10))
        assert (value(model.state[1].Liq_log_gamma_lc["H2O"]) ==
                pytest.approx(0, abs=1e-10))
        assert (value(model.state[1].Liq_log_gamma_pdh["H2O"] +
                      model.state[1].Liq_log_gamma_lc["H2O"]) ==
                pytest.approx(0, abs=1e-10))

    @pytest.mark.unit
    def test_pure_NaCl(self, model):
        # Test pure NaCl
        model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(1e-12)
        model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(0.5)
        model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(0.5)

        # Check mixing expressions
        assert value(model.state[1].Liq_X["H2O"]) == pytest.approx(
            1e-12, rel=1e-8)
        assert value(model.state[1].Liq_X["Na+"]) == pytest.approx(
            0.5, rel=1e-8)
        assert value(model.state[1].Liq_X["Cl-"]) == pytest.approx(
            0.5, rel=1e-8)
        assert value(model.state[1].Liq_X["H2O"]) == pytest.approx(
            value(model.state[1].Liq_X_ref["H2O"]), rel=1e-8)
        assert value(model.state[1].Liq_X["Na+"]) == pytest.approx(
            value(model.state[1].Liq_X_ref["Na+"]), rel=1e-8)
        assert value(model.state[1].Liq_X["Cl-"]) == pytest.approx(
            value(model.state[1].Liq_X_ref["Cl-"]), rel=1e-8)

        for v in model.state[1].Liq_Y.values():
            assert value(v) == pytest.approx(1, rel=1e-8)

        for k, v in model.state[1].Liq_alpha.items():
            if k == ("H2O", "H2O"):
                assert value(v) == 0.3
            elif k in [("Na+", "Cl-"), ("Cl-", "Na+")]:
                assert value(v) == 0
            else:
                assert value(v) == 0.2

        assert value(model.state[1].Liq_G["H2O", "H2O"]) == 1
        assert value(model.state[1].Liq_G["H2O", "Na+"]) == exp(-0.2*8.885)
        assert value(model.state[1].Liq_G["Na+", "H2O"]) == exp(-0.2*-4.549)
        assert value(model.state[1].Liq_G["H2O", "Cl-"]) == exp(-0.2*8.885)
        assert value(model.state[1].Liq_G["Cl-", "H2O"]) == exp(-0.2*-4.549)
        assert value(model.state[1].Liq_G["Na+", "Cl-"]) == 1
        assert value(model.state[1].Liq_G["Cl-", "Na+"]) == 1

        assert value(model.state[1].Liq_tau["H2O", "H2O"]) == 0
        assert value(model.state[1].Liq_tau["H2O", "Na+"]) == pytest.approx(
            8.885, rel=1e-8)
        assert value(model.state[1].Liq_tau["Na+", "H2O"]) == pytest.approx(
            -4.549, rel=1e-8)
        assert value(model.state[1].Liq_tau["H2O", "Cl-"]) == pytest.approx(
            8.885, rel=1e-8)
        assert value(model.state[1].Liq_tau["Cl-", "H2O"]) == pytest.approx(
            -4.549, rel=1e-8)
        assert value(model.state[1].Liq_tau["Na+", "Cl-"]) == 0
        assert value(model.state[1].Liq_tau["Cl-", "Na+"]) == 0

        assert (value(model.state[1].Liq_log_gamma_pdh["Na+"]) ==
                pytest.approx(0, abs=1e-10))

        assert (value(model.state[1].Liq_log_gamma_lc_I["Na+"]) ==
                pytest.approx(0, abs=1e-10))
        assert (value(model.state[1].Liq_log_gamma_lc_I0["Na+"]) ==
                pytest.approx(0, abs=1e-10))
        assert (value(model.state[1].Liq_log_gamma_lc_I["Na+"]) ==
                pytest.approx(value(model.state[1].Liq_log_gamma_lc_I0["Na+"]),
                              abs=1e-12))
        assert (value(model.state[1].Liq_log_gamma_lc["Na+"]) ==
                pytest.approx(0, abs=1e-10))

        assert (value(model.state[1].Liq_log_gamma["Na+"]) ==
                pytest.approx(0, abs=1e-10))
        assert (value(model.state[1].Liq_log_gamma_pdh["Na+"] +
                      model.state[1].Liq_log_gamma_lc["Na+"]) ==
                pytest.approx(model.state[1].Liq_log_gamma["Na+"], abs=1e-10))

        assert (value(model.state[1].Liq_log_gamma_pdh["Cl-"]) ==
                pytest.approx(0, abs=1e-10))

        assert (value(model.state[1].Liq_log_gamma_lc_I["Cl-"]) ==
                pytest.approx(0, abs=1e-10))
        assert (value(model.state[1].Liq_log_gamma_lc_I0["Cl-"]) ==
                pytest.approx(0, abs=1e-10))
        assert (value(model.state[1].Liq_log_gamma_lc_I["Cl-"]) ==
                pytest.approx(value(model.state[1].Liq_log_gamma_lc_I0["Cl-"]),
                              abs=1e-12))
        assert (value(model.state[1].Liq_log_gamma_lc["Cl-"]) ==
                pytest.approx(0, abs=1e-10))

        assert (value(model.state[1].Liq_log_gamma["Cl-"]) ==
                pytest.approx(0, abs=1e-10))
        assert (value(model.state[1].Liq_log_gamma_pdh["Cl-"] +
                      model.state[1].Liq_log_gamma_lc["Cl-"]) ==
                pytest.approx(model.state[1].Liq_log_gamma["Cl-"], abs=1e-10))
