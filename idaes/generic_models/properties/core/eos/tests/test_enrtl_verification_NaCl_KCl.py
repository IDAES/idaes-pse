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

Reference:

[1] Local Composition Model for Excess Gibbs Energy of Electrolyte Systems, Pt 1.
Chen, C.-C., Britt, H.I., Boston, J.F., Evans, L.B.,
AIChE Journal, 1982, Vol. 28(4), pgs. 588-596

[2] Song, Y. and Chen, C.-C., Symmetric Electrolyte Nonrandom Two-Liquid Activity
Coefficient Model, Ind. Eng. Chem. Res., 2009, Vol. 48, pgs. 7788–7797

Figures digitized using WebPlotDigitizer, https://apps.automeris.io/wpd/,
May 2021

Author: Andrew Lee
"""
import pytest
from math import exp, log

from pyomo.environ import (ConcreteModel,
                           units as pyunits,
                           value)

from idaes.core import (AqueousPhase,
                        Solvent,
                        Apparent,
                        Anion,
                        Cation)
from idaes.generic_models.properties.core.eos.enrtl import ENRTL
from idaes.generic_models.properties.core.eos.enrtl_reference_states import \
    Unsymmetric
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
                    "relative_permittivity_liq_comp": 73.41965}},
        "NaCl": {"type": Apparent,
                 "dissociation_species": {"Na+": 1, "Cl-": 1}},
        "KCl": {"type": Apparent,
                "dissociation_species": {"K+": 1, "Cl-": 1}},
        "Na+": {"type": Cation,
                "charge": +1},
        "K+": {"type": Cation,
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
            ("H2O", "Na+, Cl-"): 9.0234,
            ("Na+, Cl-", "H2O"): -4.5916,
            ("H2O", "K+, Cl-"): 8.1354,
            ("K+, Cl-", "H2O"): -4.1341}}}


class TestSymmetric_0KCl(object):
    # Test case for having parameters for a second salt with 0 concentration
    # Results should be the same as for the single salt case
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.params = GenericParameterBlock(default=configuration)

        m.state = m.params.build_state_block([1])

        # Need to set a value of T for checking expressions later
        m.state[1].temperature.set_value(313)

        # Set parameters to those used in 1982 paper
        m.params.Liq.tau["H2O", "Na+, Cl-"].set_value(8.885)
        m.params.Liq.tau["Na+, Cl-", "H2O"].set_value(-4.549)

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
    def test_log_gamma(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 6 [1]
        data = {0.01620941361168171: 1.4776225827856462,
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
                0.9946128753251389: -0.19481723605693535}

        # Need to correct for different reference state
        OFFSET = -1.4592

        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(x)
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(
                (1-x)/2)
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(
                (1-x)/2)

            assert pytest.approx(OFFSET, rel=0.051) == value(
                model.state[1].Liq_log_gamma["Na+"] - g)
            assert pytest.approx(OFFSET, rel=0.051) == value(
                model.state[1].Liq_log_gamma["Cl-"] - g)


class TestUnsymmetric_0KCl(object):
    @pytest.fixture(scope="class")
    def model(self):
        config = dict(configuration)
        eos_opt = config["phases"]["Liq"]["equation_of_state_options"] = {}
        eos_opt["reference_state"] = Unsymmetric

        m = ConcreteModel()
        m.params = GenericParameterBlock(default=config)

        m.state = m.params.build_state_block([1])

        # Need to set a value of T for checking expressions later
        m.state[1].temperature.set_value(313)

        # Set parameters to those used in 1982 paper
        m.params.Liq.tau["H2O", "Na+, Cl-"].set_value(8.885)
        m.params.Liq.tau["Na+, Cl-", "H2O"].set_value(-4.549)

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

        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(x)
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(
                (1-x)/2)
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(
                (1-x)/2)

            assert pytest.approx(g, abs=0.06) == value(
                model.state[1].Liq_log_gamma_pdh["Na+"])
            assert pytest.approx(g, abs=0.06) == value(
                model.state[1].Liq_log_gamma_pdh["Cl-"])

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

        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(x)
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(
                (1-x)/2)
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(
                (1-x)/2)

            assert pytest.approx(g, rel=3e-2, abs=2e-2) == value(
                model.state[1].Liq_log_gamma_lc["Na+"])
            assert pytest.approx(g, rel=3e-2, abs=2e-2) == value(
                model.state[1].Liq_log_gamma_lc["Cl-"])

    @pytest.mark.unit
    def test_log_gamma(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 6 [1]
        data = {0.01620941361168171: 1.4776225827856462,
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
                # Error gets too big at this point
                # 0.9946128753251389: -0.19481723605693535
                }

        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(x)
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(
                (1-x)/2)
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(
                (1-x)/2)

            assert pytest.approx(g, rel=4e-2, abs=4e-2) == value(
                model.state[1].Liq_log_gamma["Na+"])
            assert pytest.approx(g, rel=4e-2, abs=4e-2) == value(
                model.state[1].Liq_log_gamma["Cl-"])

    @pytest.mark.unit
    def test_pure_water(self, model):
        # Start by setting all mole fractions to small number
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Test pure water
        model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(1)

        # Unsymmetric reference state - all ln(gammas) should be 0
        for v in model.state[1].Liq_log_gamma.values():
            assert value(v) == pytest.approx(0, abs=1.2e-5)
        for v in model.state[1].Liq_log_gamma_pdh.values():
            assert value(v) == pytest.approx(0, abs=1.2e-5)
        for v in model.state[1].Liq_log_gamma_lc.values():
            assert value(v) == pytest.approx(0, abs=1e-5)


class TestUnsymmetric_0NaCl(object):
    @pytest.fixture(scope="class")
    def model(self):
        config = dict(configuration)
        eos_opt = config["phases"]["Liq"]["equation_of_state_options"] = {}
        eos_opt["reference_state"] = Unsymmetric

        m = ConcreteModel()
        m.params = GenericParameterBlock(default=config)
        
        # Set parameters to those used in 1982 paper
        m.params.Liq.tau["H2O", "K+, Cl-"].set_value(8.064)
        m.params.Liq.tau["K+, Cl-", "H2O"].set_value(-4.107)

        m.state = m.params.build_state_block([1])

        # Need to set a value of T for checking expressions later
        m.state[1].temperature.set_value(298)

        return m

    @pytest.mark.unit
    def test_parameters(self, model):
        assert model.params.Liq.tau["H2O", "K+, Cl-"].value == 8.064
        assert model.params.Liq.tau["K+, Cl-", "H2O"].value == -4.107

    @pytest.mark.unit
    def test_log_gamma(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data from [2] - Form {molality: gamma_KCl}
        # Ignoring final points due to large error (stll <7%)
        data = {0.001: 0.965,
                0.002: 0.951,
                0.005: 0.927,
                0.01: 0.901,
                0.0201: 0.869,
                0.05: 0.816,
                0.0982: 0.768,
                0.1983: 0.717,
                0.3049: 0.687,
                0.3999: 0.665,
                0.5001: 0.649,
                0.6032: 0.636,
                0.7003: 0.626,
                0.7996: 0.617,
                0.8968: 0.61,
                0.9926: 0.604,
                1.1899: 0.594,
                1.3844: 0.586,
                1.5817: 0.58,
                1.7923: 0.576,
                1.9895: 0.573,
                2.5039: 0.568,
                2.9837: 0.568,
                3.4982: 0.571,
                # 3.994: 0.576,
                # 4.4897: 0.584,
                # 4.7909: 0.589,
                # 4.9908: 0.593
                }

        for x, g in data.items():
            w = 1000/18

            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(
                w/(w+2*x))
            model.state[1].mole_frac_phase_comp["Liq", "K+"].set_value(
                x/(w+2*x))
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(
                x/(w+2*x))

            gamma_KCl = exp(value(0.5*(model.state[1].Liq_log_gamma["K+"] +
                                       model.state[1].Liq_log_gamma["Cl-"])))
            assert pytest.approx(g, rel=5e-2) == gamma_KCl


class TestUnsymmetric_Mixed(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        config = dict(configuration)
        eos_opt = config["phases"]["Liq"]["equation_of_state_options"] = {}
        eos_opt["reference_state"] = Unsymmetric

        m.params = GenericParameterBlock(default=config)

        m.state = m.params.build_state_block([1])

        # Need to set a value of T for checking expressions later
        m.state[1].temperature.set_value(300)

        return m

    @pytest.mark.unit
    def test_log_gamma_salt(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data regressed from digitized Figs 1 and 2 [2]
        # Form x_KCl: (ln(gamma NaCl), ln(gamma KCl))
        data = {0.00000001: (-0.105600000912, -0.143200000791),
                0.1: (-0.114916, -0.151286),
                0.2: (-0.124624, -0.159724),
                0.3: (-0.134724, -0.168514),
                0.4: (-0.145216, -0.177656),
                0.5: (-0.1561, -0.18715),
                0.6: (-0.167376, -0.196996),
                0.7: (-0.179044, -0.207194),
                0.8: (-0.191104, -0.217744),
                0.9: (-0.203556, -0.228646),
                0.99999999: (-0.216399998696, -0.239899998857)}

        for x, g in data.items():
            m = 4  # molality of solution
            w = 1000/18  # mol H2O per kg H2O
            k = x*m  # mol K+
            n = m-k  # mol Na+

            # Set composition
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(
                w/(w+2*m))
            model.state[1].mole_frac_phase_comp["Liq", "K+"].set_value(
                k/(w+2*m))
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(
                n/(w+2*m))
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(
                m/(w+2*m))

            assert pytest.approx(0.06293706294, rel=1e-8) == value(
                model.state[1].mole_frac_phase_comp["Liq", "Cl-"])
            assert pytest.approx(0.874125874, rel=1e-8) == value(
                model.state[1].mole_frac_phase_comp["Liq", "H2O"])
            assert value(model.state[1].Liq_ionic_strength) == pytest.approx(
                0.06293706294, rel=1e-8)

            # Check local contribution term
            # THis should only depend on ionic strength, and should be constant
            assert pytest.approx(-0.81465175, rel=1e-6) == value(
                model.state[1].Liq_log_gamma_pdh["Na+"])
            assert pytest.approx(-0.81465175, rel=1e-6) == value(
                model.state[1].Liq_log_gamma_pdh["K+"])
            assert pytest.approx(-0.81465175, rel=1e-6) == value(
                model.state[1].Liq_log_gamma_pdh["Cl-"])

            # Calculate molality scale activity coefficients.
            lgm_Na = log(
                value(model.state[1].mole_frac_phase_comp["Liq", "Na+"]) +
                exp(value(model.state[1].Liq_log_gamma["Na+"])))
            lgm_K = log(
                value(model.state[1].mole_frac_phase_comp["Liq", "K+"]) +
                exp(value(model.state[1].Liq_log_gamma["K+"])))
            lgm_Cl = log(
                value(model.state[1].mole_frac_phase_comp["Liq", "Cl-"]) +
                exp(value(model.state[1].Liq_log_gamma["Cl-"])))

            # Calculate mean activity coefficients for salts
            lg_NaCl = value(0.5*(model.state[1].Liq_log_gamma["Na+"] +
                                 model.state[1].Liq_log_gamma["Cl-"]))
            lg_KCl = value(0.5*(model.state[1].Liq_log_gamma["K+"] +
                                model.state[1].Liq_log_gamma["Cl-"]))
            lgm_NaCl = 0.5*(lgm_Na + lgm_Cl)
            lgm_KCl = 0.5*(lgm_K + lgm_Cl)
            # print(x, g[0], lgm_NaCl, g[1], lgm_KCl, lg_NaCl, lg_KCl)
            # model.state[1].mole_frac_phase_comp.display()
            # print(value(model.state[1].Liq_log_gamma["Na+"]),
            #       value(model.state[1].Liq_log_gamma["Cl-"]))
            
            # TODO: This is WIP
