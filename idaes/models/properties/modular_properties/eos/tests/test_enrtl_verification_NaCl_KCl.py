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
from math import log

from pyomo.environ import ConcreteModel, units as pyunits, value

from idaes.core import AqueousPhase, Solvent, Apparent, Anion, Cation
from idaes.models.properties.modular_properties.eos.enrtl import ENRTL
from idaes.models.properties.modular_properties.eos.enrtl_reference_states import (
    Unsymmetric,
)
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
    StateIndex,
)
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.pure.electrolyte import (
    relative_permittivity_constant,
)


def dummy_method(b, *args, **kwargs):
    return 1000 / 18e-3 * pyunits.mol / pyunits.m**3


configuration = {
    "components": {
        "H2O": {
            "type": Solvent,
            "dens_mol_liq_comp": dummy_method,
            "relative_permittivity_liq_comp": relative_permittivity_constant,
            "parameter_data": {
                "mw": (18e-3, pyunits.kg / pyunits.mol),
                "relative_permittivity_liq_comp": 78.54,
            },
        },
        "NaCl": {"type": Apparent, "dissociation_species": {"Na+": 1, "Cl-": 1}},
        "KCl": {"type": Apparent, "dissociation_species": {"K+": 1, "Cl-": 1}},
        "Na+": {"type": Cation, "charge": +1},
        "K+": {"type": Cation, "charge": +1},
        "Cl-": {"type": Anion, "charge": -1},
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
        "Liq_tau": {
            ("H2O", "Na+, Cl-"): 9.0234,
            ("Na+, Cl-", "H2O"): -4.5916,
            ("H2O", "K+, Cl-"): 8.1354,
            ("K+, Cl-", "H2O"): -4.1341,
        }
    },
}


class TestSymmetric_0KCl(object):
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
    def test_log_gamma_h2o(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 6 [1]
        data = {
            0.0166772823326004: -1.32412227463648,
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
            0.999998: 0,
        }

        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(x)
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value((1 - x) / 2)
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value((1 - x) / 2)

            assert pytest.approx(g, rel=4.5e-2, abs=0.01) == value(
                model.state[1].Liq_log_gamma["H2O"]
            )

    @pytest.mark.unit
    def test_log_gamma_pdh(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 6 [1]
        data = {
            0.0102202097143692: -0.934567107147103,
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
            0.994543298706282: -0.295478442814278,
        }

        # Need to correct for different reference state
        OFFSET = 0.971

        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(x)
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value((1 - x) / 2)
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value((1 - x) / 2)

            assert pytest.approx(g + OFFSET, rel=0.04, abs=0.07) == value(
                model.state[1].Liq_log_gamma_pdh["Na+"]
            )
            assert pytest.approx(g + OFFSET, rel=0.04, abs=0.07) == value(
                model.state[1].Liq_log_gamma_pdh["Cl-"]
            )

    @pytest.mark.unit
    def test_log_gamma_lc(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 6 [1]
        data = {
            0.00642910940949107: 2.44391469691802,
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
            0.99266736005891: 0.0264587394852946,
        }

        # Need to correct for different reference state
        OFFSET = 2.414

        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(x)
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value((1 - x) / 2)
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value((1 - x) / 2)

            assert pytest.approx(g - OFFSET, rel=4e-2, abs=6e-2) == value(
                model.state[1].Liq_log_gamma_lc["Na+"]
            )
            assert pytest.approx(g - OFFSET, rel=4e-2, abs=6e-2) == value(
                model.state[1].Liq_log_gamma_lc["Cl-"]
            )

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


class TestUnsymmetric_0KCl(object):
    @pytest.fixture(scope="class")
    def model(self):
        config = dict(configuration)
        eos_opt = config["phases"]["Liq"]["equation_of_state_options"] = {}
        eos_opt["reference_state"] = Unsymmetric

        m = ConcreteModel()
        m.params = GenericParameterBlock(**config)

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
    def test_log_gamma_h2o(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 6 [1]
        data = {
            0.0166772823326004: -1.32412227463648,
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
            0.999998: 0,
        }

        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(x)
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value((1 - x) / 2)
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value((1 - x) / 2)

            assert pytest.approx(g, rel=4.5e-2, abs=0.01) == value(
                model.state[1].Liq_log_gamma["H2O"]
            )

    @pytest.mark.unit
    def test_log_gamma_pdh(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 6 [1]
        data = {
            0.0102202097143692: -0.934567107147103,
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
            0.994543298706282: -0.295478442814278,
        }

        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(x)
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value((1 - x) / 2)
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value((1 - x) / 2)

            assert pytest.approx(g, rel=0.04, abs=0.07) == value(
                model.state[1].Liq_log_gamma_pdh["Na+"]
            )
            assert pytest.approx(g, rel=0.04, abs=0.07) == value(
                model.state[1].Liq_log_gamma_pdh["Cl-"]
            )

    @pytest.mark.unit
    def test_log_gamma_lc(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 6 [1]
        data = {
            0.00642910940949107: 2.44391469691802,
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
            0.99266736005891: 0.0264587394852946,
        }

        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(x)
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value((1 - x) / 2)
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value((1 - x) / 2)

            assert pytest.approx(g, rel=3e-2, abs=2e-2) == value(
                model.state[1].Liq_log_gamma_lc["Na+"]
            )
            assert pytest.approx(g, rel=3e-2, abs=2e-2) == value(
                model.state[1].Liq_log_gamma_lc["Cl-"]
            )

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
            # Error gets too big at this point
            # 0.9946128753251389: -0.19481723605693535
        }

        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(x)
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value((1 - x) / 2)
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value((1 - x) / 2)

            assert pytest.approx(g, rel=3e-2, abs=6e-2) == value(
                model.state[1].Liq_log_gamma["Na+"]
            )
            assert pytest.approx(g, rel=3e-2, abs=6e-2) == value(
                model.state[1].Liq_log_gamma["Cl-"]
            )


class TestUnsymmetric_0NaCl(object):
    @pytest.fixture(scope="class")
    def model(self):
        config = dict(configuration)
        eos_opt = config["phases"]["Liq"]["equation_of_state_options"] = {}
        eos_opt["reference_state"] = Unsymmetric

        m = ConcreteModel()
        m.params = GenericParameterBlock(**config)

        # Set parameters to those used in 1982 paper [1]
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

        # Data from [3] - Form {molality: gamma_KCl}
        data = {
            0.001: 0.965,
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
            3.994: 0.576,
            4.4897: 0.584,
            4.7909: 0.589,
            4.9908: 0.593,
        }

        w = 1000 / 18
        for x, g in data.items():
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(w / (w + 2 * x))
            model.state[1].mole_frac_phase_comp["Liq", "K+"].set_value(x / (w + 2 * x))
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(x / (w + 2 * x))

            conv = log(1 + 18 * 2 * x / 1000)  # Convert mole frac to molal basis
            assert pytest.approx(log(g), rel=3e-2, abs=6e-3) == value(
                model.state[1].Liq_log_gamma["K+"] - conv
            )
            assert pytest.approx(log(g), rel=3e-2, abs=6e-3) == value(
                model.state[1].Liq_log_gamma["Cl-"] - conv
            )


class TestUnsymmetric_Mixed_tau0(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        config = dict(configuration)
        eos_opt = config["phases"]["Liq"]["equation_of_state_options"] = {}
        eos_opt["reference_state"] = Unsymmetric

        m.params = GenericParameterBlock(**config)

        m.state = m.params.build_state_block([1])

        # Need to set a value of T for checking expressions later
        m.state[1].temperature.set_value(298.15)

        return m

    @pytest.mark.unit
    def test_binary_parameters(self, model):
        # Make sure binary parameters are correct
        assert model.params.Liq.tau["H2O", "Na+, Cl-"].value == 9.0234
        assert model.params.Liq.tau["Na+, Cl-", "H2O"].value == -4.5916
        assert model.params.Liq.tau["H2O", "K+, Cl-"].value == 8.1354
        assert model.params.Liq.tau["K+, Cl-", "H2O"].value == -4.1341

    @pytest.mark.unit
    def test_log_gamma_NaCl(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 1 [2]
        # Form x_KCl: log(gamma_m NaCl)
        # Important notes, data appears to be in log_10 and on a molal basis
        # Thus, need to convert model output to this form
        data = {
            1e-12: -0.1056,
            0.0751258278145695: -0.113328198570951,
            0.0976953642384106: -0.114607787343188,
            0.120623147272153: -0.11695157101671,
            0.142834437086092: -0.119035007793955,
            0.164568064753495: -0.121020719322441,
            0.187973509933774: -0.123402355074649,
            0.211378955114054: -0.125424751806415,
            0.232276674025018: -0.127697473900636,
            0.2560403658152: -0.129986897040709,
            0.279923473142016: -0.132425202153624,
            0.302493009565857: -0.134818431889228,
            0.324226637233259: -0.137175357072171,
            0.347632082413539: -0.139533043556349,
            0.368529801324503: -0.141913537356703,
            0.390741091138442: -0.144253573434263,
            0.413668874172185: -0.14674789536397,
            0.435402501839587: -0.149056922010853,
            0.457972038263429: -0.15155792345259,
            0.48054157468727: -0.154058924894327,
            0.501558709134868: -0.15629235614116,
            0.520330831493745: -0.158378332401811,
            0.541562913907284: -0.160888304509771,
            0.564132450331125: -0.163197711807272,
            0.586701986754966: -0.166034003001421,
            0.609510354252076: -0.168453001424376,
            0.63100515084621: -0.170915878894129,
            0.653932933879953: -0.173605216292076,
            0.674950068327551: -0.176119196107253,
            0.697041942604856: -0.17853786825825,
            0.720447387785136: -0.18123084449484,
            0.7413451066961: -0.183934653413591,
            0.758732008830022: -0.186191407214401,
            0.778125091979396: -0.188585585404461,
            0.801769368232944: -0.19140868071053,
            0.823264164827078: -0.193922878040346,
            0.838310522442972: -0.195941468265935,
            0.890805592347314: -0.202460414614303,
            0.910198675496688: -0.20481148412191,
            0.932409965310627: -0.207521023191923,
            0.950322295805739: -0.209499862617294,
        }

        m = 4  # molality of solution
        w = 1000 / 18  # mol H2O per kg H2O

        for x, g in data.items():
            k = x * m  # mol K+
            n = m - k  # mol Na+

            # Set composition
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(w / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "K+"].set_value(k / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(n / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(m / (w + 2 * m))

            assert pytest.approx(0.06293706294, rel=1e-8) == value(
                model.state[1].mole_frac_phase_comp["Liq", "Cl-"]
            )
            assert pytest.approx(0.874125874, rel=1e-8) == value(
                model.state[1].mole_frac_phase_comp["Liq", "H2O"]
            )
            assert value(model.state[1].Liq_ionic_strength) == pytest.approx(
                0.06293706294, rel=1e-8
            )

            # Calculate mean activity coefficient for NaCl
            # Natural log, mole fraction basis
            ln_pm = value(
                0.5
                * (
                    model.state[1].Liq_log_gamma["Na+"]
                    + model.state[1].Liq_log_gamma["Cl-"]
                )
            )

            # Convert to log base 10 and molal basis
            # Eqn 28 [1] - note that m=4 appears to be based on total salt
            log_m = ln_pm / log(10) - log((w + 2 * m) / w, 10)
            assert pytest.approx(g, rel=1e-3, abs=1e-3) == log_m

    @pytest.mark.unit
    def test_log_gamma_KCl(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 2 [2]
        # Form x_KCl: log(gamma_m KCl)
        # Important notes, data appears to be in log_10 and on a molal basis
        # Thus, need to convert model output to this form
        data = {
            0.0129962981042131: -0.145305327383741,
            0.030732184190065: -0.146667293238219,
            0.0507470065983459: -0.147845261727147,
            0.0745635482458493: -0.149617513015021,
            0.0963530760430225: -0.151219807408623,
            0.11799240133926: -0.15327982779779,
            0.139427357933593: -0.154892684744564,
            0.16286430422693177: -0.156754262749827,
            0.18592090645731912: -0.158512000605861,
            0.2072038463717374: -0.160114525371398,
            0.22937553241084108: -0.16228871809103,
            0.2524319330657862: -0.163994449482976,
            0.2728287473362676: -0.165701390327575,
            0.2962677957735029: -0.168105321458338,
            0.31805925293562287: -0.170205392008222,
            0.3402301326729579: -0.172171558871497,
            0.3624012139857349: -0.174189732198862,
            0.3854588240933333: -0.176207502375342,
            0.4064907449892218: -0.178397072421566,
            0.4289156673767195: -0.180556290965491,
            0.45197347905976004: -0.182626067606061,
            0.4732576284268315: -0.184540631156132,
            0.4954291128904931: -0.186662817411674,
            0.5172202820876957: -0.18868859301286,
            0.5388865609002058: -0.19116762539409,
            0.5623249182216393: -0.193393248647977,
            0.5823430010327078: -0.195412401056061,
            0.5991933977607408: -0.197042944808064,
            0.652406290871447: -0.202479434484355,
            0.6710309546417759: -0.204421210719124,
            0.6933255729041333: -0.205636941007559,
            0.7153739235690989: -0.20866558323021,
            0.7366604918414769: -0.21120422434935,
            0.7588315731542539: -0.213222397676715,
            0.7804979095597475: -0.215716289047685,
            0.8022900290411777: -0.217987237979576,
            0.823576798888998: -0.220577885562805,
            0.8453684576265602: -0.222729962576778,
            0.8670352547759218: -0.225342725865667,
            0.8892071423904676: -0.227568925049388,
            0.9120150005251731: -0.23050063828083,
            0.9317857251371712: -0.234051327390149,
            0.9512961168980646: -0.235786092975788,
            0.9672582363606549: -0.236965934458349,
        }

        m = 4  # molality of solution
        w = 1000 / 18  # mol H2O per kg H2O

        for x, g in data.items():
            k = x * m  # mol K+
            n = m - k  # mol Na+

            # Set composition
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(w / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "K+"].set_value(k / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(n / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(m / (w + 2 * m))

            assert pytest.approx(0.06293706294, rel=1e-8) == value(
                model.state[1].mole_frac_phase_comp["Liq", "Cl-"]
            )
            assert pytest.approx(0.874125874, rel=1e-8) == value(
                model.state[1].mole_frac_phase_comp["Liq", "H2O"]
            )
            assert value(model.state[1].Liq_ionic_strength) == pytest.approx(
                0.06293706294, rel=1e-8
            )

            # Calculate mean activity coefficients for KCl
            ln_pm = value(
                0.5
                * (
                    model.state[1].Liq_log_gamma["K+"]
                    + model.state[1].Liq_log_gamma["Cl-"]
                )
            )

            # Convert to log base 10 and molal basis
            # Eqn 28 [1] - note that m=4 appears to be based on total salt
            log_m = ln_pm / log(10) - log((w + 2 * m) / w, 10)
            assert pytest.approx(g, rel=1e-3, abs=2e-3) == log_m


class TestUnsymmetric_Mixed_tau1(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        config = dict(configuration)
        eos_opt = config["phases"]["Liq"]["equation_of_state_options"] = {}
        eos_opt["reference_state"] = Unsymmetric

        m.params = GenericParameterBlock(**config)

        m.state = m.params.build_state_block([1])

        # Need to set a value of T for checking expressions later
        m.state[1].temperature.set_value(298.15)

        # Set salt-salt interaction parameter
        m.params.Liq.tau["K+, Cl-", "Na+, Cl-"].set_value(-1)
        m.params.Liq.tau["Na+, Cl-", "K+, Cl-"].set_value(1)

        return m

    @pytest.mark.unit
    def test_binary_parameters(self, model):
        # Make sure binary parameters are correct
        assert model.params.Liq.tau["H2O", "Na+, Cl-"].value == 9.0234
        assert model.params.Liq.tau["Na+, Cl-", "H2O"].value == -4.5916
        assert model.params.Liq.tau["H2O", "K+, Cl-"].value == 8.1354
        assert model.params.Liq.tau["K+, Cl-", "H2O"].value == -4.1341

    @pytest.mark.unit
    def test_log_gamma_NaCl(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 1 [2]
        # Form x_KCl: log(gamma_m NaCl)
        # Important notes, data appears to be in log_10 and on a molal basis
        # Thus, need to convert model output to this form
        data = {
            1e-12: -0.1056,
            0.0308226637233259: -0.105016851831747,
            0.0542281089036056: -0.104718802398099,
            0.0776335540838852: -0.104420752964451,
            0.099844843897824: -0.103835856101046,
            0.123608535688006: -0.103393947923243,
            0.147013980868285: -0.102952202881418,
            0.170419426048565: -0.102366762231417,
            0.193824871228845: -0.102068712797769,
            0.217230316409124: -0.10177066336412,
            0.242098601913171: -0.101795262910288,
            0.257473352254809: -0.101423890651408,
            0.277248565121413: -0.101168550217065,
            0.299149374540103: -0.101158577170881,
            0.322913066330285: -0.101106699929557,
            0.345124356144224: -0.10106579358282,
            0.367335645958163: -0.101004359292059,
            0.499767476085357: -0.10135461223899,
            0.56747608535688: -0.100892692714428,
            0.589687375170819: -0.100985218003855,
            0.613451066961001: -0.101015452538631,
            0.63685651214128: -0.101004794321336,
            0.66026195732156: -0.100994136104041,
            0.68366740250184: -0.101270869103099,
            0.707072847682119: -0.101547602102156,
            0.72964238410596: -0.101537324535479,
            0.751376011773363: -0.101919324732628,
            0.768762913907285: -0.101519510086571,
            0.788155997056659: -0.102085461424947,
            0.811561442236939: -0.102074803207652,
            0.834966887417219: -0.102638927423062,
            0.858372332597498: -0.102628269205767,
            0.881777777777778: -0.103228317323222,
            0.905183222958057: -0.103397278616147,
            0.928588668138337: -0.103674011615205,
            0.970384105960265: -0.103439435672056,
        }

        m = 4  # molality of solution
        w = 1000 / 18  # mol H2O per kg H2O

        for x, g in data.items():
            k = x * m  # mol K+
            n = m - k  # mol Na+

            # Set composition
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(w / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "K+"].set_value(k / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(n / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(m / (w + 2 * m))

            assert pytest.approx(0.06293706294, rel=1e-8) == value(
                model.state[1].mole_frac_phase_comp["Liq", "Cl-"]
            )
            assert pytest.approx(0.874125874, rel=1e-8) == value(
                model.state[1].mole_frac_phase_comp["Liq", "H2O"]
            )
            assert value(model.state[1].Liq_ionic_strength) == pytest.approx(
                0.06293706294, rel=1e-8
            )

            # Calculate mean activity coefficient for NaCl
            # Natural log, mole fraction basis
            ln_pm = value(
                0.5
                * (
                    model.state[1].Liq_log_gamma["Na+"]
                    + model.state[1].Liq_log_gamma["Cl-"]
                )
            )

            # Convert to log base 10 and molal basis
            # Eqn 28 [1] - note that m=4 appears to be based on total salt
            log_m = ln_pm / log(10) - log((w + 2 * m) / w, 10)
            assert pytest.approx(g, rel=2e-3, abs=3e-3) == log_m

    @pytest.mark.unit
    def test_log_gamma_KCl(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 2 [2]
        # Form x_KCl: log(gamma_m KCl)
        # Important notes, data appears to be in log_10 and on a molal basis
        # Thus, need to convert model output to this form
        data = {
            0.0306504344206857: -0.313263578684857,
            0.0517245274885355: -0.311486988458052,
            0.0742677280058269: -0.309128617009572,
            0.0941260080120701: -0.307139252576374,
            0.118594245876905: -0.304858253566582,
            0.196254305187035: -0.296515421571863,
            0.2169991155507: -0.294755198591704,
            0.239827272254305: -0.291969046833168,
            0.258001144581447: -0.290602009755382,
            0.283488892357318: -0.288662119044999,
            0.304765620935435: -0.286683170132394,
            0.326548938289222: -0.284972048892172,
            0.349978669163935: -0.282933582719038,
            0.36948233702721: -0.281188983019767,
            0.393418656677592: -0.279288150511607,
            0.416468445970553: -0.277257123996039,
            0.439518235263514: -0.275590640701214,
            0.460794963841632: -0.273663769391572,
            0.484509650902659: -0.271541607070818,
            0.502993808854898: -0.270199567242454,
            0.528171271005671: -0.268612241904133,
            0.551221060298632: -0.267049913815235,
            0.571611258519328: -0.265331352917446,
            0.595040989394041: -0.263858300719342,
            0.616824306747828: -0.262467084754465,
            0.639874096040789: -0.261008911871493,
            0.66292388533375: -0.259550738988521,
            0.685087144269289: -0.258144643708512,
            0.706363872847407: -0.256582315619614,
            0.729413662140368: -0.254863754721825,
            0.751576921075907: -0.253582645688928,
            0.774626710368868: -0.252311952176624,
            0.797676499661829: -0.251270400117358,
            0.820726288954789: -0.250020537646239,
            0.842509606308577: -0.248592123393532,
            0.86593933718329: -0.247312502292148,
            0.887216065761407: -0.246375105438809,
            0.907960876125072: -0.245146074008875,
        }

        m = 4  # molality of solution
        w = 1000 / 18  # mol H2O per kg H2O

        for x, g in data.items():
            k = x * m  # mol K+
            n = m - k  # mol Na+

            # Set composition
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(w / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "K+"].set_value(k / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(n / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(m / (w + 2 * m))

            assert pytest.approx(0.06293706294, rel=1e-8) == value(
                model.state[1].mole_frac_phase_comp["Liq", "Cl-"]
            )
            assert pytest.approx(0.874125874, rel=1e-8) == value(
                model.state[1].mole_frac_phase_comp["Liq", "H2O"]
            )
            assert value(model.state[1].Liq_ionic_strength) == pytest.approx(
                0.06293706294, rel=1e-8
            )

            # Calculate mean activity coefficients for KCl
            ln_pm = value(
                0.5
                * (
                    model.state[1].Liq_log_gamma["K+"]
                    + model.state[1].Liq_log_gamma["Cl-"]
                )
            )

            # Convert to log base 10 and molal basis
            # Eqn 28 [1] - note that m=4 appears to be based on total salt
            log_m = ln_pm / log(10) - log((w + 2 * m) / w, 10)
            assert pytest.approx(g, rel=1e-3, abs=2e-3) == log_m


class TestUnsymmetric_Mixed_tau2(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        config = dict(configuration)
        eos_opt = config["phases"]["Liq"]["equation_of_state_options"] = {}
        eos_opt["reference_state"] = Unsymmetric

        m.params = GenericParameterBlock(**config)

        m.state = m.params.build_state_block([1])

        # Need to set a value of T for checking expressions later
        m.state[1].temperature.set_value(298.15)

        # Set salt-salt interaction parameter
        m.params.Liq.tau["K+, Cl-", "Na+, Cl-"].set_value(1)
        m.params.Liq.tau["Na+, Cl-", "K+, Cl-"].set_value(-1)

        return m

    @pytest.mark.unit
    def test_binary_parameters(self, model):
        # Make sure binary parameters are correct
        assert model.params.Liq.tau["H2O", "Na+, Cl-"].value == 9.0234
        assert model.params.Liq.tau["Na+, Cl-", "H2O"].value == -4.5916
        assert model.params.Liq.tau["H2O", "K+, Cl-"].value == 8.1354
        assert model.params.Liq.tau["K+, Cl-", "H2O"].value == -4.1341

    @pytest.mark.unit
    def test_log_gamma_NaCl(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 1 [2]
        # Form x_KCl: log(gamma_m NaCl)
        # Important notes, data appears to be in log_10 and on a molal basis
        # Thus, need to convert model output to this form
        data = {
            1e-12: -0.1056,
            0.0391817512877115: -0.116508693979689,
            0.0575717439293598: -0.121412096818314,
            0.0759617365710081: -0.126158740811656,
            0.0943517292126563: -0.131166649547137,
            0.112741721854305: -0.13603521708681,
            0.131131714495953: -0.140868949327531,
            0.149521707137601: -0.145693972743515,
            0.167911699779249: -0.150571249107926,
            0.186301692420898: -0.15535272840022,
            0.204691685062546: -0.160225650352262,
            0.223081677704194: -0.164889560510593,
            0.241471670345843: -0.170372846664985,
            0.252802877933121: -0.173974305960385,
            0.269892568064754: -0.177160750193646,
            0.288282560706402: -0.181995353316841,
            0.30667255334805: -0.18656259552058,
            0.325062545989698: -0.191349300107717,
            0.336765268579838: -0.193873013702975,
            0.388591611479029: -0.208149303584423,
            0.406981604120677: -0.21261726518615,
            0.425371596762325: -0.217625173921631,
            0.443761589403974: -0.222438004982981,
            0.462151582045622: -0.227101044258839,
            0.48054157468727: -0.232058441810839,
            0.498931567328918: -0.23687506743154,
            0.520665194996321: -0.242245171409142,
            0.539055187637969: -0.247122447773553,
            0.557445180279617: -0.25176458587004,
            0.575835172921266: -0.256532131042752,
            0.594225165562914: -0.26110111501144,
            0.612615158204562: -0.265913075190316,
            0.631005150846211: -0.27060746623523,
            0.649395143487859: -0.275210414620395,
            0.667785136129507: -0.279826426242668,
            0.686175128771155: -0.284546943761795,
            0.704565121412804: -0.28911070243564,
            0.722955114054452: -0.293752840532126,
            0.774781456953643: -0.30697536262028,
            0.793171449595291: -0.311634918366242,
            0.811561442236939: -0.31626834763799,
            0.829951434878587: -0.320771144538669,
            0.848341427520236: -0.325361029686727,
            0.866731420161884: -0.330055420731641,
            0.885121412803532: -0.334645305879699,
            0.90351140544518: -0.339209064553544,
            0.921901398086829: -0.343825076175816,
            0.940291390728477: -0.34851946722073,
            0.958681383370125: -0.353057099420361,
            0.974563649742458: -0.356916585242283,
        }

        m = 4  # molality of solution
        w = 1000 / 18  # mol H2O per kg H2O

        for x, g in data.items():
            k = x * m  # mol K+
            n = m - k  # mol Na+

            # Set composition
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(w / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "K+"].set_value(k / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(n / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(m / (w + 2 * m))

            assert pytest.approx(0.06293706294, rel=1e-8) == value(
                model.state[1].mole_frac_phase_comp["Liq", "Cl-"]
            )
            assert pytest.approx(0.874125874, rel=1e-8) == value(
                model.state[1].mole_frac_phase_comp["Liq", "H2O"]
            )
            assert value(model.state[1].Liq_ionic_strength) == pytest.approx(
                0.06293706294, rel=1e-8
            )

            # Calculate mean activity coefficient for NaCl
            # Natural log, mole fraction basis
            ln_pm = value(
                0.5
                * (
                    model.state[1].Liq_log_gamma["Na+"]
                    + model.state[1].Liq_log_gamma["Cl-"]
                )
            )

            # Convert to log base 10 and molal basis
            # Eqn 28 [1] - note that m=4 appears to be based on total salt
            log_m = ln_pm / log(10) - log((w + 2 * m) / w, 10)

            assert pytest.approx(g, rel=2e-2, abs=1e-3) == log_m

    @pytest.mark.unit
    def test_log_gamma_KCl(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data digitized from Fig 2 [2]
        # Form x_KCl: log(gamma_m KCl)
        # Important notes, data appears to be in log_10 and on a molal basis
        # Thus, need to convert model output to this form
        data = {
            0.0263950887050621: -0.0202162787518713,
            0.0458987565683367: -0.0247896391575565,
            0.0654024244316112: -0.0295050293895052,
            0.0849060922948858: -0.034021577864685,
            0.10440976015816: -0.038708562131381,
            0.123913428021435: -0.0432819225370662,
            0.143417095884709: -0.0478836889080041,
            0.162920763747984: -0.0525706731747001,
            0.182424431611259: -0.0571440335803853,
            0.201928099474533: -0.0617173939860705,
            0.221431767337808: -0.0664043782527665,
            0.240935435201082: -0.0714417027575792,
            0.260616409135841: -0.0755084901162578,
            0.279942770927631: -0.080124459469822,
            0.299446438790906: -0.0847830377712654,
            0.31895010665418: -0.0892711802811924,
            0.338453774517455: -0.0937593227911196,
            0.349092138806514: -0.0956625224630507,
            0.395191717392435: -0.106996502598879,
            0.41469538525571: -0.11131420931729,
            0.434199053118985: -0.115802351827217,
            0.453702720982259: -0.120148464510881,
            0.473206388845534: -0.124475639884376,
            0.492710056708808: -0.129234585929712,
            0.505786379480776: -0.132487897458466,
            0.522852088861142: -0.135680846311059,
            0.542355756724416: -0.140117858083531,
            0.561859424587691: -0.144369284216352,
            0.581363092450965: -0.148668053624595,
            0.60086676031424: -0.15324141403028,
            0.620370428177514: -0.157672744609701,
            0.639874096040789: -0.162018857293365,
            0.659377763904064: -0.166364969977029,
            0.678881431767338: -0.170597458799681,
            0.698385099630613: -0.175042992361729,
            0.717888767493887: -0.179460119958525,
            0.737392435357162: -0.18373454139655,
            0.756718797148952: -0.188085760390635,
            0.776399771083711: -0.192242804322241,
            0.788811196087613: -0.194818278505153,
            0.856187503251652: -0.209712472952653,
            0.875691171114927: -0.213888149844801,
            0.895194838978201: -0.218262668493717,
            0.914698506841476: -0.222438345385864,
            0.93420217470475: -0.226670834208517,
            0.951932781853182: -0.230439358932042,
        }

        m = 4  # molality of solution
        w = 1000 / 18  # mol H2O per kg H2O

        for x, g in data.items():
            k = x * m  # mol K+
            n = m - k  # mol Na+

            # Set composition
            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(w / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "K+"].set_value(k / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(n / (w + 2 * m))
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(m / (w + 2 * m))

            assert pytest.approx(0.06293706294, rel=1e-8) == value(
                model.state[1].mole_frac_phase_comp["Liq", "Cl-"]
            )
            assert pytest.approx(0.874125874, rel=1e-8) == value(
                model.state[1].mole_frac_phase_comp["Liq", "H2O"]
            )
            assert value(model.state[1].Liq_ionic_strength) == pytest.approx(
                0.06293706294, rel=1e-8
            )

            # Calculate mean activity coefficients for KCl
            ln_pm = value(
                0.5
                * (
                    model.state[1].Liq_log_gamma["K+"]
                    + model.state[1].Liq_log_gamma["Cl-"]
                )
            )

            # Convert to log base 10 and molal basis
            # Eqn 28 [1] - note that m=4 appears to be based on total salt
            log_m = ln_pm / log(10) - log((w + 2 * m) / w, 10)
            assert pytest.approx(g, rel=1e-3, abs=2.5e-3) == log_m
