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
import os
import json
from pyomo.common.fileutils import this_file_dir


def load_generic_ccs_costing_dictionary(path=None):

    """
    Custom dictionaries have been added as a way to add new scaling equations
    that are not based on the Bituminous Baseline report.
    New cost accounts include:

        - 5.1.a.epri can be used to cost the Cansolv carbon capture system
          between 90% and 97.7% capture rate. It is only valid for a flue gas CO2
          concentration of 4%. The additional cost component of the CO2 system is
          in account 5.1.b and remains unchanged.

        - 6.1 can be used to cost the MEA solvent-based carbon capture system,
          using the reference: Kangkang Li, Ashleigh Cousins, Hai You, Paul Feron,
          Weilang Luo, Jian Chen (2016). Systematic study of aqueous
          monoethanolamine (MEA)-based CO2 capture process: Techno-economic
          assessment of the MEA process and its improvements. Applied Energy, 165,
          648-659.
    """
    if path is None:
        directory = this_file_dir()
    else:
        directory = path

    generic_ccs_costing_exponents = {
        "6": {
            "5.1.a.epri": {
                "Account Name": "Cansolv CO2 Removal System",
                "Exponent": 2.788,
                "Process Parameter": "CO2 Flowrate",
            },
            "6.1.ccs": {
                "Account Name": "MEA solvent capture system absorber",
                "Exponent": 0.6,
                "Process Parameter": "Absorber volume",
            },
            "6.2.ccs": {
                "Account Name": "MEA solvent capture system absorber packing",
                "Exponent": 0.6,
                "Process Parameter": "Absorber packing volume",
            },
            "6.3.ccs": {
                "Account Name": "MEA solvent capture system stripper",
                "Exponent": 0.6,
                "Process Parameter": "Stripper volume",
            },
            "6.4.ccs": {
                "Account Name": "MEA solvent capture system stripper packing",
                "Exponent": 0.6,
                "Process Parameter": "Stripper packing volume",
            },
            "6.5.ccs": {
                "Account Name": "MEA solvent capture system stripper condenser",
                "Exponent": 0.6,
                "Process Parameter": "Stripper condenser area",
            },
            "6.6.ccs": {
                "Account Name": "MEA solvent capture system stripper reboiler",
                "Exponent": 0.6,
                "Process Parameter": "Stripper reboiler area",
            },
            "6.7.ccs": {
                "Account Name": "MEA solvent capture system lean rich heat exchanger",
                "Exponent": 0.6,
                "Process Parameter": "Lean rich heat exchanger area",
            },
            "6.8.ccs": {
                "Account Name": "MEA solvent capture system lean solvent cooler",
                "Exponent": 0.6,
                "Process Parameter": "Lean solvent cooler area",
            },
            "6.9.1.ccs": {
                "Account Name": "MEA solvent capture system flue gas blower cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.9.2.ccs": {
                "Account Name": "MEA solvent capture system flue gas blower cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.10.1.ccs": {
                "Account Name": "MEA solvent capture system flue gas direct contact cooler cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.10.2.ccs": {
                "Account Name": "MEA solvent capture system flue gas direct contact cooler cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.11.1.ccs": {
                "Account Name": "MEA solvent capture system flue gas direct contact cooler packing cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.11.2.ccs": {
                "Account Name": "MEA solvent capture system flue gas direct contact cooler packing cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.12.1.ccs": {
                "Account Name": "MEA solvent capture system pretreatment pump cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.12.2.ccs": {
                "Account Name": "MEA solvent capture system pretreatment pump cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.13.1.ccs": {
                "Account Name": "MEA solvent capture system pretreatment cooler cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.13.2.ccs": {
                "Account Name": "MEA solvent capture system pretreatment cooler cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.14.1.ccs": {
                "Account Name": "MEA solvent capture system pretreatment tank cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.14.2.ccs": {
                "Account Name": "MEA solvent capture system pretreatment tank cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.15.1.ccs": {
                "Account Name": "MEA solvent capture system washing column cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.15.2.ccs": {
                "Account Name": "MEA solvent capture system washing column cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.16.1.ccs": {
                "Account Name": "MEA solvent capture system washing column packing cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.16.2.ccs": {
                "Account Name": "MEA solvent capture system washing column packing cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.17.1.ccs": {
                "Account Name": "MEA solvent capture system washing solvent cooler cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.17.2.ccs": {
                "Account Name": "MEA solvent capture system washing solvent cooler cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.18.1.ccs": {
                "Account Name": "MEA solvent capture system washing solvent pump cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.18.2.ccs": {
                "Account Name": "MEA solvent capture system washing solvent pump cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.19.1.ccs": {
                "Account Name": "MEA solvent capture system condenser pump cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.19.2.ccs": {
                "Account Name": "MEA solvent capture system condenser pump cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.20.1.ccs": {
                "Account Name": "MEA solvent capture system stripper reflux drum cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.20.2.ccs": {
                "Account Name": "MEA solvent capture system stripper reflux drum cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.21.1.ccs": {
                "Account Name": "MEA solvent capture system lean solvent pump cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.21.2.ccs": {
                "Account Name": "MEA solvent capture system lean solvent pump cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.22.1.ccs": {
                "Account Name": "MEA solvent capture system solvent storage tank cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.22.2.ccs": {
                "Account Name": "MEA solvent capture system solvent storage tank cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.23.1.ccs": {
                "Account Name": "MEA solvent capture system washing solvent tank cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.23.2.ccs": {
                "Account Name": "MEA solvent capture system washing solvent tank cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.24.1.ccs": {
                "Account Name": "MEA solvent capture system solvent stripper reclaimer cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.24.2.ccs": {
                "Account Name": "MEA solvent capture system solvent stripper reclaimer cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.25.1.ccs": {
                "Account Name": "MEA solvent capture system solvent reclaimer cooler cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.25.2.ccs": {
                "Account Name": "MEA solvent capture system solvent reclaimer cooler cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
            "6.26.1.ccs": {
                "Account Name": "MEA solvent capture system solvent filtration cost component 1",
                "Exponent": 0.6,
                "Process Parameter": "CO2 product flow rate",
            },
            "6.26.2.ccs": {
                "Account Name": "MEA solvent capture system solvent filtration cost component 2",
                "Exponent": 0.6,
                "Process Parameter": "Flue gas inlet to absorber",
            },
        }
    }

    generic_ccs_costing_params = {
        "6": {
            "B": {
                "5.1.a.epri": {
                    "BEC": 224191.4,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.18,
                    "Project Contingency": 0.2,
                    "RP Value": 493587.88,
                    "Units": "lb/hr",
                },
                "6.1.ccs": {
                    "BEC": 8376000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 1074.424688,
                    "Units": "m**3",
                },
                "6.2.ccs": {
                    "BEC": 6117000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 791.6813487,
                    "Units": "m**3",
                },
                "6.3.ccs": {
                    "BEC": 2273000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 419.6971436,
                    "Units": "m**3",
                },
                "6.4.ccs": {
                    "BEC": 1464000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 309.2505268,
                    "Units": "m**3",
                },
                "6.5.ccs": {
                    "BEC": 305000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 800,
                    "Units": "m**2",
                },
                "6.6.ccs": {
                    "BEC": 3183000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 4250,
                    "Units": "m**2",
                },
                "6.7.ccs": {
                    "BEC": 1502000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 9100,
                    "Units": "m**2",
                },
                "6.8.ccs": {
                    "BEC": 522000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 1200,
                    "Units": "m**2",
                },
                "6.9.1.ccs": {
                    "BEC": 0.6 * 812000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.9.2.ccs": {
                    "BEC": 0.4 * 812000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.10.1.ccs": {
                    "BEC": 0.6 * 2875000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.10.2.ccs": {
                    "BEC": 0.4 * 2875000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.11.1.ccs": {
                    "BEC": 0.6 * 1883000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.11.2.ccs": {
                    "BEC": 0.4 * 1883000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.12.1.ccs": {
                    "BEC": 0.6 * 108000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.12.2.ccs": {
                    "BEC": 0.4 * 108000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.13.1.ccs": {
                    "BEC": 0.6 * 202000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.13.2.ccs": {
                    "BEC": 0.4 * 202000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.14.1.ccs": {
                    "BEC": 0.6 * 93000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.14.2.ccs": {
                    "BEC": 0.4 * 93000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.15.1.ccs": {
                    "BEC": 0.6 * 2842000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.15.2.ccs": {
                    "BEC": 0.4 * 2842000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.16.1.ccs": {
                    "BEC": 0.6 * 2060000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.16.2.ccs": {
                    "BEC": 0.4 * 2060000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.17.1.ccs": {
                    "BEC": 0.6 * 69000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.17.2.ccs": {
                    "BEC": 0.4 * 69000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.18.1.ccs": {
                    "BEC": 0.6 * 17000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.18.2.ccs": {
                    "BEC": 0.4 * 17000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.19.1.ccs": {
                    "BEC": 0.6 * 35000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.19.2.ccs": {
                    "BEC": 0.4 * 35000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.20.1.ccs": {
                    "BEC": 0.6 * 46000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.20.2.ccs": {
                    "BEC": 0.4 * 46000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.21.1.ccs": {
                    "BEC": 0.6 * 159000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.21.2.ccs": {
                    "BEC": 0.4 * 159000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.22.1.ccs": {
                    "BEC": 0.6 * 342000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.22.2.ccs": {
                    "BEC": 0.4 * 342000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.23.1.ccs": {
                    "BEC": 0.6 * 46000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.23.2.ccs": {
                    "BEC": 0.4 * 46000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.24.1.ccs": {
                    "BEC": 0.6 * 181000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.24.2.ccs": {
                    "BEC": 0.4 * 181000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.25.1.ccs": {
                    "BEC": 0.6 * 180000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.25.2.ccs": {
                    "BEC": 0.4 * 180000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
                "6.26.1.ccs": {
                    "BEC": 0.6 * 866000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 262349.78,
                    "Units": "lb/hr",
                },
                "6.26.2.ccs": {
                    "BEC": 0.4 * 866000 / 1e3,
                    "Eng Fee": 0.2,
                    "Process Contingency": 0.12,
                    "Project Contingency": 0.2,
                    "RP Value": 672504.35,
                    "Units": "m**3/hr",
                },
            }
        }
    }

    # this is a way to "zip together" two generic_ccs_costing dictionaries in the formats above
    # the BEC_units happen to be in thousands of 2013 USD for these accounts
    # for the power plant methods to work, BEC_units must be "$year", "K$year" or
    # "$Myear" for USD, thousands USD and millions USD, respectively

    if not os.path.exists(
        os.path.join(directory, "generic_ccs_costing_data.json")
    ):  # make the dictionary

        generic_ccs_costing_data = generic_ccs_costing_params

        for tech in generic_ccs_costing_data.keys():  # do one technology at a time
            for ccs in ["A", "B"]:
                if (
                    ccs in generic_ccs_costing_data[tech]
                ):  # check if CCS = A, for indexing
                    accounts_dict = generic_ccs_costing_data[tech][ccs]  # shorter alias
                    for account in accounts_dict.keys():  # do one account at a time
                        if account == "5.1.a.epri":
                            accounts_dict[account][
                                "BEC_units"
                            ] = "K$2018"  # add BEC units as thousands of 2018 USD
                        else:
                            accounts_dict[account][
                                "BEC_units"
                            ] = "K$2013"  # add BEC units as thousands of 2013 USD
                        for accountkey in generic_ccs_costing_exponents[tech][
                            account
                        ].keys():  # get one " exponents"account property at a time
                            accounts_dict[account][
                                accountkey
                            ] = generic_ccs_costing_exponents[tech][account][accountkey]
                        sorted_accountkeys = sorted(
                            accounts_dict[account]
                        )  # now, sort the accountkeys alphabetically within each account
                        accounts_dict[account] = {
                            key: accounts_dict[account][key]
                            for key in sorted_accountkeys
                        }  # re-add the keys in alphabetical order
                    generic_ccs_costing_data[tech][
                        ccs
                    ] = accounts_dict  # use the alias to update the original dictionary

        with open(
            os.path.join(directory, "generic_ccs_costing_data.json"), "w"
        ) as outfile:
            json.dump(generic_ccs_costing_data, outfile)
        print("Success! New costing dictionary generated.")

    # assuming the dictionary exists, load it so it is importable when called
    with open(os.path.join(directory, "generic_ccs_costing_data.json"), "r") as file:
        generic_ccs_costing_params = json.load(file)
    return generic_ccs_costing_params
