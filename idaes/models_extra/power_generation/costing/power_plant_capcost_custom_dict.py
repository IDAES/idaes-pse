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

custom_costing_exponents = {
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
        "6.9.ccs": {
            "Account Name": "MEA solvent capture system flue gas blower",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.10.ccs": {
            "Account Name": "MEA solvent capture system flue gas direct contact cooler",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.11.ccs": {
            "Account Name": "MEA solvent capture system flue gas direct contact cooler packing",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.12.ccs": {
            "Account Name": "MEA solvent capture system pretreatment pump",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.13.ccs": {
            "Account Name": "MEA solvent capture system pretreatment cooler",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.14.ccs": {
            "Account Name": "MEA solvent capture system pretreatment tank",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.15.ccs": {
            "Account Name": "MEA solvent capture system washing column",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.16.ccs": {
            "Account Name": "MEA solvent capture system washing column packing",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.17.ccs": {
            "Account Name": "MEA solvent capture system washing solvent cooler",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.18.ccs": {
            "Account Name": "MEA solvent capture system washing solvent pump",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.19.ccs": {
            "Account Name": "MEA solvent capture system condenser pump",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.20.ccs": {
            "Account Name": "MEA solvent capture system stripper reflux drum",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.21.ccs": {
            "Account Name": "MEA solvent capture system lean solvent pump",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.22.ccs": {
            "Account Name": "MEA solvent capture system solvent storage tank",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.23.ccs": {
            "Account Name": "MEA solvent capture system washing solvent tank",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.24.ccs": {
            "Account Name": "MEA solvent capture system solvent stripper reclaimer",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.25.ccs": {
            "Account Name": "MEA solvent capture system solvent reclaimer cooler",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
        "6.26.ccs": {
            "Account Name": "MEA solvent capture system solvent filtration",
            "Exponent": 0.6,
            "Process Parameter": [
                "CO2 product flow rate",
                "Flue gas inlet to absorber",
            ],
        },
    }
}

custom_costing_params = {
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
                "BEC": 6128000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": 1074.424688,
                "Units": "m**3",
            },
            "6.2.ccs": {
                "BEC": 6040000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": 791.6813487,
                "Units": "m**3",
            },
            "6.3.ccs": {
                "BEC": 1986000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": 419.6971436,
                "Units": "m**3",
            },
            "6.4.ccs": {
                "BEC": 1438000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": 309.2505268,
                "Units": "m**3",
            },
            "6.5.ccs": {
                "BEC": 260000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": 800,
                "Units": "m**2",
            },
            "6.6.ccs": {
                "BEC": 3095000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": 4250,
                "Units": "m**2",
            },
            "6.7.ccs": {
                "BEC": 1151000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": 9100,
                "Units": "m**2",
            },
            "6.8.ccs": {
                "BEC": 465000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": 1200,
                "Units": "m**2",
            },
            "6.9.ccs": {
                "BEC": 731000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.10.ccs": {
                "BEC": 2082000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.11.ccs": {
                "BEC": 1855000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.12.ccs": {
                "BEC": 89000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.13.ccs": {
                "BEC": 165000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.14.ccs": {
                "BEC": 74000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.15.ccs": {
                "BEC": 1992000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.16.ccs": {
                "BEC": 2036000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.17.ccs": {
                "BEC": 46000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.18.ccs": {
                "BEC": 9000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.19.ccs": {
                "BEC": 25000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.20.ccs": {
                "BEC": 34000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.21.ccs": {
                "BEC": 260000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.22.ccs": {
                "BEC": 296000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.23.ccs": {
                "BEC": 34000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.24.ccs": {
                "BEC": 144000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.25.ccs": {
                "BEC": 135000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
            "6.26.ccs": {
                "BEC": 791000 / 1e3,
                "Eng Fee": 0.2,
                "Process Contingency": 0.12,
                "Project Contingency": 0.2,
                "RP Value": [262349.78, 672504.35],
                "Cost scaling fraction": [0.6, 0.4],
                "Units": ["lb/hr", "m3/hr"],
            },
        }
    }
}
