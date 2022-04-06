# -*- coding: UTF-8 -*-
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
This module contains common constants of use in process systems engineering.

All units are SI, and generally expressed to 9 signficant figures where
avaialble
"""

__author__ = "Andrew Lee"

from pyomo.environ import units

import math


class Constants:
    # -------------------------------------------------------------------------
    # General geometric relationships
    pi = math.pi

    # -------------------------------------------------------------------------
    # Constants used as fundamental definitions in SI, from
    # https://www.bipm.org/utils/common/pdf/si-brochure/SI-Brochure-9.pdf
    avogadro_number = 6.02214076e23 / units.mol

    boltzmann_constant = 1.38064900e-23 * units.joule / units.degK

    elemental_charge = 1.602176634e-19 * units.coulomb

    planck_constant = 6.62607015e-34 * units.joule * units.second

    speed_light = 299792458 * units.m / units.s  # in a vacuum

    # -------------------------------------------------------------------------
    # Constants derived from fundamental constants

    # Faraday constant = elemental charge * Avogadro's constant
    faraday_constant = 96485.33212 * units.coulomb / units.mol

    # Gas constant = Avogadro's constant * Boltzmann's constant
    gas_constant = 8.314462618 * units.joule / units.mol / units.degK

    # Stefan-Boltzmann constant
    # Function of Boltzmann constant, pi and speed of light
    stefan_constant = 5.67037442e-8 * units.watt / units.metre**2 / units.degK**4

    # -------------------------------------------------------------------------
    # Other constants - all values sourced from NIST to avaialble uncertainty
    # All values retrieved 8th Jan 2020 unless otherwise noted

    # https://physics.nist.gov/cgi-bin/cuu/Value?gn
    acceleration_gravity = 9.80665 * units.metre / units.second**2

    # https://physics.nist.gov/cgi-bin/cuu/Value?bg
    gravitational_constant = (
        6.67430e-11 * units.metre**3 / units.kg / units.second**2
    )

    # https://physics.nist.gov/cgi-bin/cuu/Value?me
    mass_electron = 9.1093837015e-31 * units.kilogram

    # https://physics.nist.gov/cgi-bin/cuu/Value?ep0
    # 8th April 2021
    vacuum_electric_permittivity = 8.8541878128e-12 * units.farad * units.metre**-1
