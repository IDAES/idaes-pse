# -*- coding: UTF-8 -*-
##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
This module contains common constants of use in process systems engineering.

All units are SI, and generally expressed to 9 signficant figures where
avaialble
"""

__author__ = "Andrew Lee"

# -----------------------------------------------------------------------------
# General geometric relationships
import math
pi = math.pi

# -----------------------------------------------------------------------------
# Constants used as fundamental definitions in SI, from
# https://www.bipm.org/utils/common/pdf/si-brochure/SI-Brochure-9.pdf

avogadro_number = 6.02214076e23  # unitless or mol^-1

boltzmann_constant = 1.38064900e-23  # J⋅K^-1

elemental_charge = 1.602176634e-19  # C

planck_constant = 6.62607015e-34  # J⋅s

speed_light = 299792458  # in a vacuum m⋅s^-1

# -----------------------------------------------------------------------------
# Constants derived from fundamental constants

# Faraday constant = elemental charge * Avogadro's constant
faraday_constant = 96485.33212  # C⋅mol^-1

# Gas constant = Avogadro's constant * Boltzmann's constant
gas_constant = 8.314462618  # J⋅mol^-1⋅K^-1

# Stefan-Boltzmann constant
# Function of Boltzmann constant, pi and speed of light
stefan_constant = 5.67037442e-8  # W⋅m^−2⋅K^−4

# -----------------------------------------------------------------------------
# Other constants - all values sourced from NIST to avaialble uncertainty
# All values retrieved 8th Jan 2020

# https://physics.nist.gov/cgi-bin/cuu/Value?gn
acceleration_gravity = 9.80665  # m⋅s^-2

# https://physics.nist.gov/cgi-bin/cuu/Value?bg
gravitational_constant = 6.67430e-11  # m^3⋅kg^−1⋅s^−2

# https://physics.nist.gov/cgi-bin/cuu/Value?me
mass_electron = 9.1093837015e-31  # kg
