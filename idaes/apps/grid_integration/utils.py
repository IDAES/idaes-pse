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


def convert_marginal_costs_to_actual_costs(power_marginal_cost_pairs):

    """
    Convert a list of power and marginal cost pairs to a list of power and actual
    cost pairs.

    Args:
        power_marginal_cost_pairs: a list of power and marginal cost pairs (tuple)

    Returns:
        list: a list of power and actual cost pairs (tuple).
    """

    # sort by power
    power_marginal_cost_pairs.sort(key=lambda pair: pair[0])

    # initialize
    actual_costs = []
    pre_power = 0
    pre_cost = 0

    for power, marginal_cost in power_marginal_cost_pairs:

        delta_p = power - pre_power
        cur_cost = pre_cost + marginal_cost * delta_p
        actual_costs.append((power, cur_cost))

        pre_power = power
        pre_cost += marginal_cost * delta_p

    return actual_costs
