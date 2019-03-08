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
Test plots (at least creation, not display) in idaes.vis.bokeh_plots
"""
# third-party
from pandas import DataFrame
# package
from idaes.vis import bokeh_plots
from idaes.vis.plotutils import HENStreamType


exchangers = [
    {'hot': 'H2', 'cold': 'C1', 'Q': 1400, 'A': 159, 'annual_cost': 28358, 'stg': 2},
    {'hot': 'H1', 'cold': 'C1', 'Q': 667, 'A': 50, 'annual_cost': 10979, 'stg': 3},
    {'hot': 'H1', 'cold': 'C1', 'Q': 233, 'A': 10, 'annual_cost': 4180, 'stg': 1},
    {'hot': 'H1', 'cold': 'C2', 'Q': 2400, 'A': 355, 'annual_cost': 35727, 'stg': 2},
    {
        'hot': 'H2',
        'cold': 'W',
        'Q': 400,
        'A': 50,
        'annual_cost': 10979,
        'stg': 3,
        'utility_type': HENStreamType.cold_utility,
    },
    {
        'hot': 'S',
        'cold': 'C2',
        'Q': 450,
        'A': 50,
        'annual_cost': 0,
        'stg': 1,
        'utility_type': HENStreamType.hot_utility,
    },
]

streams = [
    {'name': 'H2', 'temps': [423, 423, 330, 303], 'type': HENStreamType.hot},
    {'name': 'H1', 'temps': [443, 435, 355, 333], 'type': HENStreamType.hot},
    {'name': 'C1', 'temps': [408, 396, 326, 293], 'type': HENStreamType.cold},
    {'name': 'C2', 'temps': [413, 413, 353, 353], 'type': HENStreamType.cold},
]


def test_heat_exchanger():
    bokeh_plots.HeatExchangerNetwork(
        exchangers, streams, mark_temperatures_with_tooltips=True
    )

# profile_data = DataFrame(columns=[range(10),)
# def test_profile():
#     bokeh_plots.ProfilePlot()