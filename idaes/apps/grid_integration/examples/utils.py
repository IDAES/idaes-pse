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

try:
    # Pyton 3.8+
    from importlib import resources
except ImportError:
    # Python 3.7
    import importlib_resources as resources
from pathlib import Path
import os

import pandas as pd


with resources.path("idaes.tests.prescient.5bus", "__init__.py") as pkg_file:
    prescient_5bus = Path(pkg_file).parent

# rts_gmlc_generator_dataframe = pd.read_csv("gen.csv")
rts_gmlc_generator_dataframe = pd.read_csv(os.path.join(prescient_5bus, "gen.csv"))
rts_gmlc_bus_dataframe = pd.read_csv(os.path.join(prescient_5bus, "bus.csv"))

# price forecasts utils
daily_da_price_means = [
    18.12183987,
    16.95203894,
    21.35542779,
    16.28283605,
    18.05356279,
    26.75966193,
    15.02385805,
    7.97800459,
    6.75568025,
    7.35747065,
    8.10476168,
    9.42647617,
    10.78553799,
    12.45377683,
    14.79835537,
    17.10387646,
    25.31337418,
    29.34686351,
    41.59693577,
    26.85816007,
    28.68259094,
    23.51281011,
    20.69422377,
    18.55041211,
]
daily_rt_price_means = [
    646.0645177,
    688.56564529,
    767.9131885,
    734.9646116,
    642.81583885,
    1044.96253081,
    1281.25714066,
    1035.88602594,
    137.74265062,
    102.9523983,
    70.00893574,
    70.0554967,
    41.63036343,
    51.43459105,
    52.2201115,
    23.76061729,
    36.56761717,
    201.50946959,
    523.27255929,
    529.40450156,
    555.14680986,
    645.38156746,
    777.07099054,
    845.98180641,
]
daily_da_price_stds = [
    8.52102171,
    9.19761893,
    75.304164,
    9.30600699,
    8.07077925,
    92.92510996,
    13.94045029,
    9.83439913,
    9.47028416,
    9.83386367,
    10.28061207,
    10.84413467,
    11.20525267,
    11.42147398,
    11.16693943,
    10.5249286,
    53.2882641,
    52.90611362,
    117.15560661,
    6.21862006,
    54.40205825,
    5.33925754,
    7.16726919,
    8.47854486,
]
daily_rt_price_stds = [
    2305.18158192,
    2378.86464441,
    2542.18800864,
    2526.03313382,
    2329.00046949,
    2956.7654427,
    3185.27587193,
    2876.75731403,
    1049.04095182,
    909.72152018,
    743.16499996,
    742.97617584,
    355.8614282,
    480.37033221,
    509.39278934,
    110.16004212,
    400.26846307,
    1233.7527661,
    2009.23951093,
    2037.43290886,
    2067.75905496,
    2238.26562172,
    2558.9623294,
    2641.99925492,
]
