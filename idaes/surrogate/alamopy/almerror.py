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
Errors
"""
from idaes.surrogate.alamopy import writethis


def almerror(code):
    # reports errors encountered

    if (code == 'p1'):
        writethis('Inconsistent number of data points provided between input arrays')
        quit()
    elif (code == 'p2'):
        writethis('Error in structure of validation data')
        quit()
    elif (code == 'p3'):
        writethis('Inputs option not understood')
        quit()


class AlamoError(Exception):
    pass


class AlamoInputError(AlamoError):
    def __init__(self, msg):
        msg = 'Bad input to doalamo(): {}'.format(msg)
        super(AlamoInputError, self).__init__(msg)
