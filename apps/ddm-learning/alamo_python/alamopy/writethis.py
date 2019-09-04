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


def writethis(*arg):
    # This function is used to write to the terminal
    import sys
    for pnt in arg:
        if isinstance(pnt, str):
            sys.stdout.write(pnt)
        else:
            sys.stdout.write('!!! Printing Variable !!!')
            sys.stdout.write(pnt)
