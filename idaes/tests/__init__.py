###############################################################################
# Copyright
# =========
#
# Institute for the Design of Advanced Energy Systems Process Systems Engineering
# Framework (IDAES PSE Framework) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021 by the
# software owners: The Regents of the University of California, through Lawrence
# Berkeley National Laboratory,  National Technology & Engineering Solutions of
# Sandia, LLC, Carnegie Mellon University, West Virginia University Research
# Corporation, et al.  All rights reserved.
#
# NOTICE.  This Software was developed under funding from the U.S. Department of
# Energy and the U.S. Government consequently retains certain rights. As such, the
# U.S. Government has been granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to
# reproduce, distribute copies to the public, prepare derivative works, and
# perform publicly and display publicly, and to permit other to do so. Copyright
# (C) 2018-2019 IDAES - All Rights Reserved
#
###############################################################################
import logging
import os


def level_num(name, default_level):
    """Get level num from name"""
    level_mapping = dict(ERROR=logging.ERROR, INFO=logging.INFO,
                         WARNING=logging.WARNING, DEBUG=logging.DEBUG)
    return level_mapping.get(name.upper(), default_level)


# set up test logger
_log = logging.getLogger(__name__)

level = logging.WARNING
if 'IDAES_TEST_LOG_LEVEL' in os.environ:
    env_level = os.environ['IDAES_TEST_LOG_LEVEL']
    level = level_num(env_level, level)

_h = logging.StreamHandler()
_h.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] '
                                  '%(filename)s:%(lineno)d :: %(message)s'))
_log.addHandler(_h)
_log.setLevel(level)
_log.propagate = False