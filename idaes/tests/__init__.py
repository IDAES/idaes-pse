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
import logging
import os


def level_num(name, default_level):
    """Get level num from name"""
    level_mapping = dict(
        ERROR=logging.ERROR,
        INFO=logging.INFO,
        WARNING=logging.WARNING,
        DEBUG=logging.DEBUG,
    )
    return level_mapping.get(name.upper(), default_level)


# set up test logger
_log = logging.getLogger(__name__)

level = logging.WARNING
if "IDAES_TEST_LOG_LEVEL" in os.environ:
    env_level = os.environ["IDAES_TEST_LOG_LEVEL"]
    level = level_num(env_level, level)

_h = logging.StreamHandler()
_h.setFormatter(
    logging.Formatter(
        "%(asctime)s [%(levelname)s] " "%(filename)s:%(lineno)d :: %(message)s"
    )
)
_log.addHandler(_h)
_log.setLevel(level)
_log.propagate = False
