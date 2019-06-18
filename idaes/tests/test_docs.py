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
Test doc build.

This does *not* try to build the docs. Instead it looks
for the log of the build and looks for errors in it.
"""
# stdlib
import logging
import os
import re

# third-party
import pytest


_log = logging.getLogger(__name__)


@pytest.fixture
def docs_path():
    """Find the docs.

    We have to assume that the docs are somewhere
    relative to either this file or the user's current
    working directory.
    Start with the working directory.

    Returns:
        str: Doc path or None if not found.
    """
    cwd = os.path.realpath(os.getcwd())
    # look for "./docs"
    path = os.path.join(cwd, "docs")
    if os.path.exists(path):
        return path
    # look for (last) "docs/" in full current working dir
    cwd_comp = cwd.split(os.path.sep)
    n = len(cwd_comp)
    try:
        idx = list(reversed(cwd_comp)).index("docs")
        return os.path.sep.join(cwd_comp[: n - idx])
    except ValueError:
        pass
    # look for docs relative to this file
    try:
        path = os.path.join(os.path.dirname(__file__), "..", "..", "docs")
        # look for {this file}/../../docs
        if os.path.exists(path):
            return path
        # look for "docs/" in file's full path
        file_comp = path.split(os.path.sep)
        n = len(file_comp)
        try:
            idx = list(reversed(file_comp)).index("docs")
            return os.path.sep.join(file_comp[: n - idx])
        except ValueError:
            pass
    except NameError:
        pass  # __file__ not defined(?)


ERRLOG = "sphinx-errors.txt"


def test_sphinx_build_log(docs_path):
    _log.info('docs path = "{}"'.format(docs_path))
    if docs_path is None:
        _log.warning('Could not find "docs" directory')
        return
    log_path = os.path.join(docs_path, ERRLOG)
    if not os.path.exists(log_path):
        _log.warning(
            'Could not find "{}" in docs directory: {}'.format(ERRLOG, log_path)
        )
        return
    log = open(log_path)
    for line in log:
        if "WARNING: " in line:
            pass
            # assert re.search(r'duplicate label sub(module|package)s', line), \
            #    'Non-trivial warning: {}'.format(line.strip())
        elif "ERROR: " in line:
            assert False, "Error: {}".format(line.strip())
