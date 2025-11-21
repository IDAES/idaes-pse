#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
import os
import re
import sys

__author__ = "Dan Gunter"


def git_hash():
    """Get current git hash, with no dependencies on external packages."""
    # find git root (in grandparent dir to this file, if anywhere)
    git_root = os.path.realpath(os.path.join(__file__, "..", "..", ".git"))
    if not os.path.exists(git_root) or not os.path.isdir(git_root):
        raise ValueError(f"git root '{git_root}' not found")
    # get HEAD ref's file
    try:
        head = open(os.path.join(git_root, "HEAD"))
    except FileNotFoundError as err:
        raise ValueError(f"cannot open HEAD: {err}")
    # parse file looking for 'ref: <path>'
    head_ref = None
    for line in head:
        ref_match = re.match(r"ref:\s+(\S+)", line)
        if ref_match:
            head_ref = ref_match.group(1)
            break
    head.close()
    if head_ref is None:
        raise ValueError(f"no ref found in HEAD '{head}'")
    # read value of ref in <path> found previously
    ref_file_path = os.path.join(git_root, head_ref)
    try:
        ref_file = open(ref_file_path)
    except FileNotFoundError:
        raise ValueError(f"ref file '{ref_file_path}' not found")
    ref = ref_file.read().strip()
    ref_file.close()
    return ref


# Get git hash. No output unless IDAES_DEBUG is set in env
gh = None
try:
    try:
        gh = git_hash()
        if os.environ.get("IDAES_DEBUG", None):
            print(f"git hash = {gh}", file=sys.stderr)
    except ValueError as err:
        if os.environ.get("IDAES_DEBUG", None):
            print(f"git_hash() error: {err}", file=sys.stderr)
except NameError:  # eg, if invoked from setup.py
    pass
