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
This module contains utility functions for initialization of IDAES models.
"""

from pyomo.environ import Block

__author__ = "Andrew Lee, John Siirola"


# HACK, courtesy of J. Siirola
def solve_indexed_blocks(solver, blocks, **kwds):
    """
    This method allows for solving of Indexed Block components as if they were
    a single Block. A temporary Block object is created which is populated with
    the contents of the objects in the blocks argument and then solved.

    Args:
        solve : a Pyomo solver object to use when solving the Indexed Block
        blocks : an object which inherits from Block, or a list of Blocks
        kwds : a dict of argumnets to be passed to the solver

    Returns:
        A Pyomo solver results object
    """
    # Check blocks argument, and convert to a list of Blocks
    if isinstance(blocks, Block):
        blocks = [blocks]

    try:
        # Create a temporary Block
        tmp = Block(concrete=True)

        nBlocks = len(blocks)

        # Iterate over indexed objects
        for i, b in enumerate(blocks):
            # Check that object is a Block
            if not isinstance(b, Block):
                raise TypeError("Trying to apply solve_indexed_blocks to "
                                "object containing non-Block objects")
            # Append components of BlockData to temporary Block
            try:
                tmp._decl["block_%s" % i] = i
                tmp._decl_order.append((b, i+1 if i < nBlocks-1 else None))
            except Exception:
                raise Exception("solve_indexed_blocks method failed adding "
                                "components to temporary block.")

        # Set ctypes on temporary Block
        tmp._ctypes[Block] = [0, nBlocks-1, nBlocks]

        # Solve temporary Block
        results = solver.solve(tmp, **kwds)

    finally:
        # Clean up temporary Block contents so they are not removed when Block
        # is garbage collected.
        tmp._decl = {}
        tmp._decl_order = []
        tmp._ctypes = {}

    # Return results
    return results
