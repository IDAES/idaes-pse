RIPE: Reaction Identification and Parameter Estimation
=======================================================

.. warning::
  RIPE was removed from the IDAES codebase in v2.1 due to a lack of use and maintenance (3/9/2023, commit hash: 2c563ad). RIPE can be found in older releases if required, and maybe restored if suffiicent interest is shown and a maintainer can be identified (see https://github.com/IDAES/idaes-pse/wiki/Maintenance-and-Testing-of-User-Supported-Code)


The RIPE module provides tools for reaction network identification. RIPE uses reactor data consisting of concentration, or conversion, values for multiple species that are obtained dynamically, or at multiple process conditions (temperatures, flow rates, working volumes) to identify probable reaction kinetics. The RIPE module also contains tools to facilitate adaptive experimental design. The experimental design tools in RIPE require the use of the python package RBFopt. More information for RBFopt is available at www.github.com/coin-or/rbfopt

