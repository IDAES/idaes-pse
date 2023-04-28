Block Triangularization Initializer
===================================

The Block Triangularization initializer is a general purpose ``Initializer`` that can be applied to any model type. It applies the Pyomo incidence analysis toolbox to attempt to decompose the model into the maximum set of smaller sub-problems that can be solved sequentially using a block traingularization transformation. These sub-problems can then be solved by applying either a simple Newton solver (for 1x1 blocks, i.e., one variable and constraint) or a user selected solver for NxN blocks. The Block Triangularization initializer also inherently handles plug-ins, as these are decomposed as part of the overall model structure.

This approach is often faster and more reliable than heuristic methods, however it may struggle to decompose problems involving tightly coupled systems of equations such as column models, counter-current flow and vapor-liquid equilibrium.

The `BlockTraingulariaztionInitializer`` is the default ``Initializer`` assigned to all IDAES property packages (via ``model.default_initializer`` unless this is overwritten by the model developer.

.. module:: idaes.core.initialization.block_triangularization

BlockTriangularizationInitializer Class
---------------------------------------

.. autoclass:: BlockTriangularizationInitializer
  :members: precheck, initialize
