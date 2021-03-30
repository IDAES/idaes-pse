Equilibrium Reaction Forms
==========================

.. contents:: Contents 
    :depth: 2

Power Law (power_law_equil)
---------------------------

The method uses a power law form using the concentration form provided  to calculate the reaction rate.

.. math:: k_{eq} = \prod_{(p, j)}{x_{(p,j)}^{O_{(p,j)}}}

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Indices", "Description"

   ":math:`O`", "reaction_order", "phase, component", "Reaction order"

Providing a `reaction_order` dict is optional. If one is not provided, it will be assumed that this is an elementary reaction and that the reaction order is equal to the stoichiometric coefficient for all component in non-solid phases (the contribution of solid phases is assumed to be constant and included in the equilibrium constant, thus an order of zero is assumed).

Log Power Law (log_power_law_equil)
-----------------------------------

The method uses a log form of a power law using the concentration form provided  to calculate the reaction rate.

.. math:: log(k_{eq}) = \sum_{(p, j)}{O_{(p,j)}*log(x_{(p,j))}

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Indices", "Description"

   ":math:`O`", "reaction_order", "phase, component", "Reaction order"

Providing a `reaction_order` dict is optional. If one is not provided, it will be assumed that this is an elementary reaction and that the reaction order is equal to the stoichiometric coefficient for all component in non-solid phases (the contribution of solid phases is assumed to be constant and included in the equilibrium constant, thus an order of zero is assumed).
