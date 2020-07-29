Rate-Based Reaction Forms
=========================

.. contents:: Contents 
    :depth: 2

Power Law (power_law_rate)
--------------------------

The method uses a power law form using the concentration form provided to calculate the reaction rate.

.. math:: r_{rxn} = k_{rxn} \prod_{(p, j)}{C_{(p,j)}^{O_{(p,j)}}}

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Indices", "Description"

   ":math:`O`", "reaction_order", "phase, component", "Reaction order"

Providing a `reaction_order` dict is optional. If one is not provided, it will be assumed that this is an elementary reaction and that the reaction order is equal to the stoichiometric coefficient for the products (i.e. for all phase-component pairs with a *negative* stoichiometric coefficient, the reaction order is equal to the absolute value of the stoichiometric coefficient).
