Equilibrium Reaction Forms
==========================

.. contents:: Contents 
    :depth: 2

Power Law (power_law_equil)
---------------------------

The method uses a power law form using the concentration form provided to calculate the reaction equilibrium.

.. math:: k_{eq} = \prod_{(p, j)}{x_{(p,j)}^{O_{(p,j)}}}

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Indices", "Description"

   ":math:`O`", "reaction_order", "phase, component", "Reaction order"

Providing a `reaction_order` dict is optional. If one is not provided, it will be assumed that this is an elementary reaction and that the reaction order is equal to the stoichiometric coefficient for all components in non-solid phases (the contribution of solid phases is assumed to be constant and included in the equilibrium constant, thus an order of zero is assumed).

Log Power Law (log_power_law_equil)
-----------------------------------

The method uses a log form of a power law using the concentration form provided to calculate the reaction equilibrium.

.. math:: log(k_{eq}) = \sum_{(p, j)}{O_{(p,j)}*log(x_{(p,j)})}

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Indices", "Description"

   ":math:`O`", "reaction_order", "phase, component", "Reaction order"

Providing a `reaction_order` dict is optional. If one is not provided, it will be assumed that this is an elementary reaction and that the reaction order is equal to the stoichiometric coefficient for all components in non-solid phases (the contribution of solid phases is assumed to be constant and included in the equilibrium constant, thus an order of zero is assumed).

Solubility Product (solubility_product)
---------------------------------------

The method uses a complementarity formulation for solid precipitation using solubility product form.

For precipitation at equilibrium, 

.. math:: Q = k_{sp} - \prod_{(liq, j)}{x_{(liq,j)}^{O_{(liq,j)}}}

* :math:`Q = 0` only holds if solids (S) are present in the system, :math:`S \geq 0`.
* :math:`Q \geq 0` if the system is subsaturated and no solids (S) are present, :math:`S = 0`.

:math:`S` represents the solid presentation and is normalized and scaled as:

.. math:: S = C*\frac{s}{s+N}

where :math:`s` is assumed to be the sum of the flowrates of any solids formed in the reaction. 

Only one of :math:`Q` and :math:`S` can be greater than zero at any time, which can be written in the form of a complementarity constraint as:

.. math:: Q - max(0, Q-S) == 0

The :math:`\max()` function is provided as an IDAES utility which provides a smooth max expression with mutable smoothing parameter, :math:`\epsilon`.

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Indices", "Default", "Description"

   ":math:`O`", "reaction_order", "liq, component", "\-", "Reaction order"    
    ":math:`C`", "s_scale", "\-", "10", "Scaling factor for the solid term"    
    ":math:`N`", "s_norm", "\-", "1e-4 (`k_eq_ref` if provided)", "Normalization parameter for the solid term" 
    ":math:`\epsilon`", "eps", "\-", "1e-4", "Smoothing factor"

Providing a `reaction_order` dict is optional. If one is not provided, it will be assumed that this is an elementary reaction and that the reaction order is equal to the stoichiometric coefficient for all components in non-solid phases (the contribution of solid phases is assumed to be constant and included in the equilibrium constant, thus an order of zero is assumed).

Log Solubility Product (log_solubility_product)
-----------------------------------------------

The method uses a complementarity formulation for solid precipitation using a log form of solubility product.

For precipitation at equilibrium, 

.. math:: Q = log(k_{sp}) - \sum_{(liq, j)}{O_{(liq,j)}*log(x_{(liq,j)})}

* :math:`Q = 0` only holds if solids (S) are present in the system, :math:`S \geq 0`.
* :math:`Q \geq 0` if the system is subsaturated and no solids (S) are present, :math:`S = 0`.

:math:`S` represents the solid presentation and is normalized and scaled as:

.. math:: S = C\frac{s}{s+N}

where :math:`s` is assumed to be the sum of the flowrates of any solids formed in the reaction. 

Only one of :math:`Q` and :math:`S` can be greater than zero at any time, which can be written in the form of a complementarity constraint as:

.. math:: Q - max(0, Q-S) == 0

The :math:`\max()` function is provided as an IDAES utility which provides a smooth max expression with mutable smoothing parameter, :math:`\epsilon`.

**Parameters**

.. csv-table::
   :header: "Symbol", "Parameter Name", "Indices", "Default", "Description"

   ":math:`O`", "reaction_order", "liq, component", "\-", "Reaction order"    
    ":math:`C`", "s_scale", "\-", "1", "Scaling factor for the solid term"    
    ":math:`N`", "s_norm", "\-", "1e-4 (`k_eq_ref` if provided)", "Normalization parameter for the solid term" 
    ":math:`\epsilon`", "eps", "\-", "1e-4", "Smoothing factor"

Providing a `reaction_order` dict is optional. If one is not provided, it will be assumed that this is an elementary reaction and that the reaction order is equal to the stoichiometric coefficient for all component in non-solid phases (the contribution of solid phases is assumed to be constant and included in the equilibrium constant, thus an order of zero is assumed).

