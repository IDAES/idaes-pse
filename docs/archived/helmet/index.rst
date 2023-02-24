HELMET: HELMholtz Energy Thermodynamics
========================================

.. warning::

    HELMET was removed from the IDAES codebase in v2.1 (24th February 2023) due to lack of support an maintenance. Git hash: 717deec3681b87752e9fa32819f251f634d0f9e8

The purpose of HELMET (HELMholtz Energy Thermodynamics) is to provide a framework for regressing multiparameter equations of state that identify an equation for Helmholtz energy and multiple thermodynamic properties simultaneously. HELMET uses best subset selection to simultaneously model various thermodynamic properties based on the properties thermodynamic relation to Helmholtz energy. The generated model is a function of reduced density and inverse reduced temperature and uses partial derivatives to calculate the different properties. Constraints are placed on the regression to maintain thermodynamically feasible values and improve extrapolation and behavior of the model based on physical restrictions.

