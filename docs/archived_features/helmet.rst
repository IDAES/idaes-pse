HELMET: HELMholtz Energy Thermodynamics
========================================

.. warning::
  HELMET was removed from the IDAES codebase in v2.1 due to a lack of maintenance (3/9/2023, commit hash: 653fdcc). HELMET can be found in older releases if required, and maybe restored if sufficient interest is shown and a maintainer can be identified (see https://github.com/IDAES/idaes-pse/wiki/Maintenance-and-Testing-of-User-Supported-Code)


The purpose of HELMET (HELMholtz Energy Thermodynamics) is to provide a framework for regressing multiparameter equations of state that identify an equation for Helmholtz energy and multiple thermodynamic properties simultaneously. HELMET uses best subset selection to simultaneously model various thermodynamic properties based on the properties thermodynamic relation to Helmholtz energy. The generated model is a function of reduced density and inverse reduced temperature and uses partial derivatives to calculate the different properties. Constraints are placed on the regression to maintain thermodynamically feasible values and improve extrapolation and behavior of the model based on physical restrictions.

