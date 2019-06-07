Vapor-Liquid Equilibrium Property Models (Ideal Gas - Non-ideal Liquids)
=========

This property package supports phase equilibrium calucations with a smooth phase transition formulation that makes it amenable for equation oriented optimization. The gas phase is assumed to be ideal and for the liquid phase,
the package supports an ideal liquid or a non-ideal liquid using an activity
coefficient model. To compute the activity coefficient, the package currently supports the Non Random Two Liquid Model (NRTL) or the
Wilson model. Therefore, this property package supports the following combinations for gas-liquid mixtures for VLE calculations:

1. Ideal - Ideal
2. Ideal - NRTL
3. Ideal - Wilson

**Flow basis**: Molar

**Units**: SI units

**State Variables**: Total molar flow rate, temperature, pressure and mixture mole fraction

Support for other combinations of state variables will be made available in the future.



Inputs
------

When instantiating the parameter block that uses this particular state block, 2 arguments can be passed:

- "valid_phase" - Valid options: "Liq" or "Vap" or ("Liq", "Vap") or ("Vap", "Liq")
- "activity_coeff_model" - Valid options: "Ideal" or "NRTL" or "Wilson"

The 'valid_phase' argument denotes the phases a user expects for a given set of state variables. For example, if the user knows with certainty that the valid phase will be only Liquid or Vapor, then it is best not include the complex flash equilibrium constraints in the model. If the user does not specify any option, then the package defaults to a 2 phase assumption meaning that the phase equilibrium will be computed. 

The 'activity_coeff_model' denotes the liquid phase assumption to be used. If the user does not specify any option, then the package defaults to asuming an ideal liquid assumption. 


Degrees of Freedom
------------------
The number of degrees of freedom that need to be fixed to yield a square problem (i.e. degrees of freedom = 0) depends on the options selected. If an ideal-ideal option is selected, the only variables that need to be fixed are the state variables mentioned above and this totals to 3 + number of components.

If ideal-NRTL is chosen then in addition to the state variables mentioned above, the user will have to also fix the following variables required for the NRTL activity coefficient model:

- alpha[i, j] - (number of components)^2
- tau[i, j] - (number of components)^2 (Please refer to reference 3 for recommended values)

If ideal-Wilson is chosen then in addition to the state variables mentioned above, the user will have to also fix the following variables required for the Wilson coefficient model:

- alpha[i, j] - (number of components)^2
- tau[i, j] - (number of components)^2 (Please refer to reference 3 for recommended values)






VLE Model with Smooth Phase Transition
--------------------------------------





Activity Coefficient Model - NRTL
-----------------------------------




Activity Coefficient Model - Wilson
-----------------------------------



Initialization
--------------


References:

1. "The properties of gases and liquids by Robert C. Reid"
2. "Perry's Chemical Engineers Handbook by Robert H. Perry".
3. H. Renon and J.M. Prausnitz, "Local compositions in thermodynamic excess
   functions for liquid mixtures.", AIChE Journal Vol. 14, No.1, 1968.
4. AP Burgard, JP Eason, JC Eslick, JH Ghouse, A Lee, LT Biegler, DC Miller. "A Smooth, Square Flash Formulation for Equation Oriented Flowsheet Optimization",
   Computer Aided Chemical Engineering 44, 871-876, 2018


