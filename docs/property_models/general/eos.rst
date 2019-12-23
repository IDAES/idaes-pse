Equations of State
==================

Equations of State (or equivalent methods) describe the relationship between different thermophysical properties and ensure that the behavior of these are thermodynamically consistent. A wide range of equations of state have been develop for different applications and levels of rigor. Equations of state generally start with ideal pure component properties, and provide a set of relationships which describe how these are combined and deviate from ideality in real mixtures. Equation of state packages within the IDAES Generic Property Package Framework need to implement equations (either `Constraints` or `Expressions`) for all of the mixture properties of interest to the user relating these to the pure component properties and state variables.

The IDAES Generic Property Package Framework provides a number of prebuilt equation of state packages for users to use, which are listed below.

Equation of State Libraries
---------------------------

.. toctree::
    :maxdepth: 1

    eos/ideal
