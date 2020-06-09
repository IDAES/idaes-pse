Developer Standards
===================
.. contents:: Contents 
    :depth: 3


Model Formatting and General Standards
--------------------------------------
The section describes the recommended formatting used within the IDAES framework. Users are 
strongly encouraged to follow these standards in developing their models in order to improve 
readability of their code.

Headers and Meta-data
^^^^^^^^^^^^^^^^^^^^^
Model developers are encouraged to include some documentation in the header of their model 
files which provides a brief description of the purpose of the model and how it was developed. 
Some suggested information to include is:

* Model name,
* Model publication date,
* Model author
* Any necessary licensing and disclaimer information (see below).
* Any additional information the modeler feels should be included.

Coding Standard
^^^^^^^^^^^^^^^
All code developed as part of IDAES should conform to the PEP-8 standard.

Model Organization
^^^^^^^^^^^^^^^^^^
Whilst the overall IDAES modeling framework enforces a hierarchical structure on models, model 
developers are still encouraged to arrange their models in a logical fashion to aid other users 
in understanding the model. Model constraints should be grouped with similar constraints, and 
each grouping of constraints should be clearly commented. 

For property packages, it is recommended that all the equations necessary for calculating a 
given property be grouped together, clearly separated and identified by using comments.

Additionally, model developers are encouraged to consider breaking their model up into a number 
of smaller methods where this makes sense. This can facilitate modification of the code by 
allowing future users to inherit from the base model and selectively overload sub-methods where 
desired.

Commenting
^^^^^^^^^^
To help other modelers and users understand the how a model works, model builders are strongly 
encouraged to comment their code. It is suggested that every constraint should be commented 
with a description of the purpose of the constraint, and if possible/necessary a reference to a 
source or more detailed explanation. Any deviations from standard units or formatting should be 
clearly identified here. Any initialization procedures, or other procedures required to get the 
model to converge should be clearly commented and explained where they appear in the code. 
Additionally, modelers are strongly encouraged to add additional comments explaining how their 
model works to aid others in understanding the model.

Units of Measurement and Reference States
-----------------------------------------
Due to the flexibility provided by the IDAES modeling framework, there is no standard set of 
units of measurement or standard reference state that should be used in models. This places the 
onus on the user to understand the units of measurement being used within their models and to 
ensure that they are consistent.

The standard units and reference states are described in the 
:ref:`user guide<user_guide/conventions:Units of Measurement and Reference States>`.

Standard Variable Names
-----------------------
The standard variable names are described in the
:ref:`user guide<user_guide/conventions:Standard Variable Names>`.

Testing
-------
The testing standards are included :ref:`here<advanced_user_guide/developer/testing:Testing>`.

