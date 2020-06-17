Modeling Extensions
===================
The IDAES platform includes several modeling extensions that provide additional capabilities 
including surrogate modeling, material design, and control.

.. contents:: :local:

Surrogate Modeling
------------------
The surrogate modeling tools included on the IDAES platform are detailed below.

ALAMOPY: ALAMO Python
^^^^^^^^^^^^^^^^^^^^^
:ref:`ALAMOPY<user_guide/modeling_extensions/surrogate/alamopy/index:ALAMOPY: ALAMO Python>`
provides a wrapper for the software ALAMO which generates algebraic surrogate models of 
black-box systems for which a simulator or experimental setup is available.

RIPE: Reaction Identification and Parameter Estimation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:ref:`RIPE<user_guide/modeling_extensions/surrogate/ripe/index:RIPE: Reaction Identification and  Parameter Estimation>`
provides tools for reaction network identification. RIPE uses reactor data consisting of 
concentration, or conversion, values for multiple species that are obtained dynamically, or at 
multiple process conditions (temperatures, flow rates, working volumes) to identify probable 
reaction kinetics. The RIPE module also contains tools to facilitate adaptive experimental 
design.

HELMET: HELMholtz Energy Thermodynamics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:ref:`HELMET<user_guide/modeling_extensions/surrogate/helmet/index:HELMET: HELMholtz Energy Thermodynamics>`
provides a framework for regressing multiparameter equations of state that identify an equation 
for Helmholtz energy and multiple thermodynamic properties simultaneously. 

PySMO: Python-based Surrogate Modelling Objects
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:ref:`PySMO<user_guide/modeling_extensions/surrogate/pysmo/index:PySMO: Python-based Surrogate Modelling Objects>`
provides tools for generating different types of reduced order models. It provides IDAES users 
witha set of surrogate modeling tools which supports flowsheeting and direct integration into 
an equation-oriented modeling framework. It allows users to directly integrate reduced order 
models with algebraic high-fidelity process models within an single IDAES flowsheet.

.. image:: ../../_images/pysmo-logo.png
    :width: 500px
    :align: center
    
MatOpt: Nanomaterials Optimization
----------------------------------
:ref:`MatOpt<user_guide/modeling_extensions/matopt/index:MatOpt: Nanomaterials Optimization>`
provides tools for nanomaterials design using Mathematical Optimization. MatOpt can be used to 
design crystalline nanostructured materials, including but not limited to particles, wires, 
surfaces, and periodic bulk structures.

.. image:: ../../_images/matopt_logo_full.png
    :width: 500px
    :align: center

Caprese
-------
:ref:`Caprese<user_guide/modeling_extensions/caprese/index:Caprese>`
is a module for the simulation of IDAES flowsheets with nonlinear program (NLP)-based control 
and estimation strategies, namely Nonlinear Model Predictive Control (NMPC) and Moving Horizon 
Estimation (MHE).

.. image:: ../../_images/logocappresse-01.png
    :width: 500px
    :align: center


.. toctree::
    :glob:
    :hidden:
    
    */index