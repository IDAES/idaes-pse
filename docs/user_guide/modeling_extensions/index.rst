Modeling Extensions
===================
The IDAES platform includes several modeling extensions that provide additional capabilities 
including surrogate modeling, material design, and control. A brief description of each is 
provided below.

.. toctree::
    :maxdepth: 2
    
    surrogate/index

.. toctree::
    :maxdepth: 1
    
    matopt/index
    caprese/index
    uncertainty_propagation/index
    diagnostics/index

.. rubric:: ALAMOPY: ALAMO Python

:ref:`ALAMOPY<user_guide/modeling_extensions/surrogate/alamopy/index:ALAMOPY: ALAMO Python>`
provides a wrapper for the software ALAMO which generates algebraic surrogate models of 
black-box systems for which a simulator or experimental setup is available.

.. rubric:: RIPE: Reaction Identification and Parameter Estimation

:ref:`RIPE<user_guide/modeling_extensions/surrogate/ripe/index:RIPE: Reaction Identification and  Parameter Estimation>`
provides tools for reaction network identification. RIPE uses reactor data consisting of 
concentration, or conversion, values for multiple species that are obtained dynamically, or at 
multiple process conditions (temperatures, flow rates, working volumes) to identify probable 
reaction kinetics. The RIPE module also contains tools to facilitate adaptive experimental 
design.

.. rubric:: HELMET: HELMholtz Energy Thermodynamics

:ref:`HELMET<user_guide/modeling_extensions/surrogate/helmet/index:HELMET: HELMholtz Energy Thermodynamics>`
provides a framework for regressing multiparameter equations of state that identify an equation 
for Helmholtz energy and multiple thermodynamic properties simultaneously. 

.. rubric:: PySMO: Python-based Surrogate Modelling Objects

:ref:`PySMO<user_guide/modeling_extensions/surrogate/pysmo/index:PySMO: Python-based Surrogate Modelling Objects>`
provides tools for generating different types of reduced order models. It provides IDAES users 
with a set of surrogate modeling tools which supports flowsheeting and direct integration into 
an equation-oriented modeling framework. It allows users to directly integrate reduced order 
models with algebraic high-fidelity process models within an single IDAES flowsheet.

.. image:: /images/pysmo-logo.png
    :width: 500px
    :align: center
    
.. rubric:: MatOpt: Nanomaterials Optimization

:ref:`MatOpt<user_guide/modeling_extensions/matopt/index:MatOpt: Nanomaterials Optimization>`
provides tools for nanomaterials design using Mathematical Optimization. MatOpt can be used to 
design crystalline nanostructured materials, including but not limited to particles, wires, 
surfaces, and periodic bulk structures.

.. image:: /images/matopt_logo_full.png
    :width: 500px
    :align: center

.. rubric:: Caprese

:ref:`Caprese<user_guide/modeling_extensions/caprese/index:Caprese>`
is a module for the simulation of IDAES flowsheets with nonlinear program (NLP)-based control 
and estimation strategies, namely Nonlinear Model Predictive Control (NMPC) and Moving Horizon 
Estimation (MHE).

.. image:: /images/logocappresse-01.png
    :width: 500px
    :align: center

.. rubric:: Uncertainty Propagation Toolbox

:ref:`uncertainty_propagation<user_guide/modeling_extensions/uncertainty_propagation/index:Uncertainty Propagation Toolbox>`
is a module for quantifying and propagating parametric uncertainty through an optimization or simulation problem based on an IDAES model.

.. rubric:: Degeneracy Hunter

:ref:`Degeneracy Hunter<user_guide/modeling_extensions/diagnostics/index:Degeneracy Hunter>`
is coming soon!
