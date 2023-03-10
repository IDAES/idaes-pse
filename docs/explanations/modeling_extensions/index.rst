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

.. rubric:: PySMO: Python-based Surrogate Modeling Objects

:ref:`PySMO: Python-based Surrogate Modeling Objects<explanations/modeling_extensions/surrogate/api/pysmo/index:PySMO: Python-based Surrogate Modeling Objects>`
provides tools for generating different types of reduced order models. It provides IDAES users 
with a set of surrogate modeling tools which supports flowsheeting and direct integration into 
an equation-oriented modeling framework. It allows users to directly integrate reduced order 
models with algebraic high-fidelity process models within an single IDAES flowsheet.

.. image:: /images/pysmo-logo.png
    :width: 500px
    :align: center

.. rubric:: OMLT: Optimization and Machine Learning Toolkit

:ref:`OMLT: Optimization and Machine Learning Toolkit<explanations/modeling_extensions/surrogate/api/omlt-keras/index:OMLT: Optimization and Machine Learning Toolkit>`
provides tools for generating machine learning, neural network and gradient-boosted tree models with the Pyomo optimization environment. It provides IDAES users optimization formulations for ML models (e.g. full-space, reduced-space, MILP) and an interface to import Sequential Keras and general ONNX models.


.. note::
    OMLT is an external tool that requires additional dependencies. See the OMLT documentation for details on installation and usage.


.. image:: /images/omlt-keras-logo.png
    :width: 500px
    :align: center

.. rubric:: ALAMOPY: ALAMO Python

:ref:`ALAMOPY: ALAMO Python<explanations/modeling_extensions/surrogate/api/alamopy/index:ALAMOPY: ALAMO Python>`
provides a wrapper for the software ALAMO which generates algebraic surrogate models of 
black-box systems for which a simulator or experimental setup is available.


.. note::
    ALAMO is an external tool to needs to be installed separately. See the ALAMOPy documentation for details on ALAMO installation and usage.


.. rubric:: MatOpt: Nanomaterials Optimization

:ref:`MatOpt: Nanomaterials Optimization<explanations/modeling_extensions/matopt/index:MatOpt: Nanomaterials Optimization>`
provides tools for nanomaterials design using Mathematical Optimization. MatOpt can be used to 
design crystalline nanostructured materials, including but not limited to particles, wires, 
surfaces, and periodic bulk structures.

.. image:: /images/matopt_logo_full.png
    :width: 500px
    :align: center

.. rubric:: Caprese

:ref:`Caprese<explanations/modeling_extensions/caprese/index:Caprese>`
is a module for the simulation of IDAES flowsheets with nonlinear program (NLP)-based control 
and estimation strategies, namely Nonlinear Model Predictive Control (NMPC) and Moving Horizon 
Estimation (MHE).

.. image:: /images/logocappresse-01.png
    :width: 500px
    :align: center

.. rubric:: Uncertainty Propagation Toolbox

:ref:`Uncertainty Propagation Toolbox<explanations/modeling_extensions/uncertainty_propagation/index:Uncertainty Propagation Toolbox>`
is a module for quantifying and propagating parametric uncertainty through an optimization or simulation problem based on an IDAES model.

.. rubric:: Degeneracy Hunter

:ref:`Degeneracy Hunter<explanations/modeling_extensions/diagnostics/index:Degeneracy Hunter>`
is coming soon!
