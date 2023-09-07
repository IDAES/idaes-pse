Surrogate Modeling
==================

.. toctree::
    :maxdepth: 1

    api/index
    sampling/index
    plotting/index

IDAES contains several surrogate modeling tools, including the IDAES Surrogates API
which enables integrating ALAMO, PySMO or Keras surrogate models into IDAES
flowsheets.

ALAMOpy is a data driven machine learning (ddm-learning) tool
which leverage ALAMO, a regression application for generating equation-based
surrogates. The provided tools allow access to ALAMO and other solvers through the Python
API, and an interface for accessing ALAMO solvers (via an external license)
within IDAES flowsheets.

.. image:: /images/pysmo-logo.png
    :width: 300px
    :align: center

Python-based Surrogate Modeling Objects (PySMO) is a framework for general-purpose
surrogate modeling techniques, integrated within the Pyomo mathematical optimization
framework (on which IDAES is also based). The provided tools include an interface for
accessing PySMO methods (via IDAES-internal scripts) within IDAES flowsheets.

.. image:: /images/omlt-keras-logo.png
    :width: 300px
    :align: center

Keras is a deep learning framework that integrates with TensorFlow's structure for
building and training artificial neural networks, and minimizes the number of user
actions required to construct accurate networks. OMLT (Optimization and Machine
Learning Toolkit) provides an interface to formulate machine learning models and
import Keras or ONNX models as Pyomo blocks. The provided tools include an interface
for accessing the Keras module (via the publicly available Python package) within
IDAES flowsheets.
