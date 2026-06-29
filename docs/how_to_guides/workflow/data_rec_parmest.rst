Data Reconciliation and Parameter Estimation
============================================

This workflow generally describes features of the IDAES framework that are useful
for data reconciliation and parameter estimation.  Many of these features can
be used for any task where plant data is to be used in conjunction with
a process model. It is assumed that the user is familiar with the IDAES modeling
platform and Pyomo. See the :ref:`General Workflow <how_to_guides/workflow/general:General Workflow>`
for more information on how to set up a model.

This provides general information about IDAES functionality laid out in
terms of typical use cases. See
:ref:`﻿Tutorials and Examples<tutorials/tutorials_examples:Tutorials and Examples>`, for
specific complete examples of data reconciliation and parameter estimation examples.
Relevant tutorials can be found in ``tutorials/advanced/data_recon_and_parameter_estimation``.

.. warning::

   A previous version of this document contained additional sections referring to the Data Management Framework (DMF).
   Starting with IDAES 2.6, the DMF is no lonnger supported, and those sections have been removed.

Parameter Estimation
--------------------

For parameter estimation and data reconciliation, it is recommended to use
`Paramest <https://pyomo.readthedocs.io/en/stable/contributed_packages/parmest/>`_.
If more control is needed a user can also set up their own parameter estimation
problem, by combining multiple process models into a larger parameter estimation
model.
