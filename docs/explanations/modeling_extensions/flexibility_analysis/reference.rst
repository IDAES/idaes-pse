Flexibility Analysis Reference
==============================

Table of Contents
-----------------
Enumerations

* :py:enum:`FlexTestMethod<idaes.apps.flexibility_analysis.FlexTestMethod>`
* :py:enum:`FlexTestTermination<idaes.apps.flexibility_analysis.FlexTestTermination>`
* :py:enum:`SamplingStrategy<idaes.apps.flexibility_analysis.SamplingStrategy>`
* :py:enum:`SamplingIniStrategy<idaes.apps.flexibility_analysis.SamplingInitStrategy>`

Configuration

* :class:`FlexTestConfig<idaes.apps.flexibility_analysis.FlexTestConfig>`
* :class:`ActiveConstraintConfig<idaes.apps.flexibility_analysis.ActiveConstraintConfig>`
* :class:`SamplingConfig<idaes.apps.flexibility_analysis.SamplingConfig>`
* :class:`LinearDRConfig<idaes.apps.flexibility_analysis.LinearDRConfig>`
* :class:`ReluDRConfig<idaes.apps.flexibility_analysis.ReluDRConfig>`

Results

* :class:`FlexTestResults<idaes.apps.flexibility_analysis.FlexTestResults>`

Functions

* :meth:`solve_flextest<idaes.apps.flexibility_analysis.solve_flextest>`
* :meth:`solve_flex_index<idaes.apps.flexibility_analysis.solve_flex_index>`

Enumerations
------------

.. autoenum:: idaes.apps.flexibility_analysis.FlexTestMethod

.. autoenum:: idaes.apps.flexibility_analysis.FlexTestTermination

.. autoenum:: idaes.apps.flexibility_analysis.SamplingStrategy

.. autoenum:: idaes.apps.flexibility_analysis.SamplingInitStrategy

Configuration
-------------

.. autoclass:: idaes.apps.flexibility_analysis.FlexTestConfig

.. autoclass:: idaes.apps.flexibility_analysis.ActiveConstraintConfig

.. autoclass:: idaes.apps.flexibility_analysis.SamplingConfig

.. autoclass:: idaes.apps.flexibility_analysis.decision_rules.dr_config.DRConfig

.. autoclass:: idaes.apps.flexibility_analysis.LinearDRConfig

.. autoclass:: idaes.apps.flexibility_analysis.ReluDRConfig


Results
-------

.. autoclass:: idaes.apps.flexibility_analysis.FlexTestResults

Functions
---------

.. autofunction:: idaes.apps.flexibility_analysis.solve_flextest

.. autofunction:: idaes.apps.flexibility_analysis.solve_flex_index
