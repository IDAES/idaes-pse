Single Control Volume Initializer
=================================

The Single Control Volume Initializer is a hierarchical ``Initializer`` suitable for any unit model that involves a single control volume (either 0- or 1-dimensional control volume). This ``Initializer`` is the default for all unit models in IDAES unless the unit model overloads this with a different option.

This routine applies a hierarchical approach to initialize the `StateBlock(s)` within the model followed by calling a solver to converge the full model. The following steps are used to initialize the model:

1. Determine if the control volume is 0- or 1-dimensional.

    a. For 0-D control volumes, initialize the inlet ``StateBlock`` using the appropriate ``Initializer`` followed by the outlet ``StateBlock``. The outlet ``StateBlock`` may be initialized using an ``Initializer`` or by copying the state of the inlet ``StateBlock``.
    b. For 1-D control volumes, initialize the indexed ``StateBlock`` using the appropriate ``Initializer``.

2. If present, initialize the ``ReactionBlock`` using the appropriate ``Initializer``.
3. Call user-specified solver to converge the full model.

This ``Initializer`` also supports plug-ins to the main model (sub-models which are attached to the main model after construction, such as costing calculations). If plug-ins are present, they are initialized according to the default sequence which is as follows:

1. After fixing the degrees of freedom by before initializing the main model, iterate through ``model.initialization_order`` and call ``SubmodelInitializer.plugin_prepare`` for all plug-ins found.
2. Iterate through ``model.initialization_order`` again and either call the main model initialization routine or `SubmodelInitializer.plugin_initialize`` as appropriate.
3. Call user specified solver to converge the full model with all plug-ins.
4. iterate through ``model.initialization_order`` in reverse and either call ``SubmodelInitializer.plugin_finalize`` for all plug-ins.

Users can specify ``Initializers`` for each sub-model (``StateBlocks`` and plug-ins); if users do not specify an ``Initializer`` for a sub-model then the routine will default to either i) the default ``Initializer`` specified when the main model ``Initializer`` was instantiated or ii) the default ``Initializer`` for the sub-model (specified in ``submodel.default_initializer``).

.. module:: idaes.core.initialization.general_hierarchical

SingleControlVolumeUnitInitializer Class
----------------------------------------

.. autoclass:: SingleControlVolumeUnitInitializer
  :members: precheck, initialize
