Phase Class
===========

.. module:: idaes.core.base.phases

.. autoclass:: Phase
    :members:

Phase Type Enum
---------------

In some cases, it is useful to be able to indicate a given type of phase, rather than an instance specific `Phase` class; an example would be indicating the set of valid phases for a given chemical species. In these cases, the `PhaseType` `Enum` can be used, which enumerates the different types of phases recognized by the IDAES framework.

The `PhaseType` `Enum` has the following possible values:

* `liquidPhase` (1)
* `vaporPhase` (2)
* `solidPhase` (3)
