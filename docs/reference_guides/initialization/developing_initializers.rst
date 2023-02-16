Developing Custom Initializers
==============================

To assist users with developing custom ``Initializer`` objects, the IDAES-IP provides two base classes which define the standard API for ``Initializer`` objects and methods for performing common activities.

.. contents:: Contents
    :depth: 2

The ``InitializerBase`` class defines the standard API for ``Initializer`` objects, whilst the ``ModularInitializerBase`` class extends this with some methods useful for defining hierarchical initialization routines.

Standard API
------------

All ``Initializer`` objects are expected to define an ``initialize`` method which will be called in order to run the initialization routine. This method may in turn call other supporting methods, and the base classes provide pre-defined methods for common tasks such as fixing and restoring degrees of freedom and performing pre- and post-initialization checks. Custom routines are not required to make use of these methods, and may overload these with custom methods if desired.

``Initializers`` intended for “plug-in” type models (i.e., models attached to other models after the parent model has been constructed) and used as part of  hierarchical initialization routine should also implement the ``plugin_prepare``, ``plugin_initialize`` and ``plugin_finalize`` methods as necessary.

.. module:: idaes.core.initialization.initializer_base

InitializerBase Class
---------------------

.. autoclass:: InitializerBase
  :members:

ModularInitializerBase Class
----------------------------

.. autoclass:: ModularInitializerBase
  :members:
