Initializing Models
===================

Good initialization of variable values is critical for reliably solving non-linear models. To assist users with developing initial solutions for their problems, the IDAES-IP provides a number on initialization routines for different models. The routines are implemented as ``Initializer`` objects which contain an ``initialize`` method which will execute the initialization routine.

Usage
-----

.. code:: python

    import MyInitializer

    initializer = MyInitializer(**configuration)
    initializer.initialize(model, **arguments)

Available Initializers
----------------------

The following ``Initializers`` are available in the core libraries. Common use ``Initializers`` are generally applicable to any model type, while model specific ``Initializers`` are suitable for specific types of models or applications.

.. toctree::
    :maxdepth: 2

    general_initializers
    specific_initializers
    developing_initializers

