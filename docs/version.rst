IDAES Versioning
================

The IDAES Python package is versioned according to the general guidelines
of `semantic versioning <https://semver.org/>`_, following the recommendations of
`PEP 440 <https://www.python.org/dev/peps/pep-0440/>`_ with respect to
extended versioning descriptors (alpha, beta, release candidate, etc.).

Basic usage
-----------

You can see the version of the package at any time interactively by
printing out the `__version__` variable in the top-level package::

    import idaes
    print(idaes.__version__)
    # prints a version like "1.2.3"

.. testcode::
    :hide:
    
    # This test ensures that the above code works without needing to be
    # updated every time the IDAES version changes
    import idaes
    idaes.__version__

Advanced usage
--------------

This section describes the module's variables and classes.

Overview
^^^^^^^^

.. automodule:: idaes.ver
    :members: package_version, __version__
    :noindex:

Version class
^^^^^^^^^^^^^
The versioning semantics are encapsulated in a class called `Version`.

.. autoclass:: idaes.ver.Version
    :members: __init__, __iter__, __str__
    :noindex:

HasVersion class
^^^^^^^^^^^^^^^^
For adding versions to other classes in a simple and standard way,
you can use the `HasVersion` mixin class.

.. autoclass:: idaes.ver.HasVersion
    :members: __init__
    :noindex:
