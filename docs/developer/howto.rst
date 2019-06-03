
.. _idaes-contrib-guide:

IDAES contributor guide
========================

About
-----
This page tries to give all the essential information needed
to contribute software to the IDAES project. It is designed
to be useful to both internal and external collaborators.

Code and other file locations
-----------------------------
Source code
    The main Python package is under the `idaes/` directory.
    Sub-directories, aka subpackages, should be documented elsewhere.
    If you add a new directory in this tree, be sure to add a `__init__.py` in that directory
    so Python knows it is a subpackage with Python modules.
    Code that is not part of the core package is under `apps/`. This code can have any
    layout that the creator wants.

Documentation
    The documentation for the core package is under `docs`.
    The documentation for the `apps/` directory is not (currently) being built automatically.

Examples
    Examples are under the `examples/` directory.
    Tutorials from workshops are under the `examples/workshops/` subdirectory.
    

Code style
------------
The code style is not entirely consistent. But some general guidelines are:

* follow the `PEP8`_ style (or variants such as `Black`_)
* use `Google-style`_ docstrings on classes, methods, and functions
* format your docstrings as `reStructuredText`_ so they can be nicely rendered as HTML by Sphinx
* add logging to your code by creating and using a global log object named
  for the module, which can be created like: ``_log = logging.getLogger(__name__)``
* take credit by adding a global author variable: ``__author__ = 'yourname'``

.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _Black: https://github.com/python/black
.. _Google-style: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
.. _reStructuredText: http://docutils.sourceforge.net/rst.html

Tests
-----
Coming soon.


Documentation
--------------
The documentation is built from its sources with a tool called Sphinx.
The sources for the documentation are:

* hand-written text files, under `docs/`, with the extension ".rst" for `reStructuredText`_.
* the Python source code
* selected Jupyter Notebooks 

Building documentation
^^^^^^^^^^^^^^^^^^^^^^
To build the documentation locally, there is a "Makefile" in the `docs/` directory::

    cd docs
    make allclean
    make all

The above commands will do a completely clean build to create HTML output.
They will also attempt to execute the tutorials. During development, more
specific Makefile targets may save time:

``make html``
    Only build the HTML from the existing `.rst` files and generated API docs.
    Does not rebuild the tutorials or regenerate the API docs.

``make apidoc``
    Just regenerate API documentation source from the Python code. This does
    not change the HTML output.

``make tutorials``
    Generate HTML web pages from the Jupyter Notebook tutorials

Like any other Makefile, you can use these targets together.
So, if you are editing source code and want to preview the generated documentation,
you should run: ``make apidoc html``. This will regenerate `.rst` files from the
source code, then build those files together with hand-edited files into the
HTML output.

Previewing documentation
^^^^^^^^^^^^^^^^^^^^^^^^
The generated documentation can be previewed locally by opening
the generated HTML files in a web browser. The files are under the `docs/build/`
directory, so you can open the file ``docs/build/index.html`` to get started.
