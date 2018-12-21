Model State Serialization
=========================

.. module:: idaes.core.util.model_serializer

The IDAES framework has some utility functions for serializing the state of a
Pyomo model. These functions can save and load attributes of Pyomo components,
but cannot reconstruct the Pyomo objects (it is not a replacement for pickle).

Below are a few example use cases for this module.

* Some models are very large and complex they may take over a minutes to
initialize.  Once a model is initialization it's state can be saved and reloaded
when running the model in the future, so the initialization procedure does not
need to be run.
* Saving simulation or optimization results is another use.  Results can be
stored for later evaluation without needed to rerun the model. These results can
be archived in a data management system if needed later.
* These functions may be useful in writing initialization procedures. For example,
a model may be constructed and ready to run but first may need to be initialized.
Which these functions, which components are active and with variables
are fixed can be stored.  The initialization can change witch variables fixed and
which components are active.  The original state can be read back after
initialization, but were only values of variables that were originally fixed are
read back in.  This is an easy way to ensure whatever the initialization
procedure may do, you can always end up with exactly the same problem (with only
better initial values for unfixed variables).
* These functions can be used to send model data and results to user interface
components.

Examples
--------

This section provides a few very simple examples of how to use these functions.

Example Models
~~~~~~~~~~~~~~

This section provides some boilerplate and functions to create a couple simple
test models.  The second model is a little more complicated and includes suffixes.

.. code-block:: python
  from pyomo.environ import *
  from idaes.core.util import to_json, from_json, StoreSpec

  def setup_model01():
      model = ConcreteModel()
      model.b = Block([1,2,3])
      a = model.b[1].a = Var(bounds=(-100, 100), initialize=2)
      b = model.b[1].b = Var(bounds=(-100, 100), initialize=20)
      model.b[1].c = Constraint(expr=b==10*a)
      a.fix(2)
      return model

  def setup_model02():
      model = ConcreteModel()
      a = model.a = Param(default=1, mutable=True)
      b = model.b = Param(default=2, mutable=True)
      c = model.c = Param(initialize=4)
      x = model.x = Var([1,2], initialize={1:1.5, 2:2.5}, bounds=(-10,10))
      model.f = Objective(expr=(x[1] - a)**2 + (x[2] - b)**2)
      model.g = Constraint(expr=x[1] + x[2] - c >= 0)
      model.dual = Suffix(direction=Suffix.IMPORT)
      model.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
      model.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
      return model

Serialization
~~~~~~~~~~~~~

These examples can be appended to the boilerplate code above.


The first example creates a model, saves the state, changes a value, then reads
back the initial state.

.. codeblock:: python
  model = setup_model01()
  to_json(model, fname="ex.json.gz", gzip=True, human_read=True)
  model.b[1].a = 3000.4
  from_json(model, fname="ex.json.gz", gzip=True)
  print(model.b[1].a)

This next example show how to save only suffixes.

.. codeblock:: python
  model = setup_model02()
  store_spec = StoreSpec.suffix()
  to_json(model, fname="ex.json", wts=store_spec)
  # Do something and now I want my suffixes back
  from_json(model, fname="ex.json.gz", wts=store_spec)

to_json
-------

Despite the name of the ``to_json`` function it is capable of creating Python
dictionaries, json files, gzipped json files, and json strings. The function
documentation is below.  A :ref:`StoreSpec <core/util/model_serializer:StoreSpec>` object provides
the function with details on what to store and how to handle special cases of
Pyomo component attributes.

.. autofunction:: to_json

from_json
---------

The ``from_json`` function puts data from Python dictionaries, json files,
gzipped json files, and json strings back into a Pyomo model. The function
documentation is below.  A :ref:`StoreSpec <core/util/model_serializer:StoreSpec>` object provides
the function with details on what to read and how to handle special cases of
Pyomo component attributes.

.. autofunction:: from_json

StoreSpec
---------

StoreSpec is a class for objects that tell the ``to_json()`` and ``from_json()``
functions how to read and write Pyomo component attributes.  The default
initialization provides an object that would load and save attributes usually
needed to save a model state.  There are several other class methods that
provided canned object for specific uses. Through initialization arguments, the
behavior is highly customizable.  See the class documentation below.

.. autoclass:: StoreSpec
    :members:
