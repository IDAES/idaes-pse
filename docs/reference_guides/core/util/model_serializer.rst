Model State Serialization
=========================

.. module:: idaes.core.util.model_serializer

The IDAES framework has some utility functions for serializing the state of a
Pyomo model. These functions can save and load attributes of Pyomo components,
but cannot reconstruct the Pyomo objects (it is not a replacement for pickle).
It does have some advantages over pickle though. Not all Pyomo models are
picklable. Serialization and deserialization of the model state to/from json is
more secure in that it only deals with data and not executable code.  It should
be safe to use the ``from_json()`` function with data from untrusted sources,
while, unpickling an object from an untrusted source is not secure.  Storing a
model state using these functions is also probably more robust against Python
and Python package version changes, and possibly more suitable for long-term storage
of results.

Below are a few example use cases for this module.

* Some models are very complex and may take minutes to initialize.  Once a model is initialized it's state can be saved. For future runs, the initialized state can be reloaded instead of rerunning the initialization procedure.
* Results can be stored for later evaluation without needing to rerun the model. These results can be archived in a data management system if needed later.
* These functions may be useful in writing initialization procedures. For example, a model may be constructed and ready to run but first it may need to be initialized. Which components are active and which variables are fixed can be stored.  The initialization can change which variables are fixed and which components are active.  The original state can be read back after initialization, but where only values of variables that were originally fixed are read back in.  This is an easy way to ensure that whatever the initialization procedure may do, the result is exactly the same problem (with only better initial values for unfixed variables).
* These functions can be used to send and receive model data to/from JavaScript user interface components.

Examples
--------

This section provides a few very simple examples of how to use these functions.

Example Models
~~~~~~~~~~~~~~

This section provides some boilerplate and functions to create a couple simple
test models.  The second model is a little more complicated and includes suffixes.

.. testcode::

  from pyomo.environ import *
  from idaes.core.util import to_json, from_json, StoreSpec

  def setup_model01():
      model = ConcreteModel()
      model.b = Block([1,2,3])
      a = model.b[1].a = Var(bounds=(-1e5, 1e5), initialize=2)
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

.. testcode::

  model = setup_model01()
  to_json(model, fname="ex.json.gz", gz=True, human_read=True)
  model.b[1].a = 3000.4
  from_json(model, fname="ex.json.gz", gz=True)
  print(value(model.b[1].a))

.. testoutput::

  2

This next example show how to save only suffixes.

.. testcode::

  model = setup_model02()
  # Suffixes here are read back from solver, so to have suffix data,
  # need to solve first
  solver = SolverFactory("ipopt")
  solver.solve(model)
  store_spec = StoreSpec.suffix()
  to_json(model, fname="ex.json", wts=store_spec)
  # Do something and now I want my suffixes back
  from_json(model, fname="ex.json", wts=store_spec)

to_json
-------

Despite the name of the ``to_json`` function it is capable of creating Python
dictionaries, json files, gzipped json files, and json strings. The function
documentation is below.  A :ref:`StoreSpec <reference_guides/core/util/model_serializer:StoreSpec>`
object provides the function with details on what to store and how to handle
special cases of Pyomo component attributes.

.. autofunction:: to_json

from_json
---------

The ``from_json`` function puts data from Python dictionaries, json files,
gzipped json files, and json strings back into a Pyomo model. The function
documentation is below.  A :ref:`StoreSpec <reference_guides/core/util/model_serializer:StoreSpec>`
object provides the function with details on what to read and how to handle
special cases of Pyomo component attributes.

.. autofunction:: from_json

StoreSpec
---------

``StoreSpec`` is a class for objects that tell the ``to_json()`` and ``from_json()``
functions how to read and write Pyomo component attributes.  The default
initialization provides an object that would load and save attributes usually
needed to save a model state.  There are several other class methods that
provide canned objects for specific uses. Through initialization arguments, the
behavior is highly customizable. Attributes can be read or written using callback
functions to handle attributes that can not be directly read or written (e.g.
a variable lower bound is set by calling setlb()). See the class documentation below.

.. autoclass:: StoreSpec
    :members:

Structure
---------

Python dictionaries, json strings, or json files are generated, in any case the
structure of the data is the same.  The current data structure version is 3.

The example json below shows the top-level structure.  The
``"top_level_component"`` would be the name of the Pyomo component that is being
serialized. The top level component is the only place were the component name does
not matter when reading the serialized data.

.. testcode::

  {
      "__metadata__": {
          "format_version": 3,
          "date": "2018-12-21",
          "time": "11:34:39.714323",
          "other": {
          },
          "__performance__": {
              "n_components": 219,
              "etime_make_dict": 0.003}
      },
      "top_level_component":{
        "...": "..."
      },
  }

The data structure of a Pyomo component is shown below.  Here ``"attribute_1"``
and ``"attribute_2"`` are just examples the actual attributes saved depend on
the "wts" argument to ``to_json()``. Scalar and indexed components have the
same structure. Scalar components have one entry in ``"data"`` with an index of
``"None"``.  Only components derived from Pyomo's ``_BlockData``
have a ``"__pyomo_components__"`` field, and components appearing there are keyed
by their name. The data structure duplicates the hierarchical structure of the
Pyomo model.

Suffixes store extra attributes for Pyomo components that are not stored on the
components themselves. Suffixes are a Pyomo structure that comes from the AMPL
solver interface.  If a component is a suffix, keys in the data section are the
serial integer component IDs generated by ``to_json()``, and the value is the
value of the suffix for the corresponding component.

.. testcode::

    {
        "__type__": "<class 'some.class'>",
        "__id__": 0,
        "data":{
          "index_1":{
              "__type__":"<usually a component class but for params could be float, int, ...>",
              "__id__": 1,
              "__pyomo_components__":{
                "child_component_1": {
                  "...": "..."
                }
              },
              "attribute_1": "... could be any number of attributes like 'value': 1.0,",
              "attribute_2": "..."
          }
        },
        "attribute_1": "... could be any number of attributes like 'active': true,",
        "attribute_2": "..."
    }

As a more concrete example, here is the json generated for example model 2 in
:ref:`Examples <reference_guides/core/util/model_serializer:Examples>`.
This code can be appended to the
:ref:`example boilerplate above <reference_guides/core/util/model_serializer:Examples>`.
To generate the example json shown.

.. testcode::

  model = setup_model02()
  solver = SolverFactory("ipopt")
  solver.solve(model)
  to_json(model, fname="ex.json")

The resulting json is shown below. The top-level component
in this case is given as "unknown," because the model was not given a name. The
top level object name is not needed when reading back data, since the top level
object is specified in the call to ``from_json()``.  Types are not used when
reading back data, they may have some future application, but at this point they
just provide a little extra information.

.. code-block:: json

  {
    "__metadata__":{
      "format_version":3,
      "date":"2019-01-02",
      "time":"10:22:25.833501",
      "other":{
      },
      "__performance__":{
        "n_components":18,
        "etime_make_dict":0.000955581665039062
      }
    },
    "unknown":{
      "__type__":"<class 'pyomo.core.base.PyomoModel.ConcreteModel'>",
      "__id__":0,
      "active":true,
      "data":{
        "None":{
          "__type__":"<class 'pyomo.core.base.PyomoModel.ConcreteModel'>",
          "__id__":1,
          "active":true,
          "__pyomo_components__":{
            "a":{
              "__type__":"<class 'pyomo.core.base.param.SimpleParam'>",
              "__id__":2,
              "_mutable":true,
              "data":{
                "None":{
                  "__type__":"<class 'pyomo.core.base.param.SimpleParam'>",
                  "__id__":3,
                  "value":1
                }
              }
            },
            "b":{
              "__type__":"<class 'pyomo.core.base.param.SimpleParam'>",
              "__id__":4,
              "_mutable":true,
              "data":{
                "None":{
                  "__type__":"<class 'pyomo.core.base.param.SimpleParam'>",
                  "__id__":5,
                  "value":2
                }
              }
            },
            "c":{
              "__type__":"<class 'pyomo.core.base.param.SimpleParam'>",
              "__id__":6,
              "_mutable":false,
              "data":{
                "None":{
                  "__type__":"<class 'pyomo.core.base.param.SimpleParam'>",
                  "__id__":7,
                  "value":4
                }
              }
            },
            "x":{
              "__type__":"<class 'pyomo.core.base.var.IndexedVar'>",
              "__id__":8,
              "data":{
                "1":{
                  "__type__":"<class 'pyomo.core.base.var._GeneralVarData'>",
                  "__id__":9,
                  "fixed":false,
                  "stale":false,
                  "value":1.5,
                  "lb":-10,
                  "ub":10
                },
                "2":{
                  "__type__":"<class 'pyomo.core.base.var._GeneralVarData'>",
                  "__id__":10,
                  "fixed":false,
                  "stale":false,
                  "value":2.5,
                  "lb":-10,
                  "ub":10
                }
              }
            },
            "f":{
              "__type__":"<class 'pyomo.core.base.objective.SimpleObjective'>",
              "__id__":11,
              "active":true,
              "data":{
                "None":{"__type__":"<class 'pyomo.core.base.objective.SimpleObjective'>",
                "__id__":12,
                "active":true
                }
              }
            },
            "g":{
              "__type__":"<class 'pyomo.core.base.constraint.SimpleConstraint'>",
              "__id__":13,
              "active":true,
              "data":{
                "None":{
                  "__type__":"<class 'pyomo.core.base.constraint.SimpleConstraint'>",
                  "__id__":14,
                  "active":true
                }
              }
            },
            "dual":{
              "__type__":"<class 'pyomo.core.base.suffix.Suffix'>",
              "__id__":15,
              "active":true,
              "data":{
                "14":0.9999999626149493
              }
            },
            "ipopt_zL_out":{
              "__type__":"<class 'pyomo.core.base.suffix.Suffix'>",
              "__id__":16,
              "active":true,
              "data":{
                "9":2.1791814146763388e-10,
                "10":2.004834508495852e-10
              }
            },
            "ipopt_zU_out":{
              "__type__":"<class 'pyomo.core.base.suffix.Suffix'>",
              "__id__":17,
              "active":true,
              "data":{
                "9":-2.947875485096964e-10,
                "10":-3.3408951850535573e-10
              }
            }
          }
        }
      }
    }
  }
