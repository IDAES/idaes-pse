Serialize Pyomo Model States
================================

This describes the IDAES Pyomo model state JSON serializer module
(idaes_models.core.util.model_serializer).  All objects inheriting from
idaes_models.core.process_base.ProcessBlock (which inherits Pyomo Block) have a
from_json and to_json method that uses the functions in this module, passing
self for o, and the rest of the arguments remaining the same.

This module can load/save the model state from/to a JSON file or an in memory as
a Python dictionary.  The dictionary is in a form that can be dumped to JSON.

The following describes the important functions and classes in this module.

.. module:: idaes_models.core.util.model_serializer

.. autofunction:: to_json

.. autofunction:: from_json

.. autoclass:: StoreSpec
    :members:
