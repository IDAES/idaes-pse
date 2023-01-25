Initialize From Data
====================

This ``Initializer`` is intended to initialize a model from a user defined data source. This data source may take the form of a pre-saved json file or a user provided dict-like data structure. This ``Initializer`` loads variable values from the data provided and checks that the resulting values satisfy all constraints in the model.

.. note::

  When initializing from a json format, only the values of unfixed variables will be loaded. Any other information stored in the json format (e.g., whether a variable is fixed or any information regarding constraints) will be ignored. 

.. module:: idaes.core.initialization.initialize_from_data

FromDataInitializer Class
-------------------------

.. autoclass:: FromDataInitializer
  :members: precheck, initialize
