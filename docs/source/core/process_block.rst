Process Blocks
==============

ProcessBlock is used to simplify inheritance of Pyomo's Block.  The code below
provides an example of how a new ProcessBlock class can be implemented.

.. code-block:: python

  from idaes.core import ProcessBlockData, decalre_process_block_class

  @decalre_process_block_class("MyBlock", doc="My new block.")
  class MyBlockData(ProcessBlockData):
      def build():
          super(MyBlockData, self).build()
          # Add variables and constraints to the block


The build method
----------------

The core part of any IDAES Block is the build method, which contains the instructions on how to construct the variables, constraints and other components that make up the model. The build method serves as the default rule for constructing an instance of an IDAES Block, and is triggered automatically whenever an instance of an IDAES Block is created unless a custom rule is provided by the user.

ProcessBlock Class
------------------

.. module:: idaes.core.process_block

.. autofunction:: declare_process_block_class

.. autoclass:: ProcessBlock
    :members:

.. module:: idaes.core.process_base

.. autoclass:: ProcessBlockData
    :members:
