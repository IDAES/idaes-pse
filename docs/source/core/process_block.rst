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


.. module:: idaes.core.process_block

.. autofunction:: declare_process_block_class

.. autoclass:: ProcessBlock
    :members:

.. module:: idaes.core.process_base

.. autoclass:: ProcessBlockData
    :members:
