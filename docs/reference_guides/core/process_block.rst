Process Block
=============

Example
-------

ProcessBlock is used to simplify inheritance of Pyomo's Block.  The code below
provides an example of how a new ProcessBlock class can be implemented. The
new ProcessBlock class has a ConfigBlock that allows each element of the block to
be passed configuration options that affect how a block is built. ProcessBlocks
have a rule set by default that calls the build method of the contained
ProcessBlockData class.

.. testcode::

  from pyomo.environ import *
  from pyomo.common.config import ConfigValue
  from idaes.core import ProcessBlockData, declare_process_block_class

  @declare_process_block_class("MyBlock")
  class MyBlockData(ProcessBlockData):
      CONFIG = ProcessBlockData.CONFIG()
      CONFIG.declare("xinit", ConfigValue(default=1001, domain=float))
      CONFIG.declare("yinit", ConfigValue(default=1002, domain=float))
      def build(self):
          super(MyBlockData, self).build()
          self.x = Var(initialize=self.config.xinit)
          self.y = Var(initialize=self.config.yinit)

The following example demonstrates creating a scalar instance of the new class.
The ``default`` key word argument is used to pass information on the the
MyBlockData ConfigBlock.

.. testcode::

    m = ConcreteModel()
    m.b = MyBlock(default={"xinit":1, "yinit":2})

The next example creates an indexed MyBlock instance.  In this case, each block is
configured the same, using the ``default`` argument.

.. testcode:: 

    m = ConcreteModel()
    m.b = MyBlock([0,1,2,3,4], default={"xinit":1, "yinit":2})

The next example uses the ``initialize`` argument to override the configuration of
the first block. Initialize is a dictionary of dictionaries where the key of the
top level dictionary is the block index and the second level dictionary is
arguments for the config block.

.. testcode:: 

    m = ConcreteModel()
    m.b = MyBlock([0,1,2,3,4], default={"xinit":1, "yinit":2},
                  initialize={0:{"xinit":1, "yinit":2}})

The next example shows a more complicated configuration where there are three
configurations, one for the first block, one for the last block, and one for the
interior blocks.  This is accomplished by providing the ``idx_map`` argument to
MyBlock, which is a function that maps a block index to a index in the initialize
dictionary.  In this case 0 is mapped to 0, 4 is mapped to 4, and all elements
between 0 and 4 are mapped to 1.  A lambda function is used to convert the block
index to the correct index in initialize.

.. testcode:: 

    m = ConcreteModel()
    m.b = MyBlock(
        [0,1,2,3,4],
        idx_map = lambda i: 1 if i > 0 and i < 4 else i,
        initialize={0:{"xinit":2001, "yinit":2002},
                    1:{"xinit":5001, "yinit":5002},
                    4:{"xinit":7001, "yinit":7002}})

The build method
----------------

The core part of any IDAES Block is the build method, which contains the instructions on how to construct the variables, constraints and other components that make up the model. The build method serves as the default rule for constructing an instance of an IDAES Block, and is triggered automatically whenever an instance of an IDAES Block is created unless a custom rule is provided by the user.

ProcessBlock Class
------------------

.. module:: idaes.core.base.process_block

.. autofunction:: declare_process_block_class

.. autoclass:: ProcessBlock
    :members:

.. module:: idaes.core.base.process_base

.. autoclass:: ProcessBlockData
    :members:
