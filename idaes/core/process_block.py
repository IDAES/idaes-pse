##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
The process_block module simplifies inheritance of Pyomo blocks. The main reason
to subclass a Pyomo block is to create a block that comes with pre-defined model
equations. This is used in the IDAES modeling framework to to create modular
process model blocks.

"""
from __future__ import absolute_import, division, print_function

import sys
import logging
from pyomo.environ import Block

__author__ = "John Eslick"
__all__ = ['ProcessBlock', 'declare_process_block_class']

def _rule_default(b, *args):
    """Default rule for ProcessBlock, which calls build(). A different rule can
    be specified to add additional build steps, or to not call build at all
    using the normal rule argument to ProcessBlock init.
    """
    try:
        b.build()
    except Exception as e:
        logging.getLogger(__name__).exception(
            "Failure in build: {}".format(b))
        raise e

_process_block_docstring = """
Args:
    rule: (Optional) A rule function or None. Default rule calls build().
    concrete: If True, make this a toplevel model. **Default** - False.
    ctype: (Optional) Pyomo ctype of the Block.
    default: dict with default arguments to ProcessBlockData init
    initialize: dict with block index keys where the values are argument dicts
        for specific ProcessBlockData element init. If a key is missing,
        default will be used.\n"""

class _IndexedProcessBlockMeta(type):
    """Metaclass used to create an indexed model class."""

    def __new__(meta, name, bases, dct):
        def __init__(self, *args, **kwargs):
            kwargs.setdefault("rule", _rule_default)
            kwargs.setdefault("default", {})
            kwargs.setdefault("initialize", {})
            self._block_data_config_default = kwargs.pop("default")
            self._block_data_config_initialize = kwargs.pop("initialize")
            bases[0].__init__(self, *args, **kwargs)
        dct["__init__"] = __init__
        dct["__process_block__"] = "indexed"
        return type.__new__(meta, name, bases, dct)


class _ScalarProcessBlockMeta(type):
    """Metaclass used to create a scalar model class."""

    def __new__(meta, name, bases, dct):
        def __init__(self, *args, **kwargs):
            kwargs.setdefault("rule", _rule_default)
            kwargs.setdefault("default", {})
            kwargs.setdefault("initialize", {})
            self._block_data_config_default = kwargs.pop("default")
            self._block_data_config_initialize = kwargs.pop("initialize")
            bases[0].__init__(self, component=self)
            bases[1].__init__(self, *args, **kwargs)
        dct["__init__"] = __init__
        dct["__process_block__"] = "scalar"
        return type.__new__(meta, name, bases, dct)


class ProcessBlock(Block):
    """Process block.

    Process block behaves like a Pyomo Block. The important differences are
    listed below.

    * There is a default rule that calls the build() method for _BlockData
      subclass ojects, so subclass of _BlockData used in a ProcessBlock
      should have a build() method. A different rule or no rule (None) can be
      set with the usual rule argument, if additional steps are required to
      build an element of a block. A example of such a case is where
      different elements of an indexed block require addtional
      information to construct.

    * Some of the arguments to __init__, which are not expected arguments of
      Block, are split off and stored in self._block_data_config. If the
      _BlockData subclass inherits ProcessBlockData, self._block_data_config
      is sent to the self.config ConfigBlock.

    """
    def __new__(cls, *args, **kwds):
        """Create a new indexed or scalar ProcessBlock subclass instance
        depending on whether there are args.  If there are args those should be
        an indexing set."""
        if hasattr(cls, "__process_block__"):
            # __process_block__ is a class attribute created when making an
            # indexed or scalar subclass of ProcessBlock (or subclass thereof).
            # If cls dosen't have it, the indexed or scalar class has not been
            # created yet.
            #
            # You get here after creating a new indexed or scalar class in the
            # next if below. The first time in, cls is a ProcessBlock subclass
            # that is neither indexed or scalar so you go to the if below and
            # create an index or scalar subclass of cls.
            return super(Block, cls).__new__(cls)
        if args == (): # no args so make scalar class
            bname = "_Scalar{}".format(cls.__name__)
            n = _ScalarProcessBlockMeta(bname, (cls._ComponentDataClass, cls),{})
            return n.__new__(n) #calls this __new__() again with scalar class
        else: # args so make indexed class
            bname = "_Indexed{}".format(cls.__name__)
            n = _IndexedProcessBlockMeta(bname, (cls,), {})
            return n.__new__(n) #calls this __new__() again with indexed class

    @classmethod
    def base_class_name(cls):
        """Name given by the user to the ProcessBase class.

        Return:
           (str) Name of the class.
        Raises:
           AttributeError, if no base class name was set, e.g. this class
               was *not* wrapped by the `declare_process_block_class`
               decorator.

        """
        return cls._orig_name

    @classmethod
    def base_class_module(cls):
        """Return module of the associated ProcessBase class.

        Return:
           (str) Module of the class.
        Raises:
           AttributeError, if no base class module was set, e.g. this class
               was *not* wrapped by the `declare_process_block_class` decorator.

        """
        return cls._orig_module


def declare_process_block_class(name, block_class=ProcessBlock, doc=""):
    """Declare a new ProcessBlock subclass.

    This is a decorator function for a class definition, where the class is
    derived from _BlockData. It creates a ProcessBlock subclass to contain it.
    For example (where ProcessBlockData is a subclass of _BlockData):

    @declare_process_block_class(name=MyUnitBlock)
    class MyUnitBlockData(ProcessBlockData):
        # This class is a _BlockData subclass contained in a Block subclass
        # MyUnitBlock
        ....

    The only requirment is that the subclass of _BlockData contain a build()
    method.

    Args:
        name: class name for the model.
        block_class: ProcessBlock or a subclass of ProcessBlock, this allows
            you to use a subclass of ProcessBlock if needed.
        doc: Documentation for the class. This should play nice with sphinx.

    """
    def proc_dec(cls):  # Decorator function
        # create a new class called name from block_class
        try:
            cb_doc = cls.CONFIG.generate_documentation(
                block_start="", block_end="", item_start="%s: ",
                indent_spacing=4, item_body="%s", item_end="")
            #cb_doc = '\n'.join(' '*4+x for x in cb_doc.splitlines())
        except:
            cb_doc = ""
        ds = "{}\n{}{}\nReturns:\n   New {} instance"\
            .format(doc, _process_block_docstring, cb_doc, name)
        c = type(name, (block_class,),
                {"__module__": cls.__module__,
                 "_ComponentDataClass": cls,
                 "__doc__":ds})
        setattr(sys.modules[cls.__module__], name, c)
        setattr(cls, '_orig_name', name)
        setattr(cls, '_orig_module', cls.__module__)
        return cls
    return proc_dec  # return decorator function
