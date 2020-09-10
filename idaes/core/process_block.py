##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
The process_block module simplifies inheritance of Pyomo blocks. The main
reason to subclass a Pyomo block is to create a block that comes with
pre-defined model equations. This is used in the IDAES modeling framework to
create modular process model blocks.
"""

import sys
import logging

from pyomo.common.config import ConfigBlock
from pyomo.environ import Block

__author__ = "John Eslick"
__all__ = ['ProcessBlock', 'declare_process_block_class']


def _rule_default(b, *args):
    """
    Default rule for ProcessBlock, which calls build(). A different rule can
    be specified to add additional build steps, or to not call build at all
    using the normal rule argument to ProcessBlock init.
    """
    try:
        b.build()
    except Exception:
        logging.getLogger(__name__).exception(
            "Failure in build: {}".format(b))
        raise

_process_block_docstring = """
    Args:
        rule (function): A rule function or None. Default rule calls build().
        concrete (bool): If True, make this a toplevel model. **Default** - False.
        ctype (str): Pyomo ctype of the block.  **Default** - "Block"
        default (dict): Default ProcessBlockData config{}
        initialize (dict): ProcessBlockData config for individual elements. Keys
            are BlockData indexes and values are dictionaries described under the
            "default" argument above.
        idx_map (function): Function to take the index of a BlockData element and
            return the index in the initialize dict from which to read arguments.
            This can be provided to overide the default behavior of matching the
            BlockData index exactly to the index in initialize.
    Returns:
        ({}) New instance
    """

_config_block_keys_docstring = """

            ..

            Keys
{}
            ..
"""

def _process_kwargs(o, kwargs):
    kwargs.setdefault("rule", _rule_default)
    o._block_data_config_default = kwargs.pop("default", None)
    o._block_data_config_initialize = ConfigBlock(implicit=True)
    o._block_data_config_initialize.set_value(kwargs.pop("initialize", None))
    o._idx_map = kwargs.pop("idx_map", None)


class _IndexedProcessBlockMeta(type):
    """Metaclass used to create an indexed model class."""

    def __new__(meta, name, bases, dct):
        def __init__(self, *args, **kwargs):
            _process_kwargs(self, kwargs)
            bases[0].__init__(self, *args, **kwargs)
        dct["__init__"] = __init__
        dct["__process_block__"] = "indexed"
        return type.__new__(meta, name, bases, dct)


class _ScalarProcessBlockMeta(type):
    """Metaclass used to create a scalar model class."""

    def __new__(meta, name, bases, dct):
        def __init__(self, *args, **kwargs):
            _process_kwargs(self, kwargs)
            bases[0].__init__(self, component=self)
            bases[1].__init__(self, *args, **kwargs)
        dct["__init__"] = __init__
        dct["__process_block__"] = "scalar"
        return type.__new__(meta, name, bases, dct)


class ProcessBlock(Block):
    __doc__ = """
    ProcessBlock is a Pyomo Block that is part of a system to make Pyomo
    Block easier to subclass. The main difference between a Pyomo Block and
    ProcessBlock from the user perspective is that a ProcessBlock has a rule
    assigned by default that calls the build() method for the contained
    ProcessBlockData objects. The default rule can be overridden, but the new
    rule should always call build() for the ProcessBlockData object.
    """ + _process_block_docstring.format("", "ProcessBlock")

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
    """
    Declare a new ProcessBlock subclass.

    This is a decorator function for a class definition, where the class is
    derived from Pyomo's _BlockData. It creates a ProcessBlock subclass to
    contain the decorated class. The only requirment is that the subclass of
    _BlockData contain a build() method. The purpose of this decorator is to
    simplify subclassing Pyomo's block class.

    Args:
        name: name of class to create
        block_class: ProcessBlock or a subclass of ProcessBlock, this allows
            you to use a subclass of ProcessBlock if needed. The typical use
            case for Subclassing ProcessBlock is to impliment methods that
            operate on elements of an indexed block.
        doc: Documentation for the class. This should play nice with sphinx.

    Returns:
        Decorator function

    """
    def proc_dec(cls):  # Decorator function
        # create a new class called name from block_class
        try:
            cb_doc = cls.CONFIG.generate_documentation(
                block_start="", block_end="", item_start="%s\n",
                indent_spacing=4, item_body="%s", item_end="\n", width=66)
            cb_doc += "\n"
            cb_doc = '\n'.join(' '*12 + x for x in cb_doc.splitlines())
        except:
            cb_doc = ""
        if cb_doc != "":
            cb_doc = _config_block_keys_docstring.format(cb_doc)
        ds = "\n".join([doc, _process_block_docstring.format(cb_doc, name)])
        c = type(name, (block_class,),
                {"__module__": cls.__module__,
                 "_ComponentDataClass": cls,
                 "__doc__":ds})
        setattr(sys.modules[cls.__module__], name, c)
        setattr(cls, '_orig_name', name)
        setattr(cls, '_orig_module', cls.__module__)
        return cls
    return proc_dec  # return decorator function
