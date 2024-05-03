#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
The process_block module simplifies inheritance of Pyomo blocks. The main
reason to subclass a Pyomo block is to create a block that comes with
pre-defined model equations. This is used in the IDAES modeling framework to
create modular process model blocks.
"""
# TODO: Look into if this is necessary
# pylint: disable=protected-access

import sys
import logging
import inspect

from pyomo.common.config import ConfigBlock, String_ConfigFormatter
from pyomo.common.pyomo_typing import get_overloads_for
from pyomo.core.base.global_set import UnindexedComponent_set
from pyomo.environ import Block

__author__ = "John Eslick"
__all__ = ["ProcessBlock", "declare_process_block_class"]


def _rule_default(b, *args):
    """
    Default rule for ProcessBlock, which calls build(). A different rule can
    be specified to add additional build steps, or to not call build at all
    using the normal rule argument to ProcessBlock init.
    """
    try:
        b.build()
    except Exception:
        logging.getLogger(__name__).exception(f"Failure in build: {b}")
        raise


_process_block_docstring = """
    Args:
        rule (function): A rule function or None. Default rule calls build().
        concrete (bool): If True, make this a toplevel model. **Default** - False.
        ctype (class): Pyomo ctype of the block.  **Default** - pyomo.environ.Block
        {}
        initialize (dict): ProcessBlockData config for individual elements. Keys
            are BlockData indexes and values are dictionaries with config arguments 
            as keys.
        idx_map (function): Function to take the index of a BlockData element and
            return the index in the initialize dict from which to read arguments.
            This can be provided to override the default behavior of matching the
            BlockData index exactly to the index in initialize.
    Returns:
        ({}) New instance
    """

_config_block_keys_docstring = """

            ..

            Config args
{}
            ..
"""


def _get_pyomo_block_kwargs():
    """This function gets the keyword argument names used by Pyomo Block.__init__
    This list is generated when importing the module rather than a static list
    to accommodate future Pyomo interface changes.
    """
    funcs = get_overloads_for(Block.__init__)
    keywords = set()
    for func in funcs:
        keywords.update(inspect.getfullargspec(func).kwonlyargs)
    return keywords


# Get a list of init kwarg names reserved for the base Pyomo Block class
_pyomo_block_keywords = _get_pyomo_block_kwargs()


class _ScalarProcessBlockMixin(object):
    """Mixin class for scalar process block classes."""

    def __init__(self, *args, **kwargs):
        # __bases__ for the ScalarProcessBlock is
        #
        #    (_ScalarProcessBlockMixin, {process_block_data}, {process_block})
        #
        # Unfortunately, we cannot guarantee that this is being called
        # from the ScalarProcessBlock (someone could have inherited from
        # that class to make another scalar class).  We will walk up the
        # MRO to find the Scalar class (which should be the only class
        # that has this Mixin as the first base class)
        for cls in self.__class__.__mro__:
            if cls.__bases__[0] is _ScalarProcessBlockMixin:
                _mixin, _data, _block = cls.__bases__
                _data.__init__(self, component=self)
                _block.__init__(self, *args, **kwargs)
                break


class ProcessBlock(Block):
    __doc__ = """
    ProcessBlock is a Pyomo Block that is part of a system to make Pyomo
    Block easier to subclass. The main difference between a Pyomo Block and
    ProcessBlock from the user perspective is that a ProcessBlock has a rule
    assigned by default that calls the build() method for the contained
    ProcessBlockData objects. The default rule can be overridden, but the new
    rule should always call build() for the ProcessBlockData object.
    """ + _process_block_docstring.format(
        "", "ProcessBlock"
    )

    def __new__(cls, *args, **kwargs):
        """Create a new indexed or scalar Process Block subclass instance
        depending on whether there are args.  If there are args those should be
        an indexing set."""
        if hasattr(cls, "__process_block__"):
            # __process_block__ is a class attribute created when making
            # an indexed or scalar subclass of ProcessBlock (or subclass
            # thereof).  If it is present, then we can assume the
            # routing of the generic "ProcessBlock" to the specific
            # Scalar or Indexed class has already occurred and we can
            # pass constrol up to (toward) object.__new__()
            return super().__new__(cls, *args, **kwargs)
        # If cls doesn't have __process_block__, the user is attempting
        # to create the "generic" (derived) ProcessBlock.  Depending on
        # the arguments, we need to map this to either the
        # ScalarProcessBlock aor IndexedProcessBlock subclass.
        if not args or (args[0] is UnindexedComponent_set and len(args) == 1):
            return super().__new__(cls._scalar_process_block, *args, **kwargs)
        else:
            return super().__new__(cls._indexed_process_block, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        Block.__init__(self, *args, **self._process_kwargs(kwargs))

    def _process_kwargs(self, kwargs):
        "Filter to separate IDAES __init__ arguments from core Pyomo arguments"

        kwargs.setdefault("rule", _rule_default)
        self._block_data_config_initialize = ConfigBlock(implicit=True)
        self._block_data_config_initialize.set_value(kwargs.pop("initialize", None))
        self._idx_map = kwargs.pop("idx_map", None)

        _pyomo_kwargs = {}
        for arg in _pyomo_block_keywords:
            if arg in kwargs:
                _pyomo_kwargs[arg] = kwargs.pop(arg)

        self._block_data_config_default = kwargs
        return _pyomo_kwargs


def declare_process_block_class(name, block_class=ProcessBlock, doc=""):
    """
    Declare a new ProcessBlock subclass.

    This is a decorator function for a class definition, where the class is
    derived from Pyomo's _BlockData. It creates a ProcessBlock subclass to
    contain the decorated class. The only requirement is that the subclass of
    _BlockData contain a build() method. The purpose of this decorator is to
    simplify subclassing Pyomo's block class.

    Args:
        name: name of class to create
        block_class: ProcessBlock or a subclass of ProcessBlock, this allows
            you to use a subclass of ProcessBlock if needed. The typical use
            case for Subclassing ProcessBlock is to implement methods that
            operate on elements of an indexed block.
        doc: Documentation for the class. This should play nice with sphinx.

    Returns:
        Decorator function

    """

    def proc_dec(block_data):  # Decorator function
        # prepare the main docstring for the new class.
        try:
            cb_doc = block_data.CONFIG.generate_documentation(
                format=String_ConfigFormatter(
                    block_start="%s\n",
                    block_end="",
                    item_start="%s\n",
                    item_body="%s",
                    item_end="\n",
                ),
                indent_spacing=4,
                width=66,
            )
            cb_doc += "\n"
            cb_doc = "\n".join(" " * 12 + x for x in cb_doc.splitlines())
        except Exception:  # pylint: disable=W0703
            cb_doc = ""
        if cb_doc != "":
            cb_doc = _config_block_keys_docstring.format(cb_doc)
        ds = "\n".join([doc, _process_block_docstring.format(cb_doc, name)])

        # Declare the new Block component (derived from CustomBlock)
        # corresponding to the BlockData that we are decorating
        #
        # Note use of `type(block_class)` to get the metaclass that was
        # used to create the base process block class
        comp = type(block_class)(
            name,  # name of new class
            (block_class,),  # base classes
            # class body definitions (populate the new class' __dict__)
            {
                # ensure the created class is associated with the calling module
                "__module__": block_data.__module__,
                # Default IndexedComponent data object is the decorated class:
                "_ComponentDataClass": block_data,
                # Set the docstring
                "__doc__": ds,
            },
        )

        # Declare Indexed and Scalar versions of the process block.  We
        # will register them both with the calling module scope, and
        # with the new ProcessBlock (so that ProcessBlock.__new__ can route
        # the object creation to the correct class)
        comp._indexed_process_block = type(comp)(
            "Indexed" + name,
            (comp,),
            {
                # ensure the created class is associated with the calling module
                "__module__": block_data.__module__,
                # flag for detecting indexed process blocks
                "__process_block__": "indexed",
                # provide function ``base_class_module()`` to get unit
                # module, for visualizer
                "base_class_module": lambda self: block_data.__module__,
            },
        )
        comp._scalar_process_block = type(comp)(
            "Scalar" + name,
            (_ScalarProcessBlockMixin, block_data, comp),
            {
                # ensure the created class is associated with the calling module
                "__module__": block_data.__module__,
                # flag for detecting scalar process blocks
                "__process_block__": "scalar",
                # provide function ``base_class_module()`` to get unit
                # module, for visualizer
                "base_class_module": lambda self: block_data.__module__,
            },
        )

        # Register the new Block types in the same module as the BlockData
        for _cls in (comp, comp._indexed_process_block, comp._scalar_process_block):
            setattr(sys.modules[block_data.__module__], _cls.__name__, _cls)
        return block_data

    return proc_dec  # return decorator function
