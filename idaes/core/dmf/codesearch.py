#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Search through the code and index static information in the DMF.
"""
# stdlib
import glob
import importlib
import inspect
import logging
import os
import pprint
import re

__author__ = "Dan Gunter <dkgunter@lbl.gov>"

_log = logging.getLogger(__name__)


class Walker(object):
    def walk(self, visitor):
        """Interface for walkers.

        Args:
            visitor (Visitor): Class whose `visit` method will be
                               called for each item.

        Returns:
            None
        """
        pass


class ModuleClassWalker(Walker):
    """Walk modules from a given root (e.g. 'idaes'), and visit
    all classes in those modules whose name matches a given pattern.

    Example usage::

        walker = ModuleClassWalker(from_pkg=idaes,
                                   class_expr='_PropertyParameter.*')

        walker.walk(PrintMetadataVisitor())  # see below

    """

    def __init__(
        self,
        from_path=None,
        from_pkg=None,
        class_expr=None,
        parent_class=None,
        suppress_warnings=False,
        exclude_testdirs=True,
        exclude_tests=True,
        exclude_init=True,
        exclude_setup=True,
        exclude_dirs=None,
    ):
        """Constructor. Create from either a path or package root, or both.
        If just one is given, the other is deduced. At least one must be
        given.

        Args:
            from_path (str): Path to start looking for modules. Overrides
                            'from_pkg'.
            from_pkg (module): Package from which to start looking for modules.
                            The path of the package will depend on the Python
                            environment and path.
            class_expr (str): Regular expression for the classes to find.
            parent_class (class): If given, filter for this parent class.
                                  Will match class_expr OR parent_class.
            suppress_warnings (bool): If True do not log any warning messages.
            exclude_testdirs (bool): If True (the default), exclude test dirs
            exclude_tests (bool): If True (the default), exclude test files
            exclude_init (bool): If True (the default), exclude __init___.py
            exclude_setup (bool): If True (the default), exclude setup.py
            exclude_dirs (list): List of directory (regex) patterns to exclude.
                                 These will be prefixed with a directory
                                 separator (e.g. '/') but otherwise can be
                                 any path expression.
        """
        if from_path:
            self._root = from_path
            if from_pkg:
                self._pkg = from_pkg.__name__
            else:
                self._pkg = str(os.path.basename(self._root))
        elif from_pkg:
            self._pkg = from_pkg.__name__
            try:
                self._root = from_pkg.__path__[0]
            except AttributeError:
                raise IOError(
                    "Root directory cannot be deduced from"
                    ' package "{}": no __path__ attribute.'.format(self._pkg)
                )
        else:
            raise ValueError(
                'Missing arguments: either "from_pkg" or ' '"from_path" must be given'
            )
        if not os.path.isdir(self._root):
            raise IOError('Root directory "{}"'.format(self._root))
        self._expr = re.compile(class_expr) if class_expr else None
        self._warn = not suppress_warnings
        self._parent = parent_class
        # build regular expression of things to exclude
        expr_list = []
        psep = os.path.sep
        if psep == "\\":
            psep = "\\\\"
        if exclude_testdirs:
            expr_list.append(r"{sl}tests?{sl}".format(sl=psep))
        if exclude_tests:
            expr_list.append(r"tests?_.*.py")
        if exclude_init:
            expr_list.append(r"__init__\.py")
        if exclude_setup:
            expr_list.append(r"setup\.py")
        if exclude_dirs:
            for ed in exclude_dirs:
                expr_list.append("{sl}{d}".format(sl=psep, d=ed))
        self._exclude_expr = re.compile(r"|".join(expr_list))
        self._history = []
        _log.debug("exclude expr={}".format(self._exclude_expr.pattern))

    def walk(self, visitor):
        modules = self._get_modules()
        self._visit_subclasses(modules, visitor.visit)

    def get_indexed_classes(self):
        return self._history

    def _get_modules(self):
        _log.debug("getting modules from root: {}".format(self._root))
        # change file paths at 'root' to module paths from 'pkgroot'
        n, module_list = len(self._root), []
        for f in self._python_files():
            module_path = os.path.splitext(f[n + 1 :])[0]
            module_path = module_path.replace(os.path.sep, ".")
            module_list.append(self._pkg + "." + module_path)
        return module_list

    def _python_files(self):
        q = [self._root]
        while q:
            curdir = q.pop()
            for f in glob.glob(os.path.join(curdir, "*.py")):
                if not self._exclude_expr.search(f):
                    yield f
            for d in os.listdir(curdir):
                path = os.path.join(curdir, d)
                if os.path.isdir(path):
                    q.append(path)

    def _visit_subclasses(self, modules, visit):
        for modname in modules:
            _log.debug("visit module: {}".format(modname))
            try:
                mod = importlib.import_module(modname)
            except Exception:
                if self._warn:
                    _log.warn(
                        "Error during import of module: {}. Ignoring.".format(modname)
                    )
                continue
            for item in dir(mod):
                x = getattr(mod, item)
                if inspect.isclass(x):
                    fullname = modname + "." + x.__name__
                    if not self._expr and not self._parent:
                        if visit(x):
                            self._history.append(fullname)
                    elif self._expr:
                        if self._expr.match(x.__name__):
                            if visit(x):
                                self._history.append(fullname)
                    elif self._parent and issubclass(x, self._parent):
                        if visit(x):
                            self._history.append(fullname)


class Visitor(object):
    """Interface for the 'visitor' class passed to Walker subclasses'
    `walk()` method.
    """

    def visit(self, obj):
        """Visit one object.

        Args:
            obj (object): Some object to operate on.

        Returns:
            True if visit succeeded, else False
        """
        pass


class PropertyMetadataVisitor(Visitor):
    """Visit something implementing :class:`HasPropertyClassMetadata` and
    pass that metadata, as a dict, to the `visit_metadata()` method,
    which should be implemented by the subclass.
    """

    def visit(self, obj):
        """Visit one object.

        Args:
            obj (idaes.core.property_base.HasPropertyClassMetadata): The object

        Returns:
            True if visit succeeded, else False
        """
        result, meta = True, None
        try:
            meta = obj.get_metadata()
        except (AttributeError, TypeError) as err:
            module = obj.__module__
            _log.debug(
                "Cannot get metadata for {}.{}: {}".format(module, obj.__name__, err)
            )
            result = False
        except NotImplementedError:
            module = obj.__module__
            #            if not module.startswith('idaes.core.'):
            _log.warn(
                '{} in module "{}" does not define its metadata'.format(
                    obj.__name__, module
                )
            )
            result = False
        if result:
            self.visit_metadata(obj, meta)
        return result

    def visit_metadata(self, obj, meta):
        """Do something with the metadata.

        Args:
            obj (object): Object from which metadata was pulled, for context.
            meta (idaes.core.property_base.PropertyClassMetadata): The metadata

        Returns:
            None
        """
        pass


class PrintPropertyMetadataVisitor(PropertyMetadataVisitor):
    def visit_metadata(self, obj, meta):
        """Print the module and class of the object, and
        then the metadata dict, to standard output.
        """
        meta_dict = {"properties": meta.properties, "default_units": meta.default_units}
        print("----")
        print("{}: class {}".format(obj.__module__, obj.__name__))
        print("  metadata: {}".format(pprint.pformat(meta_dict)))


# This is for testing


class _TestClass(object):
    pass


class _TestClass1(_TestClass):
    @classmethod
    def get_metadata(cls):
        return {"x": 1}


class _TestClass2(_TestClass):
    @classmethod
    def get_metadata(cls):
        return {"x": 2}


class _TestClass3(_TestClass):
    @classmethod
    def get_metadata(cls):
        return {"x": 3}
