# coding: utf-8
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

import importlib
import inspect
import sys

import idaes.logger as idaeslog

_declared_beta_modules = set()
_declared_beta_module_imports = set()
_imported_beta_modules = {}


def _caller_module_name():
    # Move up the stack twice: we want the module name of the
    # function/module that called the function calling
    # _caller_module_name
    caller_module = inspect.getmodule(inspect.currentframe().f_back.f_back)
    if caller_module is None:
        return "<stdin>"
    return caller_module.__name__


def declare_beta_module(message=None):
    mod_name = _caller_module_name()
    _declared_beta_modules.add(mod_name)
    if mod_name in _declared_beta_module_imports:
        return

    if message is None:
        message = (
            "Module '%s' is in beta and must be imported using "
            "idaes.beta.import_beta()." % (mod_name,)
        )
    raise ImportError(message)


def import_beta(name, package=None):
    if package is None and name.startswith("."):
        package_info = _caller_module_name().rsplit(".", 1)
        if len(package_info) == 2:
            package = package_info[0]
    absolute_name = importlib.util.resolve_name(name, package)

    if absolute_name in _imported_beta_modules:
        # This is a beta module that we have already imported.  Do not
        # re-import it and instead return the cached imported module
        return _imported_beta_modules[absolute_name]
    else:
        _declared_beta_module_imports.add(absolute_name)
        try:
            _imported_beta_modules[absolute_name] = module = importlib.import_module(
                absolute_name
            )

            if absolute_name in _declared_beta_modules:
                # Remove this (BETA) module from sys.modules so that
                # subsequent imports will re-trigger the
                # declare_beta_module test in the module.
                print("removing %s" % (absolute_name,))
                del sys.modules[absolute_name]
            else:
                idaeslog.getLogger(absolute_name).info(
                    "Module '%s' imported module '%s' as a Beta module.  "
                    "This module is not declared beta and can be imported "
                    "using Python's normal import mechanisms."
                    % (
                        _caller_module_name(),
                        absolute_name,
                    )
                )
                # Remove this (standard) module from the
                # _imported_beta_modules cache
                _imported_beta_modules.pop(absolute_name)
            return module
        finally:
            _declared_beta_module_imports.remove(absolute_name)
