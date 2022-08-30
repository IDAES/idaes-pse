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
Find documentation for modules and classes in the
generated Sphinx documentation and return its location.
"""
# stdlib
import logging
import os
from pathlib import Path
import requests
import types
from urllib.parse import urlparse

# third-party
from lxml import html

# package
from idaes.core.dmf.dmfbase import DMFConfig

__author__ = "Dan Gunter <dkgunter@lbl.gov>"

_log = logging.getLogger(__name__)


def find_html_docs(dmf, obj=None, obj_name=None, **kw):
    """Get one or more files with HTML documentation for
    the given object, in paths referred to by the dmf instance.
    """
    if obj_name is None:
        module_, name = _get_object_name(obj)
    else:
        if obj_name.endswith("."):
            module_, name = obj_name[:-1], ""
        else:
            try:
                module_, name = obj_name.rsplit(".", 1)
            except ValueError:
                raise
    _log.debug("Find docs for object. module=/{}/ name=/{}/".format(module_, name))
    return get_html_docs(dmf, module_, name, **kw)


def get_html_docs(dmf, module_, name, sphinx_version=(1, 5, 5)):
    paths = dmf.get_doc_paths()

    if len(paths) == 0:
        # Error: no configuration, or no paths configured
        conf_path = DMFConfig.configuration_path()
        if DMFConfig.configuration_exists():
            raise ValueError(
                f"No documentation locations configured. "
                f"To set this path, set '{dmf.CONF_HELP_PATH}' in the "
                f"DMF configuration file: {conf_path}"
            )
        else:
            raise ValueError(
                f"No DMF configuration file found. " f"Expected location: {conf_path}"
            )

    _log.info(f"Get HTML docs for module={module_} class={name} on paths={paths}")
    location = None
    for p in paths:
        _log.debug(f"Examine help path '{p}'")
        parsed = urlparse(p)
        html_content = None
        is_web = parsed.scheme in ("http", "https")
        if is_web:
            _log.debug(f"Get help from online documentation")
            url = p + "/genindex.html"
            _log.debug(f"(Help) reading index file: {url}")
            try:
                response = requests.get(url)
                if response.status_code == 200:
                    html_content = response.text
                else:
                    raise requests.RequestException(
                        f"Bad response code: {response.status_code}"
                    )
            except requests.RequestException as err:
                _log.debug(f"Error getting documentation index from '{url}': {err}")
        else:
            _log.debug(f"Get help from local documentation files")
            index_path = Path(p) / "genindex.html"
            if index_path.exists():
                _log.debug(f"(Help) reading index file: {index_path}")
                html_file = str(index_path)
                html_content = open(html_file).read()
            else:
                _log.debug(f"Index 'genindex.html' not found at path '{p}'")
        if html_content:
            _log.debug("parsing index file")
            tree = html.fromstring(html_content)
            # Look for manual references first
            refs = _find_refs(tree, module_, name, sphinx_version)
            # Get full paths for relative refs
            if refs:
                ref = _pick_best_ref(refs)
                if is_web:
                    parsed_ref = urlparse(ref)
                    location = f"{p}/{parsed_ref.path}"
                else:  # file
                    if os.path.isabs(p):
                        location = os.path.join(p, ref)
                    else:
                        location = os.path.join(dmf.root, p, ref)
        if location:
            _log.debug(f"Found help documentation at: {location}")
            break  # stop once we find something
    return [location] if location else []


def _get_object_name(obj):
    if isinstance(obj, types.ModuleType):
        module, name = obj.__name__, ""
    else:
        if hasattr(obj, "_orig_module"):
            module = obj._orig_module
        else:
            module = obj.__module__
        if hasattr(obj, "_orig_name"):
            name = obj._orig_name
        elif isinstance(obj, type):
            name = obj.__name__  # class name of a class
        else:
            name = obj.__class__.__name__
    return module, name


def _find_refs(tree, module, name, sphinx_version):
    """Find object/etc. refs in the HTML."""
    _log.debug(f"looking for module {module}")
    # There are 3 types of refs we are looking for
    # 1. manually added with index directive
    xpath_expr = '//li[contains(.,"{m}")]/ul/li/a'.format(m=module)
    _log.debug("manually indexed: xpath expr={}".format(xpath_expr))
    elements = tree.xpath(xpath_expr)
    hrefs = [e.get("href") for e in elements if e.text.strip() == name]
    if hrefs:
        _log.debug("found manually indexed docs at: {}".format(hrefs[0]))
        return hrefs  # found something; stop
    if name:
        # 2a. embedded in the page
        target = "{}.{}".format(module, name)
        subtype = "embedded"  # for logging
    else:
        # 2b. generated by autodoc in the apidoc folder
        target = "module-{}".format(module)  # how Sphinx names these
        subtype = "autodoc"  # for logging
    xpath_expr = '//li/a[contains(@href,"#{}")]/@href'.format(target)
    _log.debug("{}: xpath expr={}".format(subtype, xpath_expr))
    hrefs = [p for p in tree.xpath(xpath_expr) if p.endswith(target)]
    if hrefs:
        _log.debug("found {} docs at: {}".format(subtype, hrefs[0]))
    return hrefs


def _pick_best_ref(refs):
    # trivial case:q
    if len(refs) == 1:
        return refs[0]
    # find any non-apidoc ref
    for r in refs:
        if not r.startswith("apidoc"):
            return r
    # otherwise just return first
    return refs[0]
