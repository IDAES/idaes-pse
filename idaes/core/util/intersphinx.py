#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Provide an `intersphinx_mapping` for Intersphinx_ that can be imported and used 
in the `conf.py` for Sphinx documentation in any IDAES-dependent project.

Usage (in conf.py)::

    from idaes.util.intersphinx import get_intersphinx_mapping, modify_url
    # ...
    # get mapping
    intersphinx_mapping = get_intersphinx_mapping()
    # set version a/o language for individual URLs (or all, with None for key)
    modify_url(intersphinx_mapping, "idaes", version="latest") # etc.

.. _Intersphinx: https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html

"""
from typing import Dict, Optional

__author__ = "Dan Gunter (LBNL)"


# PSE project documentation URLs
_base_mapping = {
    "idaes": ("https://idaes-pse.readthedocs.io/en/stable/", None),
    "idaes-examples": ("https://idaes-examples.readthedocs.io/en/stable/", None),
    "watertap": ("https://watertap.readthedocs.io/en/stable/", None),
    "prommis": ("https://prommis.readthedocs.io/en/stable/", None),
}


def get_intersphinx_mapping() -> Dict:
    """Get the dictionary expected by the Intersphinx extension.

    Note that the mapping returned is whatever the current version of Intersphinx
    used in the IDAES projects expects, and would change if this changed.

    Returns:
        Mapping for Intersphinx
    """
    return _base_mapping.copy()


def modify_url(
    mapping: Dict,
    key: Optional[str] = None,
    version: Optional[str] = None,
    language: Optional[str] = None,
) -> None:
    """Modify the URL of an entry in the mapping.

    URLs are of the form: `http[s]://{site}/{language}/{version}/`

    Args:
        mapping: The mapping to modify
        key: The name of the entry to modify, e.g. 'idaes', or None for all
        version: New value for version of the docs, e.g. 'stable' or 'latest'
        language: New value for language, e.g. "en" for English

    Raises:
        KeyError: If entry 'key' is not found

    Returns:
        None (modified in place)
    """

    def _modify(k):
        entry = mapping[k]
        url = entry[0]
        new_url = _modify_url(url, version=version, language=language)
        mapping[k] = (new_url, entry[1])

    if key is None:
        for key in mapping:
            _modify(key)
    else:
        _modify(key)


def _modify_url(
    url: str, version: Optional[str] = None, language: Optional[str] = None
) -> str:
    """Modify a URL in a mapping."""
    parts = url.split("/")
    trailing_slash_offset = -1 if parts[-1] == "" else 0
    if version is not None:
        parts[trailing_slash_offset - 1] = version
    if language is not None:
        parts[trailing_slash_offset - 2] = language
    modified_url = "/".join(parts)
    return modified_url
