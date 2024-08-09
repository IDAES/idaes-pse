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

    from idaes.util.intersphinx import get_intersphinx_mapping
    # ...
    intersphinx_mapping = get_intersphinx_mapping("latest", language="en")
    intersphinx_mapping.update({}) # add local mappings here

.. _Intersphinx: https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html

"""
from typing import Dict

__author__ = "Dan Gunter (LBNL)"


# IDAES-PSE project documentation URLs
_base_mapping = {
    "idaes": ("https://idaes-pse.readthedocs.io/", None),
    "idaes-examples": ("https://idaes-examples.readthedocs.io/", None),
    "watertap": ("https://watertap.readthedocs.io/", None),
    "prommis": ("https://prommis.readthedocs.io/", None),
}


def get_intersphinx_mapping(version: str = "stable", language: str = "en") -> Dict:
    """Get the dictionary expected by the Intersphinx extension, with values
    pre-processed to produce URLs corresponding to the desired documentation
    version and language. Note that the 'objects.inv' file that Intersphinx is
    looking for is only available at the full, not base, URL.

    Note that the mapping returned is whatever the current version of Intersphinx
    used in the IDAES projects expects, and would change if this changed.

    Args:
        version: Version of the docs, e.g. 'stable' or 'latest'
        language: Language, e.g. "en" for English

    Returns:
        Mapping for Intersphinx
    """
    return {
        key: (_preprocess_url(v[0], version, language), v[1])
        for key, v in _base_mapping.items()
    }


def _preprocess_url(base_url: str, version: str, language: str) -> str:
    """Modify URL with the provided version and language."""
    return f"{base_url}/{language}/{version}/"
