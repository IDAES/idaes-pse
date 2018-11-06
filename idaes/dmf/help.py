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
Find documentation for modules and classes in the
generated Sphinx documentation and return its location.
"""
import os
from lxml import html
from idaes.dmf import util

_log = util.get_logger('help')


def find_html_docs(dmf, obj, **kw):
    """Get one or more files with HTML documentation for
    the given object, in paths referred to by the dmf instance.
    """
    module_, name = _get_object_name(obj)
    return get_html_docs(dmf, module_, name, **kw)


def get_html_docs(dmf, module_, name, sphinx_version=(1, 5, 5)):
    paths = dmf.get_doc_paths()
    if not paths:
        raise ValueError('No documentation locations configured')

    _log.info('find HTML docs for module={} class={} on paths={}'
              .format(module_, name, paths))
    filenames = []
    for p in paths:
        _log.debug('examine help path "{}"'.format(p))
        html_file = os.path.join(p, 'genindex.html')
        if os.path.exists(html_file):
            _log.debug('get_html_docs: find refs in file={}'.format(html_file))
            refs = _find_refs(html_file, module_, name, sphinx_version)
            if refs:
                if os.path.isabs(p):
                    filenames = [os.path.join(p, r) for r in refs]
                else:
                    filenames = [os.path.join(dmf.root, p, r) for r in refs]
                break
    return filenames


def _get_object_name(obj):
    if hasattr(obj, '_orig_module'):
        module = obj._orig_module
    else:
        module = obj.__module__
    if hasattr(obj, '_orig_name'):
        name = obj._orig_name
    elif isinstance(obj, type):
        name = obj.__name__  # class name of a class
    else:
        name = obj.__class__.__name__
    return module, name


def _find_refs(html_file, module, name, sphinx_version):
    """Find the refs in the generated HTML.

    Right now `sphinx_version` isn't used, but assuming the
    generator changes this may be useful in the future to choose
    a different XPath expression.
    """
    html_content = open(html_file).read()
    tree = html.fromstring(html_content)
    xpath_expr = '//li[contains(.,"{m}")]/ul/li/a'.format(m=module)
    _log.debug('_find_refs: xpath expr={}'.format(xpath_expr))
    # extract all matches to the 'name' from the refs
    elements = tree.xpath(xpath_expr)
    hrefs = [e.get('href') for e in elements if e.text.strip() == name]
    return hrefs
