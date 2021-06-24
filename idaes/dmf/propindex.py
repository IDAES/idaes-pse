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
Index Property metadata
"""
# stdlib
import logging

# local
import idaes
from idaes.dmf import codesearch
from idaes.dmf import resource
from idaes.core import property_base  # noqa: F401

__author__ = "Dan Gunter <dkgunter@lbl.gov>"

_log = logging.getLogger(__name__)


def index_property_metadata(
    dmf, pkg=idaes, expr="_PropertyMetadata.*", default_version="0.0.1", **kwargs
):
    """Index all the PropertyMetadata classes in this package.

    Usually the defaults will be correct, but you can modify the package
    explored and set of classes indexed.

    When you re-index the same class (in the same module), whether or
    not that is a "duplicate" will depend on the version found in the
    containing module. If there is no version in the containing module,
    the default version is used (so it is always the same). If it is
    a duplicate, nothing is done, this is not considered an error. If a new
    version is added, it will be explicitly connected to the highest version
    of the same module/code. So, for example,

    1. Starting with (a.module.ClassName version=0.1.2)
    2. If you then find a new version (a.module.ClassName version=1.2.3)
       There will be 2 resources, and you will have the relation::

           a.module.ClassName/1.2.3 --version---> a.module.ClassName/0.1.2

    3. If you add another version (a.module.ClassName version=1.2.4), you
       will have two relations::

           a.module.ClassName/1.2.3 --version---> a.module.ClassName/0.1.2
           a.module.ClassName/1.2.4 --version---> a.module.ClassName/1.2.3

    Args:
        dmf (idaes.dmf.DMF): Data Management Framework instance in which
                             to record the found metadata.
        pkg (module): Root module (i.e. package root) from which to find
                      the classes containing metadata.
        expr (str): Regular expression pattern for the names of the classes
                    in which to look for metadata.
        default_version (str): Default version to use for modules
                    with no explicit version.
        kwargs: Other keyword arguments passed to
                      :class:`codesearch.ModuleClassWalker`.
    Returns:
        codesearch.ModuleClassWalker: Class that walked through the modules.
            You can call `.get_indexed_classes()` to see the list of classes
            walked, or `.walk()` to walk the modules again.
    Raises:
        This instantiated a `DMFVisitor` and calls its `walk()` method to
        walk/visit each found class, so any exception raised by the constructor
        or `DMFVisitor.visit_metadata()`.
    """
    wlk = codesearch.ModuleClassWalker(
        from_pkg=pkg,
        class_expr=expr,
        parent_class=idaes.core.property_meta.HasPropertyClassMetadata,
        suppress_warnings=True,
        **kwargs
    )
    vst = DMFVisitor(dmf, default_version=default_version)
    wlk.walk(vst)
    return wlk


class DMFVisitor(codesearch.PropertyMetadataVisitor):

    #: Added to resource 'tags', so easier to find later
    INDEXED_PROPERTY_TAG = "indexed-property"

    def __init__(self, dmf, default_version=None):
        """Constructor.

        Args:
            dmf (idaes.dmf.DMF): Data management framework.
            default_version (Union[None, list, tuple]): Default version to give
                 the class, if the containing module does not have
                 a `__version__` variable. If None, the absence of that
                 variable will cause an error.
        Raises:
            TypeError: if `default_version` isn't something that
                       :func:`resource.version_list` can convert.
        """
        self._dmf = dmf
        self._defver = None
        if default_version is not None:
            try:
                self._defver = resource.version_list(default_version)
            except ValueError as err:
                raise TypeError('Bad "default_version": {}'.format(err))

    def visit_metadata(self, obj, meta):
        """Called for each property class encountered during the "walk"
         initiated by `index_property_metadata()`.

        Args:
            obj (property_base.PropertyParameterBase): Property class instance
            meta (property_base.PropertyClassMetadata): Associated metadata
        Returns:
            None
        Raises:
            AttributeError: if
        """
        _log.debug(
            "Adding resource to DMF that indexes the property package "
            '"{}"'.format(".".join([obj.__module__, obj.__name__]))
        )
        r = resource.Resource(type_=resource.ResourceTypes.code)
        r.data = {"units": meta.default_units, "properties": meta.properties}
        containing_module = obj.__module__
        if hasattr(containing_module, "__version__"):
            obj_ver = resource.version_list(containing_module.__version__)
        elif self._defver is None:
            raise AttributeError(
                "No __version__ for module {}, and no "
                "default".format(containing_module)
            )
        else:
            obj_ver = self._defver
        r.v["codes"].append(
            {
                "type": "class",
                "language": "python",
                "name": ".".join([obj.__module__, obj.__name__]),
                "version": obj_ver,
            }
        )
        r.v["tags"].append(self.INDEXED_PROPERTY_TAG)
        # Search for existing indexed codes.
        # A match exists if all 3 of these are the same:
        #   codes.type == class
        #   codes.language == python
        #   codes.name == <module>.<class>
        info = {k: r.v["codes"][0][k] for k in ("type", "language", "name")}
        rsrc_list, dup_rsrc = [], None
        # Loop through all the right kind of resources
        for rsrc in self._dmf.find(
            {r.TYPE_FIELD: resource.ResourceTypes.code, "tags": ["indexed-property"]}
        ):
            # skip any resources without one code
            if len(rsrc.v["codes"]) != 1:
                continue
            code = rsrc.v["codes"][0]
            # skip any resource of wrong code type, name, lang.
            skip = False
            for k in info:
                if code[k] != info[k]:
                    skip = True
                    break
            if skip:
                continue
            # skip any resources missing the recorded metadata
            skip = False
            for data_key in r.data.keys():
                if data_key not in rsrc.data:
                    skip = True
                    break
            if skip:
                continue
            # If the version of the found code is the same as the
            # version of the one to be added, then it is a duplicate
            if code["version"] == obj_ver:
                dup_rsrc = rsrc
                break
            rsrc_list.append(rsrc)
        if dup_rsrc:
            # This is considered a normal, non-exceptional situation
            _log.debug(
                "DMFVisitor: Not adding duplicate index for "
                "{}v{}".format(info["name"], obj_ver)
            )
        else:
            # add the resource
            r.validate()
            _log.debug(
                'DMFVisitor: Adding resource for code "{}"v{} type={}'.format(
                    r.v["codes"][0]["name"],
                    r.v["codes"][0]["version"],
                    r.v["codes"][0]["type"],
                )
            )

            self._dmf.add(r)
            if rsrc_list:
                # Connect to most recent (highest) version
                rsrc_list.sort(key=lambda rs: rs.v["codes"][0]["version"])
                # for rsrc in rsrc_list:
                rsrc = rsrc_list[-1]
                rel = resource.Triple(r, resource.Predicates.version, rsrc)
                resource.create_relation(rel)
                self._dmf.update(rsrc)
                self._dmf.update(r)
