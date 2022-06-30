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
The 'experiment' is a root container for a coherent
set of 'resources'.
"""
# stdlib
from copy import deepcopy
import logging

# local
from idaes.core.dmf import resource, errors
from idaes.core.dmf.resource import Predicates, ResourceTypes

__author__ = "Dan Gunter <dkgunter@lbl.gov>"

_log = logging.getLogger(__name__)


class Experiment(resource.Resource):
    """An experiment is a way of grouping resources in a way that
    makes sense to the user.

    It is also a useful unit for passing as an argument to functions,
    since it has a standard 'slot' for the DMF instance that created it.
    """

    def __init__(self, dmf, **kwargs):
        """Constructor. Adds the new experiment to the DMF.

        Args:
            dmf (DMF): Data Management Framework instance.
            kwargs: Keyword arguments passed to parent class.
        """
        super(Experiment, self).__init__(value=kwargs)
        self.v[self.TYPE_FIELD] = ResourceTypes.experiment
        dmf.add(self)
        self._dmf = dmf

    @property
    def dmf(self):
        return self._dmf

    def add(self, rsrc):
        """Add a resource to an experiment.

        This does two things:

        1. Establishes an "experiment" type of relationship between the
           new resource and the experiment.
        2. Adds the resource to the DMF

        Args:
            rsrc (resource.Resource): The resource to add.
        Returns:
            resource.Resource: Added (input) resource, for chaining calls.
        """
        resource.create_relation(self, Predicates.contains, rsrc)
        self._dmf.update(rsrc, upsert=True)
        self._dmf.update(self)

    def copy(self, new_id=True, **kwargs):
        """Get a copy of this experiment. The returned object will
        have been added to the DMF.

        Args:
            new_id (bool): If True, generate a new unique ID for the copy.
            kwargs: Values to set in new instance after copying.
        Returns:
            Experiment: A (mostly deep) copy.

            Note that the DMF instance is just a reference to the
            same object as in the original, and they will share state.
        """
        new_exp = Experiment(self._dmf)
        new_exp.v = deepcopy(self.v)
        new_exp.v.update(kwargs)
        if new_id:
            new_exp.set_id()
        self._dmf.add(new_exp)
        return new_exp

    def update(self):
        """Update experiment to current values."""
        self._dmf.update(self, sync_relations=True)

    def remove(self):
        """Remove this experiment from the associated DMF instance."""
        # remove from the DMF
        self._dmf.remove(self.id)
        # cut the connection to the DMF instance
        self._dmf = None
        # disable known methods (via monkeypatching!)
        for m in "add", "update", "remove", "link", "copy":
            self.__dict__[m] = self._removed

    def link(self, subj, predicate=Predicates.contains, obj=None):
        """Add and update relation triple in DMF.

        Args:
            subj (resource.Resource): Subject
            predicate (str): Predicate
            obj (resource.Resource): Object

        Returns:
            None
        """
        if obj is None:
            obj = self
        resource.create_relation(subj, predicate, obj)
        self._dmf.update(subj)
        self._dmf.update(obj)

    def _removed(self, *args, **kwargs):
        raise errors.BadResourceError("This experiment has been removed")
