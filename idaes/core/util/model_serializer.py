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
Functions for saving and loading Pyomo objects to json
"""

from pyomo.environ import *
from pyomo.dae import *
from pyomo.core.base.component import ComponentData
import json
import datetime
import time
import gzip
import logging

_log = logging.getLogger(__name__)

# Some more inforation about this module
__author__ = "John Eslick"
__format_version__ = 4


def _can_serialize(o):
    try:
        json.dumps(o)
        return True
    except TypeError:
        return False


def _set_active(o, d):
    """
    Set if component is active, used for read active attribute callback.
    Args:
        o: object whoes attribute is to be set
        d: attribute value
    Returns:
        None
    """
    try:
        if d:
            o.activate()
        else:
            o.deactivate()
    except AttributeError:
        # if the object doesn't have activate and deactivate that's fine just
        # keep going
        pass


def _set_fixed(o, d):
    """
    Set if variable is fixed, used for read fixed attribute callback.
    Args:
        o: object whose attribute is to be set
        d: attribute value
    Returns:
        None
    """
    if d:
        o.fix()
    else:
        o.unfix()


def _get_value(o):
    """
    Get object value attribute callback.
    Args:
        o: object whoes attribute is to be set
        d: attribute value
    Returns:
        None
    """
    return value(o, exception=False)


def _set_value(o, d):
    """
    Set object value attribute callback. This doen't allow the value of an
    immutable paramter to be set (which would raise an exeption in Pyomo)
    Args:
        o: object whoes attribute is to be set
        d: attribute value
    Returns:
        None
    """
    if isinstance(o, Param) and not o._mutable:
        return  # ignore requests to set immutable params
    else:
        try:
            o.value = d
        except AttributeError:
            o = d  # this would be an indexed parameter


def _set_lb(o, d):
    """
    Set variable lower bound, used for read lb attribute callback.
    Args:
        o: object whoes attribute is to be set
        d: attribute value
    Returns:
        None
    """
    o.setlb(d)


def _set_ub(o, d):
    """
    Set variable upper bound, use for read ub attribute callback.
    Args:
        o: object whoes attribute is to be set
        d: attribute value
    Returns:
        None
    """
    o.setub(d)


def _value_if_not_fixed(o, d):
    """
    Returns a list of attributes to read for a variable, only whether it is
    fixed for non-fixed variables and if it is fixed and the value for fixed
    variables.  The allows you to set up a serializer that only reads fixed
    variable values.

    Args:
        o: Pyomo component being loaded
        d: State dictionary for the component o.
    Returns:
        An attribute list to read. Loads fixed for either fixed or un-fixed
        variables, but only reads in values for unfixed variables.  This is
        useful for intialization functions.
    """
    if o.fixed:
        return ()
    return ("value",)


def _only_fixed(o, d):
    """
    Returns a list of attributes to read for a variable, only whether it is
    fixed for non-fixed variables and if it is fixed and the value for fixed
    variables.  The allows you to set up a serializer that only reads fixed
    variable values.

    Args:
        o: Pyomo component being loaded
        d: State dictionary for the component o.
    Returns:
        An attribute list to read. Loads fixed for either fixed or un-fixed
        variables, but only reads in values for unfixed variables.  This is
        useful for intialization functions.
    """
    if d["fixed"]:
        return ("value", "fixed")
    else:
        return ("fixed",)


class Counter(object):
    """
    This is a counter object, which is an easy way to pass an interger pointer
    around between methods.
    """

    def __init__(self):
        self.count = 0


class StoreSpec(object):
    """
    A StoreSpec object tells the serializer functions what to read or write.
    The default settings will produce a StoreSpec configured to load/save the
    typical attributes required to load/save a model state.

    Args:
        classes: Dictionary of classes to read/write. The class is the dictionary
            key and the value is a list or tuple with two elements 1) a list or
            tuple of attributes to read/write and 2) a filter function. The
            attributes can be an empty list and the filter function can be None.
            Filter functions if present only apply to reading components and do
            not affect writing. Filter functions takes two arguments. The first
            is the object being read to and the second is the state dictionary
            being read from. The state dictionary contains keys for the
            attributes that were written. Based on the state of the component
            and the stored state of the component, the filter function returns
            a list of attributes to read.  For example, a filter function can be
            used to only read the values of varaibles that were fixed when they
            were written. The order of the class keys is important and should
            go from specific classes to general ones. For example, if a
            UnitModel is a Block, but you want to store extra attributes for a
            UnitModel, UnitModel should come before Block in the key set.
        data_classes: This takes the same form as the classes argument, the
            classes should be compoent data types.
        ignore_missing: If True will ignore a component or attribute that exists
            in the model, but not in the stored state. If false an excpetion
            will be raised for things in the model that should be loaded but
            aren't in the stored state. Extra items in the stored state will not
            raise an exception regaurdless of this argument.
        suffix: If True store suffixes and component ids.  If false, don't store
            suffixes. This is a legacy option. The prefered way to store a Suffix
            of not is just to include Suffix in classes or not.
        suffix_filter: None to store all suffixes or a list of suffixes to store.
    """

    def __init__(
        self,
        classes={
            Param: (("_mutable",), None),
            Var: ((), None),
            BooleanVar: ((), None),
            Expression: ((), None),
            Block: (("active",), None),
            Constraint: (("active",), None),
            Suffix: ((), None),
        },
        data_classes={
            Var._ComponentDataClass: (("fixed", "stale", "value", "lb", "ub"), None),
            BooleanVar._ComponentDataClass: (("fixed", "stale", "value"), None),
            pyomo.core.base.param._ParamData: (("value",), None),
            int: (("value",), None),
            float: (("value",), None),
            str: (("value",), None),
            Expression._ComponentDataClass: ((), None),
            Block._ComponentDataClass: (("active",), None),
            Constraint._ComponentDataClass: (("active",), None),
        },
        ignore_missing=True,
        suffix=None,
        suffix_filter=None,
    ):
        # convert old style list/tuple classes arg to dict if needed
        if isinstance(classes, (list, tuple)):
            self.classes = {}
            for c in classes:
                self.classes[c[0]] = c[1:]
        else:
            self.classes = classes
        # convert old style list/tuple data_classes arg to dict if needed
        if isinstance(data_classes, (list, tuple)):
            self.data_classes = {}
            for c in data_classes:
                self.data_classes[c[0]] = c[1:]
        else:
            self.data_classes = data_classes
        # Callbacks are used for attributes that cannont be directly get or set
        self.write_cbs = {  # Write callbacks (writing state so get attr)
            "value": _get_value
        }
        self.read_cbs = {  # Reads callbacks (reading in state so set attr)
            "_mutable": None,
            "active": _set_active,
            "fixed": _set_fixed,
            "lb": _set_lb,
            "ub": _set_ub,
            "value": _set_value,
        }
        # Block and BlockData are required for model structure
        if Block not in self.classes:
            self.classes[Block] = ((), None)
        if Block._ComponentDataClass not in data_classes:
            self.data_classes[Block._ComponentDataClass] = ((), None)
        # If suffix is None, deside by whether in classes, else add or remove
        # suffix based on option.  May deprecate the suffix option.
        if suffix is not None:
            if not suffix and Suffix in classes:
                del self.classes[Suffix]
            elif suffix and Suffix not in classes:
                self.classes[Suffix] = ((), None)
        # Create filter function lists, use None if not supplied
        for i, c in self.classes.items():
            if len(c) < 2:
                self.classes[i] = (c[0], None)
        for i, c in self.data_classes.items():
            if len(c) < 2:
                self.data_classes[i] = (c[0], None)
        self.ignore_missing = ignore_missing
        self.suffix_filter = suffix_filter

    def set_read_callback(self, attr, cb=None):
        """
        Set a callback to set an attribute, when reading from json or dict.
        """
        self.read_cbs[attr] = cb

    def set_write_callback(self, attr, cb=None):
        """
        Set a callback to get an attribute, when writing to json or dict.
        """
        self.write_cbs[attr] = cb

    def get_class_attr_list(self, o):
        """
        Look up what attributes to save/load for an Component object.
        Args:
            o: Object to look up attribute list for.
        Return:
            A list of attributes and a filter function for object type
        """
        alist = None  # Attributes to store
        ff = None  # Load filter function
        for cl in self.classes:
            if isinstance(o, cl):
                alist = self.classes[cl][0]
                ff = self.classes[cl][1]
                break
        if isinstance(o, Block._ComponentDataClass):
            # If you're here you are trying to serialize an element of an
            # indexed block at the top level.  We do want to allow that, so
            # we'll pretend it's a block.
            alist = self.classes[Block][0]
            ff = self.classes[Block][1]
        return (alist, ff)

    def get_data_class_attr_list(self, o):
        """
        Look up what attributes to save/load for an ComponentData object.
        Args:
            o: Object to look up attribute list for.
        Return:
            A list of attributes and a filter function for object type
        """
        alist = None  # Attributes to store
        ff = None  # Load filter function
        for cl in self.data_classes:
            if isinstance(o, cl):
                alist = self.data_classes[cl][0]
                ff = self.data_classes[cl][1]
                break
        return (alist, ff)

    @classmethod
    def bound(cls):
        """Returns a StoreSpec object to store variable bounds only."""
        return cls(
            classes={Var: ((), None)},
            data_classes={
                Var._ComponentDataClass: (("lb", "ub"), None),
            },
        )

    @classmethod
    def value(cls, only_not_fixed=False):
        """Returns a StoreSpec object to store variable values only."""
        if only_not_fixed:
            return cls(
                classes={Var: ((), None), BooleanVar: ((), None)},
                data_classes={
                    Var._ComponentDataClass: (("value",), _value_if_not_fixed),
                    BooleanVar._ComponentDataClass: (("value",), _value_if_not_fixed),
                },
            )
        return cls(
            classes={Var: ((), None), BooleanVar: ((), None)},
            data_classes={
                Var._ComponentDataClass: (("value",), None),
                BooleanVar._ComponentDataClass: (("value",), None),
            },
        )

    @classmethod
    def isfixed(cls):
        """Returns a StoreSpec object to store if variables are fixed."""
        return cls(
            classes={Var: ((), None), BooleanVar: ((), None)},
            data_classes={
                Var._ComponentDataClass: (("fixed",), None),
                BooleanVar._ComponentDataClass: (("fixed",), None),
            },
        )

    @classmethod
    def suffix(cls, suffix_filter=None):
        # need to include classes and data that can have suffixes, but don't
        # need to store any attributes
        return cls(
            classes={
                Var: ((), None),
                BooleanVar: ((), None),
                Constraint: ((), None),
                Expression: ((), None),
                Suffix: ((), None),
            },
            data_classes={
                Var._ComponentDataClass: ((), None),
                BooleanVar._ComponentDataClass: ((), None),
                Constraint._ComponentDataClass: ((), None),
                Expression._ComponentDataClass: ((), None),
            },
            suffix_filter=suffix_filter,
        )

    @classmethod
    def value_isfixed(cls, only_fixed):
        """
        Return a StoreSpec object to store variable values and if fixed.

        Args:
            only_fixed: Only load fixed variable values
        """
        if only_fixed:
            return cls(
                classes={Var: ((), None), BooleanVar: ((), None)},
                data_classes={
                    Var._ComponentDataClass: (("value", "fixed"), _only_fixed),
                    BooleanVar._ComponentDataClass: (("value", "fixed"), _only_fixed),
                },
            )
        else:
            return cls(
                classes={Var: ((), None), BooleanVar: ((), None)},
                data_classes={
                    Var._ComponentDataClass: (("value", "fixed"), None),
                    BooleanVar._ComponentDataClass: (("value", "fixed"), None),
                },
            )

    @classmethod
    def value_isfixed_isactive(cls, only_fixed):
        """
        Retur a StoreSpec object to store variable values, if variables are
        fixed and if components are active.

        Args:
            only_fixed: Only load fixed variable values
        """
        if only_fixed:
            return cls(
                classes={
                    Var: ((), None),
                    BooleanVar: ((), None),
                    Param: ((), None),
                    Constraint: (("active",), None),
                    Block: (("active",), None),
                },
                data_classes={
                    Var._ComponentDataClass: (("value", "fixed"), _only_fixed),
                    BooleanVar._ComponentDataClass: (("value", "fixed"), _only_fixed),
                    pyomo.core.base.param._ParamData: (("value",), None),
                    Constraint._ComponentDataClass: (("active",), None),
                    Block._ComponentDataClass: (("active",), None),
                },
            )
        else:
            return cls(
                classes={
                    Var: ((), None),
                    BooleanVar: ((), None),
                    Param: ((), None),
                    Constraint: (("active",), None),
                    Block: (("active",), None),
                },
                data_classes={
                    Var._ComponentDataClass: (("value", "fixed"), None),
                    BooleanVar._ComponentDataClass: (("value", "fixed"), None),
                    pyomo.core.base.param._ParamData: (("value",), None),
                    Constraint._ComponentDataClass: (("active",), None),
                    Block._ComponentDataClass: (("active",), None),
                },
            )


def _may_have_subcomponents(o):
    """
    Args:
        o: an object.
    Returns:
        True if the object has a callable component_objects method, otherwise
        False.
    """
    if hasattr(o, "component_objects"):
        if hasattr(o.component_objects, "__call__"):
            return True


def _write_component(sd, o, wts, count=None, lookup=None, suffixes=None):
    """
    Writes a component state to the save dictionary under a key given by the
    components name.

    Args:
        sd: dictionary to to save the object into, will create a key that is the
            object name (not fully qualified)
        o: object to save
        wts: a StoreSpec object indicating what object attributes to write
        count: count the number of Pyomo componets written also used for ids
        lookup: is a lookup table for compoent ids from components
        suffixes: is a list of suffixes, that we are delaying writing
    Returns:
        None
    """
    if lookup is None:
        lookup = {}
    if suffixes is None:
        suffixes = []

    # Get list of attributes to save, also returns ff, which is a filter
    # function and only used in reading stuff back in.
    alist, ff = wts.get_class_attr_list(o)
    if alist is None:
        return  # alist is none means skip this component type
    # Get the componet name, doesn't need to be fully qualified or unique because
    # we are storing the state in a hierarchy structure
    oname = o.getname(fully_qualified=False)
    # Create a dictionary for this component, if storing suffixes assign it
    # a sequential id number and create a lookup table that takes the component
    # and returns its id for use later in writing suffix data
    sd[oname] = {"__type__": str(type(o))}
    if Suffix in wts.classes:
        sd[oname]["__id__"] = count.count
    lookup[id(o)] = count.count  # used python id() here for efficency
    if count is not None:
        count.count += 1  # incriment the componet counter
    for a in alist:  # store the desired attributes
        if a in wts.write_cbs:
            if wts.write_cbs[a] is None:
                sd[oname][a] = getattr(o, a, None)
            else:
                sd[oname][a] = wts.write_cbs[a](o)
        else:
            sd[oname][a] = getattr(o, a, None)
    sd[oname]["data"] = {}  # create a dict for compoent data and subcomponents
    if isinstance(o, Suffix):  # if is a suffix, make a list and delay writing
        if wts.suffix_filter is None or oname in wts.suffix_filter:
            suffixes.append(
                {"sd": sd[oname]["data"], "o": o, "wts": wts, "lookup": lookup}
            )
    else:  # if not suffix go ahead and write the data
        _write_component_data(
            sd=sd[oname]["data"],
            o=o,
            wts=wts,
            lookup=lookup,
            count=count,
            suffixes=suffixes,
        )


def _write_component_data(sd, o, wts, count=None, lookup=None, suffixes=None):
    """
    Iterate through the component data and write to the sd dictionary. The keys
    for the data items are added to the dictionary. If the component has
    subcomponents they are written by a recursive call to _write_component under
    the __pyomo_components__ key.

    Args:
        sd: dictionary to to save the object into, will create keys that are the
            data object indexes repn.
        o: object to save
        wts: a StoreSpec object indicating what object attributes to write
        count: count the number of Pyomo componets written also used for ids
        lookup: is a lookup table for compoent ids from components
        suffixes: is a list of suffixes, that we are delaying writing
    Returns:
        None
    """
    if lookup is None:
        lookup = {}
    if suffixes is None:
        suffixes = []

    if isinstance(o, Suffix) and Suffix in wts.classes:
        # make special provision for writing suffixes.
        for key in o:
            el = o[key]
            if id(key) not in lookup:
                # didn't store these compoents so can't write suffix.
                continue
            sd[lookup[id(key)]] = el  # Asssume keys are Pyomo model components
    else:  # rest of compoents with normal componet data structure
        frst = True  # on first item when true
        try:
            item_keys = o.keys()
        except AttributeError:
            item_keys = [None]
        for key in item_keys:
            if key is None and isinstance(o, ComponentData):
                el = o
            else:
                el = o[key]
            if frst:  # assume all item are same type, use first to get alist
                alist, ff = wts.get_data_class_attr_list(el)  # get attributes
                if alist is None:
                    return  # if None then skip writing
            frst = False  # done with first only stuff
            edict = {"__type__": str(type(el))}
            if Suffix in wts.classes:  # if writing suffixes give data compoents an id
                edict["__id__"] = count.count
                lookup[id(el)] = count.count  # and add to lookup table
            if count is not None:
                count.count += 1  # inciment component counter
            sd[repr(key)] = edict  # stick item dict into component data dict
            for a in alist:  # store desired attributes
                if a in wts.write_cbs:
                    if wts.write_cbs[a] is None:
                        edict[a] = getattr(el, a)
                    else:
                        edict[a] = wts.write_cbs[a](el)
                else:
                    edict[a] = getattr(el, a)
            hascomps = False  # Has sub-components (like a Block would have)
            if _may_have_subcomponents(el):  # block or block like component
                for o2 in el.component_objects(descend_into=False):
                    # loop through sub-components
                    if not hascomps:  # if here it does have sub-components
                        cdict = {}  # so store those in __pyomo_components__
                        edict["__pyomo_components__"] = cdict
                    hascomps = True  # only make __pyomo_components__ dict once
                    # write each sub-component
                    _write_component(
                        sd=cdict,
                        o=o2,
                        wts=wts,
                        count=count,
                        lookup=lookup,
                        suffixes=suffixes,
                    )


def component_data_to_dict(o, wts):
    """
    Component data to a dict.
    """
    el = o
    alist, ff = wts.get_data_class_attr_list(el)  # get attributes
    if alist is None:
        return  # if None then skip writing
    edict = {}  # if not writing suffixes don't need ids
    for a in alist:  # store desired attributes
        edict[a] = getattr(el, a)
    hascomps = False  # Has sub-components (like a Block would have)
    if _may_have_subcomponents(el):  # block or block like component
        for o2 in el.component_objects(descend_into=False):
            # loop through sub-components
            if not hascomps:  # if here it does have sub-components
                cdict = {}  # so store those in __pyomo_components__
                edict["__pyomo_components__"] = cdict
            hascomps = True  # only make __pyomo_components__ dict once
            # write each sub-component
            _write_component(sd=cdict, o=o2, wts=wts)
    return edict


def to_json(
    o,
    fname=None,
    human_read=False,
    wts=None,
    metadata=None,
    gz=None,
    return_dict=False,
    return_json_string=False,
):
    """
    Save the state of a model to a Python dictionary, and optionally dump it
    to a json file.  To load a model state, a model with the same structure must
    exist.  The model itself cannot be recreated from this.

    Args:
        o: The Pyomo component object to save.  Usually a Pyomo model, but could
            also be a subcomponent of a model (usually a sub-block).
        fname: json file name to save model state, if None only create
            python dict
        gz: If fname is given and gv is True gzip the json file. The default is
            True if the file name ends with '.gz' otherwise False.
        human_read: if True, add indents and spacing to make the json file more
            readable, if false cut out whitespace and make as compact as
            possilbe
        metadata: A dictionary of additional metadata to add.
        wts: is What To Save, this is a StoreSpec object that specifies what
            object types and attributes to save.  If None, the default is used
            which saves the state of the compelte model state.
        metadata: additional metadata to save beyond the standard format_version,
            date, and time.
        return_dict: default is False if true returns a dictionary representation
        return_json_string: default is False returns a json string

    Returns:
        If return_dict is True returns a dictionary serialization of the Pyomo
        component.  If return_dict is False and return_json_string is True
        returns a json string dump of the dict.  If fname is given the dictionary
        is also written to a json file.  If gz is True and fname is given, writes
        a gzipped json file.
    """
    if gz is None:
        if isinstance(fname, str):
            gz = fname.endswith(".gz")
        else:
            gz = False
    if metadata is None:
        metadata = {}

    suffixes = []
    lookup = {}
    count = Counter()
    start_time = time.time()
    if wts is None:
        wts = StoreSpec()
    now = datetime.datetime.now()
    sd = {
        "__metadata__": {
            "format_version": __format_version__,
            "date": datetime.date.isoformat(now.date()),
            "time": datetime.time.isoformat(now.time()),
            "other": metadata,
        }
    }
    # first write the component
    _write_component(sd, o, wts, count, suffixes=suffixes, lookup=lookup)
    for s in suffixes:
        _write_component_data(**s)
    pdict = {}
    sd["__metadata__"]["__performance__"] = pdict
    pdict["n_components"] = count.count
    dict_time = time.time()
    pdict["etime_make_dict"] = dict_time - start_time
    # This returns the dict but if fname is specified also save to json file
    dump_kw = {"indent": 2} if human_read else {"separators": (",", ":")}
    if fname is not None:
        if gz:
            with gzip.open(fname, "w") as f:
                f.write(json.dumps(sd, **dump_kw).encode("utf-8"))
        else:
            with open(fname, "w") as f:
                json.dump(sd, f, **dump_kw)
    file_time = time.time()
    # unfortunatly I can't write how long it took to write the file in the file
    pdict["etime_write_file"] = file_time - dict_time
    if return_dict:
        # In interactive environments returning the dict can cause it to print
        # an extremely large amount of stuff.  So added this option to make sure
        # it's really what you want.
        return sd
    elif return_json_string:
        return json.dumps(sd, **dump_kw)
    else:
        return None


def _read_component(sd, o, wts, lookup=None, suffixes=None, root_name=None):
    """
    Read a component dictionary into a model
    """
    if lookup is None:
        lookup = {}
    if suffixes is None:
        suffixes = {}

    alist, ff = wts.get_class_attr_list(o)
    if alist is None:
        return
    if root_name is None:
        oname = o.getname(fully_qualified=False)
    else:
        oname = root_name
    try:
        odict = sd[oname]
    except KeyError as e:
        if wts.ignore_missing:
            return
        else:
            raise (e)
    if ff is not None:
        alist = ff(o, odict)
    if Suffix in wts.classes:
        lookup[odict["__id__"]] = o
    for a in alist:
        try:
            if a in wts.read_cbs:
                if wts.read_cbs[a] is None:
                    pass
                else:
                    wts.read_cbs[a](o, odict[a])
            else:
                setattr(o, a, odict[a])
        except KeyError as e:
            if wts.ignore_missing:
                return
            else:
                raise (e)
    if isinstance(o, Suffix):
        if wts.suffix_filter is None or oname in wts.suffix_filter:
            suffixes[odict["__id__"]] = odict["data"]  # is populated
    else:  # read nonsufix component data
        _read_component_data(odict["data"], o, wts, lookup=lookup, suffixes=suffixes)


def _read_component_data(sd, o, wts, lookup=None, suffixes=None):
    """
    Read a Pyomo component's data in from a dict.

    Args:
        sd: dictionary to read from
        o: Pyomo component whoes data to read
        wts: StoreSpec object specifying what to read in
        lookup: a lookup table for id to componet for reading suffixes
        suffixes: a list of suffixes put off reading until end

    Returns:
        None
    """
    if lookup is None:
        lookup = {}
    if suffixes is None:
        suffixes = {}

    alist = []  # list of attributes to read
    c = 0  # counter of data items in component
    try:
        item_keys = o.keys()
    except AttributeError:
        item_keys = [None]
    for key in item_keys:
        if key is None and isinstance(o, ComponentData):
            el = o
        else:
            el = o[key]
        if c == 0:  # if first data item assume all itmes are same and get alist
            alist, ff = wts.get_data_class_attr_list(el)  # ff is fileter function
            if alist is None:
                return  # skip reading this type
        c += 1
        try:  # get data from dict
            edict = sd[repr(key)]
        except KeyError as e:  # data was missing either ignore or raise except
            if wts.ignore_missing:
                return  # if ignore missing option its okay
            else:
                raise (e)  # else raise exception
        if ff is not None:  # if a filer function was given, use it to make a
            # new a list based on the model and whats stored for the state
            # this lets you contionally load things, for example only load
            # values for unfixed variables.
            alist = ff(el, edict)
        if Suffix in wts.classes:  # if loading suffixes make lookup table id
            lookup[edict["__id__"]] = el
        for a in alist:  # read in desired attributes
            try:
                if a in wts.read_cbs:
                    if wts.read_cbs[a] is None:
                        pass
                    else:
                        wts.read_cbs[a](el, edict[a])
                else:  # directly set an attribute
                    setattr(el, a, edict[a])
            except KeyError as e:  # attribute missing
                if wts.ignore_missing:
                    return  # if ignore option then is okay
                else:
                    raise (e)  # otherwise raise an exception
        if _may_have_subcomponents(el) and "__pyomo_components__" in edict:
            # read sub-components of block-like
            for o2 in el.component_objects(descend_into=False):
                # recursive read here
                _read_component(
                    edict["__pyomo_components__"],
                    o2,
                    wts,
                    lookup=lookup,
                    suffixes=suffixes,
                )


def component_data_from_dict(sd, o, wts):
    """
    Component data to a dict.
    """
    el = o
    alist = []  # list of attributes to read
    alist, ff = wts.get_data_class_attr_list(el)  # ff is fileter function
    if alist is None:
        return
    edict = sd
    if ff is not None:
        alist = ff(o, edict)
    for a in alist:  # read in desired attributes
        try:
            if a in wts.read_cbs:  # use a callback
                wts.read_cbs[a](el, edict[a])
            else:  # directly set an attribute
                setattr(el, a, edict[a])
        except KeyError as e:  # attribute missing
            if wts.ignore_missing:
                return  # if ignore option then is okay
            else:
                raise (e)  # otherwise raise an exception
    if _may_have_subcomponents(el):  # read sub-components of block-like
        for o2 in el.component_objects(descend_into=False):
            # recursive read here
            _read_component(edict["__pyomo_components__"], o2, wts)


def _read_suffixes(lookup, suffixes):
    """
    Go through the list of suffixes and read the data back in.

    Args:
        lookup: a lookup table to go from id to component
        suffixes: a dictionary with suffix id keys and value dict value
    Returns:
        None
    """
    for uid in suffixes:
        d = suffixes[uid]
        s = lookup[uid]  # suffixes keys are ids, so get suffix component
        for key in d:  # set values from value dict
            try:
                kc = lookup[int(key)]  # use int because json turn keys to string
            except KeyError:
                continue
            s[kc] = d[key]


def from_json(o, sd=None, fname=None, s=None, wts=None, gz=None, root_name=None):
    """
    Load the state of a Pyomo component state from a dictionary, json file, or
    json string.  Must only specify one of sd, fname, or s as a non-None value.
    This works by going through the model and loading the state of each
    sub-compoent of o. If the saved state contains extra information, it is
    ignored.  If the save state doesn't contain an enetry for a model component
    that is to be loaded an error will be raised, unless ignore_missing = True.

    Args:
        o: Pyomo component to for which to load state
        sd: State dictionary to load, if None, check fname and s
        fname: JSON file to load, only used if sd is None
        s: JSON string to load only used if both sd and fname are None
        wts: StoreSpec object specifying what to load
        gz: If True assume the file specified by fname is gzipped. The default is
            True if fname ends with '.gz' otherwise False.

    Returns:
        Dictionary with some perfomance information. The keys are
        "etime_load_file", how long in seconds it took to load the json file
        "etime_read_dict", how long in seconds it took to read models state
        "etime_read_suffixes", how long in seconds it took to read suffixes
    """
    if gz is None:
        if isinstance(fname, str):
            gz = fname.endswith(".gz")
        else:
            gz = False

    # keeping track of elapsed time.  want to make sure I don't do anything
    # that's too slow.
    start_time = time.time()
    # Get the model state dict from one of three sources
    if sd is not None:  # Existing Python dict (for in-memory stuff).
        pass
    elif fname is not None:  # Read in from a json file
        if gz:
            with gzip.open(fname, "r") as f:
                fr = f.read()
                sd = json.loads(fr)
        else:
            with open(fname, "r") as f:
                sd = json.load(f)  # json file
    elif s is not None:  # Use a json string (not really sure if useful)
        sd = json.loads(s)  # json string
    else:  # Didn't specify at least one source
        raise Exception("Need to specify a data source to load from")
    dict_time = time.time()  # To calculate how long it took to read file
    if wts is None:  # if no StoreSpec object given use the default, which should
        wts = StoreSpec()  # be the typlical save everything important
    lookup = {}  # A dict to use for a lookup tables
    suffixes = {}  # A list of suffixes delayed to end so lookup is complete
    # Read toplevel componet (is recursive)
    if root_name is None:
        for k in sd:
            if k.startswith("__") and k.endswith("__"):
                # This is metadata or maybe some similar future addition.
                continue
            else:
                root_name = k
                break  # should be one root, use it's name
    _read_component(sd, o, wts, lookup=lookup, suffixes=suffixes, root_name=root_name)
    read_time = time.time()  # to calc time to read model state minus suffixes
    # Now read in the suffixes
    _read_suffixes(lookup, suffixes)
    suffix_time = time.time()  # to calculate time to read suffixes
    pdict = {}  # return some perfomance information, to make sure not too slow
    pdict["etime_load_file"] = dict_time - start_time
    pdict["etime_read_dict"] = read_time - dict_time
    pdict["etime_read_suffixes"] = suffix_time - read_time
    return pdict
