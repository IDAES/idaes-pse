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
from abc import abstractmethod
from itertools import product

from pyomo.core.base.param import SimpleParam
from pyomo.opt.results import SolutionStatus

from .pyomo_modeling import *
from ..materials.design import Design


class IndexedElem(object):
    """Base class for indexed MatOpt objects.

    Should not be necessary for users to instantiate, but developers will
    utilize the constructors (especially fromComb) and mask method when
    creating new classes of Expression and Rule objects.

    Attributes:
        sites (list<int>): The sites that the object is indexed over
        bonds (list<int>): The bonds that the object is indexed over
        site_types (list<int>): The site types that the object is indexed over
        bond_types (list<int>): The bond types that the object is indexed over
        confs (list<int>): The conformations that the object is indexed over
    """

    # === STANDARD CONSTRUCTOR
    def __init__(
        self, sites=None, bonds=None, site_types=None, bond_types=None, confs=None
    ):
        """Standard constructor of IndexedElem.

        Args:
            sites (list<int>): The sites that the object is indexed over
            bonds (list<tuple<int,int>>): The bonds that the object is indexed
                over
            site_types (list<BBlock>): The site types that the object is indexed
                over
            bond_types (list<tuple<BBlock,BBlock>>): The bond types that the
                object is indexed over
            confs (list<int>): The conformations that the object is indexed over
        """

        self._sites = sites
        self._bonds = bonds
        self._site_types = site_types
        self._bond_types = bond_types
        self._confs = confs

    # === CONSTRUCTOR - From combinations of other IndexedElem objecs
    @classmethod
    def fromComb(cls, *args):
        """
        Constructor of IndexedElem from other IndexedElem objects.

        This constructor creates a new indexing object by combining two
        other indexed objects. The various attributes are mixed according
        to three simple cases for each type of index:

        Case #1: Both objects are indexed over that index.
            The intersection of the two lists of indices is used.
        Case #2: Only one object has that index.
            The resulting object gains the index.
        Case #3: Both objects are not indexed over that index.
            The resulting object is also not indexed over that index.

        This is simply meant to reproduce the way we use set notation for
        writing products of expressions and variables in constraints.

        Args:
            *args (list<IndexedElem>): Objects to combine indices for.

        Returns:
            (IndexedElem) Index object from combination of indices.
        """
        LHS, RHS, *args = args
        Comb = cls._fromComb2(LHS, RHS)
        if len(args) > 0:
            return cls.fromComb(Comb, *args)
        else:
            return Comb

    # === AUXILIARY METHODS
    @classmethod
    def _fromComb2(cls, LHS, RHS):
        """Constructor of an IndexedElem from two other IndexedElem objects.

        Utilized in the recursive fromComb method, above.

        Args:
            LHS (IndexedElem): Object to combine indices for.
            RHS (IndexedElem): Object to combine indices for.

        Returns:
            (IndexedElem) Index object from combination of indices.
        """

        if LHS.sites is not None and RHS.sites is not None:
            sites = list(set(LHS.sites) & set(RHS.sites))
        else:
            sites = LHS.sites if LHS.sites is not None else RHS.sites
        if LHS.bonds is not None and RHS.bonds is not None:
            bonds = list(set(LHS.bonds) & set(RHS.bonds))
        else:
            bonds = LHS.bonds if LHS.bonds is not None else RHS.bonds
        if LHS.site_types is not None and RHS.site_types is not None:
            site_types = list(set(LHS.site_types) & set(RHS.site_types))
        else:
            site_types = (
                LHS.site_types if LHS.site_types is not None else RHS.site_types
            )
        if LHS.bond_types is not None and RHS.bond_types is not None:
            bond_types = list(set(LHS.bond_types) & set(RHS.bond_types))
        else:
            bond_types = (
                LHS.bond_types if LHS.bond_types is not None else RHS.bond_types
            )
        if LHS.confs is not None and RHS.confs is not None:
            confs = list(set(LHS.confs) & set(RHS.confs))
        else:
            confs = LHS.confs if LHS.confs is not None else RHS.confs
        return cls(
            sites=sites,
            bonds=bonds,
            site_types=site_types,
            bond_types=bond_types,
            confs=confs,
        )

    def mask(self, index, Comb):
        """
        Method to identify the indexes relevant to this object.

        Given an instance of index that was generated by another IndexedElem
        object (Comb), we identify which parts of that index were relevant
        to this object.

        Example::

            VarIndexes = IndexedElem(sites=[1,2])
            CoefIndexes = IndexedElem(site_types=['A','B'])
            Comb = IndexedElem.fromComb(VarIndexes,CoefIndexes)
            for k in Comb.keys():
                site = VarIndexes.mask(k,Comb)
                site_type = CoefIndexes.mask(k,Comb)

        Args:
            index (tuple<int/BBlock>): index from which to identify relevant parts
            Comb (IndexedElem): object from which the index was generated

        Returns:
            (tuple<int/BBlock>) index with indices relevant to this object remaining
        """
        if Comb.sites is not None:
            i, *index = index
        if Comb.bonds is not None:
            i, j, *index = index
        if Comb.site_types is not None:
            k, *index = index
        if Comb.bond_types is not None:
            k, l, *index = index
        if Comb.confs is not None:
            c, *index = index
        result = []
        if self.sites is not None:
            result.append(i)
        if self.bonds is not None:
            result.append(i)
            result.append(j)
        if self.site_types is not None:
            result.append(k)
        if self.bond_types is not None:
            result.append(k)
            result.append(l)
        if self.confs is not None:
            result.append(c)
        if not result:
            result = [None]
        return tuple(result)

    @property
    def dims(self):
        """Relevant dimensions of indices.

        Returns:
            (list<bool>) flags to indicate which index sets are relevant.
        """

        return (
            self.sites is not None,
            self.bonds is not None,
            self.site_types is not None,
            self.bond_types is not None,
            self.confs is not None,
        )

    @property
    def index_sets(self):
        """Sets (actually lists) of indices.

        Note that in the cases that there are no relevant indices, a
        dummy list [[None]] is returned to allow the [None] key to be
        included.

        Returns:
            (list<list<int/BBlock>>) lists of indices of each relevant type.
        """

        result = [
            s
            for s in (
                self.sites,
                self.bonds,
                self.site_types,
                self.bond_types,
                self.confs,
            )
            if s is not None
        ]
        if not result:
            result = [[None]]
        return result

    @property
    def index_dict(self):
        """Dictionary of relevant attributes.

        Returns:
            (dict<string:list<int/BBlock>>) attributes of this IndexedElem.
        """

        return {
            "sites": self.sites,
            "bonds": self.bonds,
            "site_types": self.site_types,
            "bond_types": self.bond_types,
            "confs": self.confs,
        }

    @property
    def sites(self):
        """List of sites relevant to this object."""
        return self._sites

    @property
    def bonds(self):
        """List of bonds relevant to this object."""
        return self._bonds

    @property
    def site_types(self):
        """List of site types relevant to this object."""
        return self._site_types

    @property
    def bond_types(self):
        """List of bond types relevant to this object."""
        return self._bond_types

    @property
    def confs(self):
        """List of conformation types relevant to this object."""
        return self._confs

    def keys(self):
        """Method creating a generator for the keys relevant to this object.

        Note that the [None] key is returned in the case of an object
        not indexed by any of the possible index types. In other words,
        [None] is the key for a scalar variable/expression/coefficient.

        Returns:
            (generator<tuple<int/BBlock>>) keys generator (similar to dict.keys()).
        """

        index_sets = self.index_sets
        if len(index_sets) > 1:
            return product(*self.index_sets)
        elif len(index_sets) == 1:
            return (k for k in index_sets[0])
        else:
            raise NotImplementedError("There should always be at least " "a [None] key")


class Coef(IndexedElem):
    """A class for coefficient data indexed over material index sets.

    This class is useful for representing indexed data and automatically
    generating complex indexed expressions. For example, the multiplication
    of bond-type coefficients with bond-indexed variables would not be
    easily representable via standard Python objects.

    The key benefit is that these objects have the useful IndexedElem methods
    and also allow access of data via getitem operators.

    Attributes:
    vals (dict/list) data structure of coefficient values.
    (index information inherited from IndexedElem)
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, vals, **kwargs):
        """Standard constructor of indexed coefficients.

        Args:
            vals (list<float>/dict/other): Any data structure that supports the
                __gettitem__ method for keys generated by this object's
                IndexedElem.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """

        self.vals = vals
        IndexedElem.__init__(self, **kwargs)

    # === BASIC QUERY METHODS
    def __getitem__(self, k):
        """Method to access coefficient values by bracket operator"""
        return self.vals[k]


class Expr(IndexedElem):
    """An abstract class for representing expressions when building rules.

    The key benefit of this class is that we utilize the useful methods of
    IndexedElem and establish the interface for derived expressions.

    Expressions can be generated over a subset of the design space
    (i.e., only for some combinations of sites, site-types, etc.) by
    providing keywords that are passed to the constructor of IndexedElem.
    Else, the relevant indexes are infered from the expression components.

    Attributes:
        (index information inherited from IndexedElem)
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, **kwargs):
        """Standard constructor for abstract class Expr.

        Args:
        **kwargs: Optional, index information passed to IndexedElem if
            interested in a subset of indices
            Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        IndexedElem.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    @abstractmethod
    def _pyomo_expr(self, index=None):
        """Abstract interface for generating Pyomo expressions.

        Args:
        index (list): Optional, index to to create an instance of a Pyomo
            expression. In the case of a scalar, the valid index is None.

        Returns:
        An instance of a Pyomo expression.
        """
        raise NotImplementedError


class LinearExpr(Expr):
    """A class for representing simple expressions of site descriptors.

    The use of this class is to generate expressions from multiplication and
    summation of coefficients and descriptors. Importantly, expressions
    of this type maintain the same indexing of their component descriptors
    and coefficients. Summation is taken across multiple descriptors, not
    multiple instances of indexes.

    Attributes:
        coefs (float/list<float>): coefficient to multiply each descriptor by
        descs (Descriptor/list<Descriptor>): descriptors to add together
        offset (float/int): scalar value to add to the expression
            (index information inherited from IndexedElem)
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, descs=None, coefs=1.0, offset=0.0, **kwargs):
        """Standard constructor for linear expression objects.

        Args:
            coefs (float/list<float>): Optional, coefficient to multiply
                each descriptor by. Default: 1.0
            descs (Descriptor/list<Descriptor>): descriptors to add
            offset (float/int): Optional, scalar value to add to the expression
                Default: 0.0
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.coefs = coefs if type(coefs) is list else [coefs]
        self.descs = descs if type(descs) is list else [descs]
        self.offset = offset
        if descs is not None:
            if type(descs) is list:
                Comb = IndexedElem.fromComb(*descs)
                kwargs = {**Comb.index_dict, **kwargs}
            else:
                kwargs = {**descs.index_dict, **kwargs}
        Expr.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_expr(self, index=None):
        """Interface for generating Pyomo expressions.

        Args:
            index (list): Optional, index to to create an instance of a Pyomo
                expression. In the case of a scalar, the valid index is None.

        Returns:
            An instance of a Pyomo expression.
        """
        result = self.offset
        for it, desc in enumerate(self.descs):
            if desc is not None and self.coefs[it] is not None:
                result += self.coefs[it] * desc._pyomo_var[index]
        return result


class SiteCombination(Expr):
    """A class for representing summations of descriptors at two sites.

    Attributes:
        coefi (float/list<float>): coefficients at the first site
        desci (Descriptor/Expr): descriptor or expression for the first site
        coefj (float/list<float>): coefficients at the second site
        descj (Descriptor/Expr): descriptor or expression for the second site
        offset (float): scalar coefficient to add to the rest of the expression
        symmetric_bonds (bool): flag to indicate if site combinations should be
            considered symmetric (and therefore, should only generate half as
            many terms) (index information inherited from IndexedElem)
    """

    # === STANDARD CONSTRUCTOR
    def __init__(
        self,
        coefi,
        desci,
        coefj=None,
        descj=None,
        offset=0.0,
        symmetric_bonds=False,
        **kwargs
    ):
        """Standard constructor for site combination expressions.

        Args:
            coefi (float/list<float>): coefficients at the first site.
                Default: 1.0
            desci (Descriptor/Expr): term for first site
            coefj (float/list<float>): Optional, coefficients at the second site.
                Default: Equal to the coefficient for the first site.
            descj (Descriptor/Expr): Optional, term for the second site.
                Default: Equal to the descriptor for the first site.
            offset (float): Optional, scalar coefficient to add to expression.
                Default: 1.0
            symmetric_bonds (bool): Optional, flag to indicate if combinations
                are symmetric.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.coefi = coefi
        self.desci = desci
        self.coefj = coefj if coefj is not None else coefi
        self.descj = descj if descj is not None else desci
        self.offset = offset
        self.symmetric_bonds = symmetric_bonds
        if "bonds" not in kwargs:
            kwargs["bonds"] = [
                (i, j)
                for i in desci.sites
                for j in desci.canv.NeighborhoodIndexes[i]
                if (j is not None and (not symmetric_bonds or j > i))
            ]
        if "bond_types" in kwargs:
            pass  # use the kwargs bond_types
        elif coefi.bond_types is not None:
            kwargs["bond_types"] = coefi.bond_types
        elif desci.site_types is not None:
            kwargs["bond_types"] = [
                (k, l) for k in desci.site_types for l in desci.site_types
            ]
        Expr.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_expr(self, index=None):
        """Interface for generating Pyomo expressions.

        Args:
            index (list): Optional, index to to create an instance of a Pyomo
                expression. In the case of a scalar, the valid index is None.

        Returns:
            An instance of a Pyomo expression.
        """
        if len(index) == 4:
            i, j, k, l = index
        elif len(index) == 2:
            i, j, k, l = index[0], index[1], (), ()
        else:
            raise NotImplementedError(
                "Decide how to split the extra " "indices in this case..."
            )
        if (
            type(self.coefi) is float
            or type(self.coefi) is int
            or type(self.coefi) is SimpleParam
        ):
            ci = self.coefi
        else:
            coefi_index = self.coefi.mask(
                (i, i, j, k, k, l),
                IndexedElem(
                    sites=[i], bonds=[(i, j)], site_types=[k], bond_types=[(k, l)]
                ),
            )
            ci = self.coefi[coefi_index]

        if (
            type(self.coefj) is float
            or type(self.coefj) is int
            or type(self.coefj) is SimpleParam
        ):
            cj = self.coefj
        else:
            coefj_index = self.coefj.mask(
                (j, j, i, l, l, k),
                IndexedElem(
                    sites=[j], bonds=[(j, i)], site_types=[l], bond_types=[(l, k)]
                ),
            )
            cj = self.coefj[coefj_index]
        desci_index = self.desci.mask(
            (i, i, j, k, k, l),
            IndexedElem(sites=[i], bonds=[(i, j)], site_types=[k], bond_types=[(k, l)]),
        )
        descj_index = self.descj.mask(
            (j, j, i, l, l, k),
            IndexedElem(sites=[j], bonds=[(j, i)], site_types=[l], bond_types=[(l, k)]),
        )
        di = self.desci._pyomo_expr(index=desci_index)
        dj = self.descj._pyomo_expr(index=descj_index)
        return self.offset + ci * di + cj * dj


class SumNeighborSites(Expr):
    """A class for expressions for summation across neighbor sites.

    Attributes:
        desc (Descriptor): descriptors to sum around a site
        coefs (float/list<float>): Optional, coefficients to multiple each
            neighbor term by. Default=1.0
        offset: Optional, term to add to the expression.
            Default=0.0 (index information inherited from IndexedElem)
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, desc, coefs=1.0, offset=0.0, **kwargs):
        """Standard constructor for expressions of neighbor summations.

        Args:
            desc (Descriptor): descriptors to sum around a site
            coefs (float/list<float>): Optional, coefficients to multiple each
                neighbor term by. Default=1.0
            offset: Optional, term to add to the expression.
                Default=0.0
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.desc = desc
        self.coefs = coefs
        self.offset = offset
        kwargs = {**desc.index_dict, **kwargs}
        Expr.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_expr(self, index=None):
        """Interface for generating Pyomo expressions.

        Args:
        index (list): Optional, index to to create an instance of a Pyomo
            expression. In the case of a scalar, the valid index is None.

        Returns:
        An instance of a Pyomo expression.
        """
        i, *index = index
        if index == (None,):
            index = ()
        result = self.offset
        for n, j in enumerate(self.desc.canv.NeighborhoodIndexes[i]):
            if j is not None:
                result += (
                    self.coefs
                    if (
                        type(self.coefs) is float
                        or type(self.coefs) is int
                        or type(self.coefs) is SimpleParam
                    )
                    else self.coefs[n]
                ) * self.desc._pyomo_var[(j, *index)]
        return result


class SumNeighborBonds(Expr):
    """A class for expressions from summation of neighbor bond descriptors.

    Attributes:
        desc (Descriptor/Expr): descriptors to sum over
        coefs (float/list<float>): coefficients to multiply bonds to neighbor sites
        offset (float): coefficient to add to the expression
        symmetric_bonds (bool): flag to indicate if bond variables should be
            added in a symmetric way (index information inherited from IndexedElem)
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, desc, coefs=1.0, offset=0.0, symmetric_bonds=False, **kwargs):
        """Standard constructor for summation of neighboring bonds

        Args:
            desc (Descriptor): descriptors to sum around a site
            coefs (float/list<float>): Optional, coefficients to multiple each
                neighbor term by. Default=1.0
            offset (float): Optional, coefficient to add to the expression.
                Default=0.0
            symmetric_bonds (bool): Optional, flag to indicate if bond variables
                should be considered symmetric.
                Default=False
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.desc = desc
        self.coefs = coefs
        self.offset = offset
        self.symmetric_bonds = symmetric_bonds
        kwargs = {**desc.index_dict, **kwargs}
        Expr.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_expr(self, index=None):
        """Interface for generating Pyomo expressions.

        Args:
        index (list): Optional, index to to create an instance of a Pyomo
            expression. In the case of a scalar, the valid index is None.

        Returns:
        An instance of a Pyomo expression.
        """
        i, *index = index
        if index == (None,):
            index = ()
        result = self.offset
        for n, j in enumerate(self.desc.canv.NeighborhoodIndexes[i]):
            if j is not None:
                if self.symmetric_bonds:
                    i, j = min(i, j), max(i, j)
                result += (
                    self.coefs
                    if (
                        type(self.coefs) is float
                        or type(self.coefs) is int
                        or type(self.coefs) is SimpleParam
                    )
                    else self.coefs[n]
                ) * self.desc._pyomo_expr(index=(i, j, *index))
        return result


class SumSites(Expr):
    """A class for expressions formed by summation over canvas sites.

    Attributes:
        desc (Descriptor/Expr): descriptors or expressions to sum over
        coefs (float/list<float>): coefficients to multiply contributions
            from each site
        offset (float): coefficient to add to the expression
        sites_to_sum (list<int>): sites to consider in the summation
            (index information inherited from IndexedElem)
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, desc, coefs=1.0, offset=0.0, sites_to_sum=None, **kwargs):
        """Standard constructor for summation of site contributions.

        Args:
            desc (Descriptor): descriptors or expressions to sum across all sites
            coefs (float/list<float>): Optional, coefficients to multiple each
                site term by. Default=1.0
            offset (float): Optional, coefficient to add to the expression.
                Default=0.0
            sites_to_sum (list<int>): Optional, subset of canvas sites to sum.
                Default=None, meaning all sites in the desc object are considered.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.desc = desc
        self.coefs = coefs
        self.offset = offset
        self.sites_to_sum = sites_to_sum if sites_to_sum is not None else desc.sites
        kwargs = {**desc.index_dict, **kwargs}
        kwargs.pop("sites")
        Expr.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_expr(self, index=None):
        """Interface for generating Pyomo expressions.

        Args:
            index (list): Optional, index to to create an instance of a Pyomo
                expression. In the case of a scalar, the valid index is None.

        Returns:
            An instance of a Pyomo expression.
        """
        if index == (None,):
            index = ()
        result = self.offset
        for i in self.sites_to_sum:
            result += (
                self.coefs
                if (
                    type(self.coefs) is float
                    or type(self.coefs) is int
                    or type(self.coefs) is SimpleParam
                )
                else self.coefs[(i, *index)]
            ) * self.desc._pyomo_var[(i, *index)]
        return result


class SumBonds(Expr):
    """A class for expressions formed by summation over canvas bonds.

    Attributes:
        desc (Descriptor/Expr): descriptors or expressions to sum over
        coefs (float/list<float>): coefficients to multiply contributions
            from each bond
        offset (float): coefficient to add to the expression
        bonds_to_sum (list<int>): bonds to consider in the summation
            (index information inherited from IndexedElem)
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, desc, coefs=1.0, offset=0.0, bonds_to_sum=None, **kwargs):
        """Standard constructor for summation of bond contributions.

        Args:
            desc (Descriptor): descriptors or expressions to sum across all bonds
            coefs (float/list<float>): Optional, coefficients to multiple each
                bond term by. Default=1.0
            offset (float): Optional, coefficient to add to the expression.
                Default=0.0
            bonds_to_sum (list<int>): Optional, subset of canvas bonds
                (i.e., neighbor connections) to sum.
                Default=None, meaning all bonds in the desc object are considered.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.desc = desc
        self.coefs = coefs
        self.offset = offset
        self.bonds_to_sum = bonds_to_sum if bonds_to_sum is not None else desc.bonds
        kwargs = {**desc.index_dict, **kwargs}
        kwargs.pop("bonds")
        Expr.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_expr(self, index=None):
        """Interface for generating Pyomo expressions.

        Args:
            index (list): Optional, index to to create an instance of a Pyomo
                expression. In the case of a scalar, the valid index is None.

        Returns:
            An instance of a Pyomo expression.
        """
        if index == (None,):
            index = ()
        result = self.offset
        for i, j in self.bonds_to_sum:
            result += (
                self.coefs
                if (
                    type(self.coefs) is float
                    or type(self.coefs) is int
                    or type(self.coefs) is SimpleParam
                )
                else self.coefs[(i, j, *index)]
            ) * self.desc._pyomo_var[(i, j, *index)]
        return result


class SumSiteTypes(Expr):
    """A class for expressions formed by summation over building block types.

    Attributes:
        desc (Descriptor/Expr): descriptors or expressions to sum over
        coefs (float/list<float>): coefficients to multiply contributions
            from each building block type
        offset (float): coefficient to add to the expression
        site_types_to_sum (list<BBlock>): building block types to consider in
            the summation (index information inherited from IndexedElem)
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, desc, coefs=1.0, offset=0.0, site_types_to_sum=None, **kwargs):
        """Standard constructor for summation of contributions by site-type.

        Args:
            desc (Descriptor): descriptors or expressions to sum across site types
            coefs (float/list<float>): Optional, coefficients to multiple each
                site-type term by. Default=1.0
            offset (float): Optional, coefficient to add to the expression.
                Default=0.0
            bonds_types_to_sum (list<int>): Optional, subset of site types
                to sum. Default=None, meaning all site-types in the desc object are
                considered.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.desc = desc
        self.coefs = coefs
        self.offset = offset
        self.site_types_to_sum = (
            site_types_to_sum if site_types_to_sum is not None else desc.site_types
        )
        kwargs = {**desc.index_dict, **kwargs}
        kwargs.pop("site_types")
        Expr.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_expr(self, index=None):
        """Interface for generating Pyomo expressions.

        Args:
            index (list): Optional, index to to create an instance of a Pyomo
                expression. In the case of a scalar, the valid index is None.

        Returns:
            An instance of a Pyomo expression.
        """
        assert index is not None
        i, *index = index
        if index == (None,):
            index = ()
        result = self.offset
        for k in self.site_types_to_sum:
            result += (
                self.coefs
                if (
                    type(self.coefs) is float
                    or type(self.coefs) is int
                    or type(self.coefs) is SimpleParam
                )
                else self.coefs[(i, k, *index)]
            ) * self.desc._pyomo_var[(i, k, *index)]
        return result


class SumBondTypes(Expr):
    """A class for expressions formed by summation over building block types.

    Attributes:
        desc (Descriptor/Expr): descriptors or expressions to sum over
        coefs (float/list<float>): coefficients to multiply contributions
            from each pair of building block types
        offset (float): coefficient to add to the expression
        bond_types_to_sum (list<tuple<BBlock,BBlock>>): building block pairs
            to consider in the summation
            (index information inherited from IndexedElem)
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, desc, coefs=1.0, offset=0.0, bond_types_to_sum=None, **kwargs):
        """Standard constructor for summation of contributions by bond-type.

        Args:
            desc (Descriptor): descriptors or expressions to sum across bond types
            coefs (float/list<float>): Optional, coefficients to multiple each
                bond-type term by. Default=1.0
            offset (float): Optional, coefficient to add to the expression.
                Default=0.0
            bonds_types_to_sum (list<tuple<BBlock,BBlock>>): Optional, subset
                of bond types to sum.
                Default=None, meaning all bond-types in the desc object are
                considered.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.desc = desc
        self.coefs = coefs
        self.offset = offset
        self.bond_types_to_sum = (
            bond_types_to_sum if bond_types_to_sum is not None else desc.bond_types
        )
        kwargs = {**desc.index_dict, **kwargs}
        kwargs.pop("bond_types")
        Expr.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_expr(self, index=None):
        """Interface for generating Pyomo expressions.

        Args:
        index (list): Optional, index to to create an instance of a Pyomo
            expression. In the case of a scalar, the valid index is None.

        Returns:
        An instance of a Pyomo expression.
        """
        assert index is not None
        i, j, *index = index
        if index == (None,):
            index = ()
        result = self.offset
        for k, l in self.bond_types_to_sum:
            result += (
                self.coefs
                if (
                    type(self.coefs) is float
                    or type(self.coefs) is int
                    or type(self.coefs) is SimpleParam
                )
                else self.coefs[(i, j, k, l, *index)]
            ) * self.desc._pyomo_var[(i, j, k, l, *index)]
        return result


class SumSitesAndTypes(Expr):
    """A class for expressions formed by summation over sites and building
    block types.

    Attributes:
        desc (Descriptor/Expr): descriptors or expressions to sum over
        coefs (float/list<float>): coefficients to multiply contributions
            from each building block type
        offset (float): coefficient to add to the expression
        sites_to_sum (list<int>): sites to consider in the summation
        site_types_to_sum (list<BBlock>): building block types to consider in
            the summation (index information inherited from IndexedElem)
    """

    # === STANDARD CONSTRUCTOR
    def __init__(
        self,
        desc,
        coefs=1.0,
        offset=0.0,
        sites_to_sum=None,
        site_types_to_sum=None,
        **kwargs
    ):
        """Standard constructor for summation of site contributions.

        Args:
            desc (Descriptor): descriptors or expressions to sum across all
                sites and site types
            coefs (float/list<float>): Optional, coefficients to multiple each
                site and site-type term by. Default=1.0
            offset (float): Optional, coefficient to add to the expression.
                Default=0.0
            sites_to_sum (list<int>): Optional, subset of canvas sites to sum.
                Default=None, meaning all sites in the desc object are considered.
            site_types_to_sum (list<BBlock>): Optional, subset of site types to
                sum. Default=None, meaning all site types in the desc object are
                considered.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.desc = desc
        self.coefs = coefs
        self.offset = offset
        self.sites_to_sum = sites_to_sum if sites_to_sum is not None else desc.sites
        self.site_types_to_sum = (
            site_types_to_sum if site_types_to_sum is not None else desc.site_types
        )
        Expr.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_expr(self, index=None):
        """Interface for generating Pyomo expressions.

        Args:
            index (list): Optional, index to to create an instance of a Pyomo
                expression. In the case of a scalar, the valid index is None.

        Returns:
            An instance of a Pyomo expression.
        """
        if index == (None,):
            index = ()
        result = self.offset
        for i in self.sites_to_sum:
            for k in self.site_types_to_sum:
                result += (
                    self.coefs
                    if (
                        type(self.coefs) is float
                        or type(self.coefs) is int
                        or type(self.coefs) is SimpleParam
                    )
                    else self.coefs[(i, k, *index)]
                ) * self.desc._pyomo_var[(i, k, *index)]
        return result


class SumBondsAndTypes(Expr):
    """A class for expressions formed by summation over bonds and bond types.

    Attributes:
        desc (Descriptor/Expr): descriptors or expressions to sum over
        coefs (float/list<float>): coefficients to multiply contributions
            from each combination of bond and building block type
        offset (float): coefficient to add to the expression
        bonds_to_sum (list<tuple<int,int>>): bonds to consider in the summation
        bond_types_to_sum (list<tuple<BBlock,BBlock>>): building block types to
            consider in the summation
            (index information inherited from IndexedElem)
    """

    # === STANDARD CONSTRUCTOR
    def __init__(
        self,
        desc,
        coefs=1.0,
        offset=0.0,
        bonds_to_sum=None,
        bond_types_to_sum=None,
        **kwargs
    ):
        """Standard constructor for summation of contributions by bond-type.

        Args:
            desc (Descriptor): descriptors or expressions to sum across bonds
                and bond types
            coefs (float/list<float>): Optional, coefficients to multiple each
                term by. Default=1.0
            offset (float): Optional, coefficient to add to the expression.
                Default=0.0
            bonds_to_sum (list<int>): Optional, subset of bonds to sum.
                Default=None, meaning all bonds in the desc object are considered.
            bonds_types_to_sum (list<int>): Optional, subset of bond types
                to sum. Default=None, meaning all bond-types in the desc object are
                considered.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.desc = desc
        self.coefs = coefs
        self.offset = offset
        self.bonds_to_sum = bonds_to_sum if bonds_to_sum is not None else desc.bonds
        self.bond_types_to_sum = (
            bond_types_to_sum if bond_types_to_sum is not None else desc.bond_types
        )
        Expr.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_expr(self, index=None):
        """Interface for generating Pyomo expressions.

        Args:
        index (list): Optional, index to to create an instance of a Pyomo
            expression. In the case of a scalar, the valid index is None.

        Returns:
        An instance of a Pyomo expression.
        """
        if index == (None,):
            index = ()
        result = self.offset
        for i, j in self.bonds_to_sum:
            for k, l in self.bond_types_to_sum:
                result += (
                    self.coefs
                    if (
                        type(self.coefs) is float
                        or type(self.coefs) is int
                        or type(self.coefs) is SimpleParam
                    )
                    else self.coefs[i, j, k, l]
                ) * self.desc._pyomo_var[i, j, k, l]
        return result


class SumConfs(Expr):
    """A class for expressions formed by summation over conformation types.

    Attributes:
        desc (Descriptor/Expr): descriptors or expressions to sum over
        coefs (float/list<float>): coefficients to multiply contributions
            from each building block type
        offset (float): coefficient to add to the expression
        confs_to_sum (list<int>): conformations to consider in the summation
            (index information inherited from IndexedElem)
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, Zic, coefs=1.0, offset=0.0, confs_to_sum=None, **kwargs):
        """Standard constructor for summation of bond contributions.

        Args:
            Zic (Descriptor): descriptors or expressions to sum across
                conformations
            coefs (float/list<float>): Optional, coefficients to multiple each
                conformation term by. Default=1.0
            offset (float): Optional, coefficient to add to the expression.
                Default=0.0
            confs_to_sum (list<int>): Optional, subset of conformations to sum
                Default=None, meaning all conformations in the Zic object are
                considered.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.Zic = Zic
        self.coefs = coefs
        self.offset = offset
        self.confs_to_sum = confs_to_sum if confs_to_sum is not None else Zic.confs
        kwargs = {**Zic.index_dict, **kwargs}
        kwargs.pop("confs")
        Expr.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_expr(self, index=None):
        """Interface for generating Pyomo expressions.

        Args:
            index (list): Optional, index to to create an instance of a Pyomo
                expression. In the case of a scalar, the valid index is None.

        Returns:
            An instance of a Pyomo expression.
        """
        i, *index = index
        if index == (None,):
            index = ()
        result = self.offset
        for c in self.confs_to_sum:
            result += (
                self.coefs
                if (
                    type(self.coefs) is float
                    or type(self.coefs) is int
                    or type(self.coefs) is SimpleParam
                )
                else self.coefs[(i, c, *index)]
            ) * self.Zic._pyomo_var[(i, c, *index)]
        return result


class SumSitesAndConfs(Expr):
    """A class for expressions formed by summation over sites and building
    block types.

    Attributes:
        desc (Descriptor/Expr): descriptors or expressions to sum over
        coefs (float/list<float>): coefficients to multiply contributions
            from each building block type
        offset (float): coefficient to add to the expression
        sites_to_sum (list<int>): sites to consider in the summation
        confs_to_sum (list<BBlock>): conformations to consider in the summation
            (index information inherited from IndexedElem)
    """

    # === STANDARD CONSTRUCTOR
    def __init__(
        self, Zic, coefs=1.0, offset=0.0, sites_to_sum=None, confs_to_sum=None, **kwargs
    ):
        """Standard constructor for summation of site and conformation
           contributions.

        Args:
            Zic (Descriptor): descriptors or expressions to sum across all sites
                and conformations.
            coefs (float/list<float>): Optional, coefficients to multiple each
                site and conformation term by. Default=1.0
            offset (float): Optional, coefficient to add to the expression.
                Default=0.0
            sites_to_sum (list<int>): Optional, subset of canvas sites to sum.
                Default=None, meaning all sites in the desc object are considered.
            confs_to_sum (list<int>): Optional, subset of conformations to sum.
                Default=None, meaning all conformations in the desc object are
                considered.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.Zic = Zic
        self.coefs = coefs
        self.offset = offset
        self.sites_to_sum = sites_to_sum if sites_to_sum is not None else Zic.sites
        self.confs_to_sum = confs_to_sum if confs_to_sum is not None else Zic.confs
        Expr.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_expr(self, index=None):
        """Interface for generating Pyomo expressions.

        Args:
            index (list): Optional, index to to create an instance of a Pyomo
                expression. In the case of a scalar, the valid index is None.

        Returns:
            An instance of a Pyomo expression.
        """
        if index == (None,):
            index = ()
        result = self.offset
        for i in self.sites_to_sum:
            for c in self.confs_to_sum:
                result += (
                    self.coefs
                    if (
                        type(self.coefs) is float
                        or type(self.coefs) is int
                        or type(self.coefs) is SimpleParam
                    )
                    else self.coefs[(i, c, *index)]
                ) * self.Zic._pyomo_var[(i, c, *index)]
        return result


class DescriptorRule(IndexedElem):
    """An abstract base class for rules to define material descriptors.

    This class is only the abstract interface for other well-defined rules.
    Rules get attached to MaterialDescriptor objects and are intended to be
    interpretable in relation to the descriptor that they are attached to.

    Examples:
        'Coordination number is equal to the sum of variables for the presence
            of any bond to the neighbors of a site.'
            ->  m.CNi.rules.append(EqualTo(SumNeighborBonds(m.Bondij)))
       'Size of a cluster is equal to the number of of atoms'
           ->  m.ClusterSize.rules.append(EqualTo(SumSites(m.Yi)))

    Rules can be applied over a subset of the design space (i.e., only for
    some combinations of sites, site-types, etc.) by providing keywords
    that are passed to the constructor of IndexedElem. Else, the relevant
    indexes are infered from the rule components.

    Attributes:
        (index information inherited from IndexedElem)
    """

    def __init__(self, **kwargs):
        """Standard constructor for DescriptorRule base class.

        Args:
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        IndexedElem.__init__(self, **kwargs)

    @abstractmethod
    def _pyomo_cons(self, var):
        """Abstract method to define necessary interface for descriptor rules.

        Args:
        var (MaterialDescriptor): Variable that the rule is attached to.
            (i.e., the variable that should be read before the class name to
            interpret the rule)

        Returns:
        (list<Constraint>) list of Pyomo constraint objects.
        """
        raise NotImplementedError


class SimpleDescriptorRule(DescriptorRule):
    """An base class for simple rules with a left and right hand side.

    This class is just intended to create a common interface for the
    EqualTo, LessThan, and GreaterThan rules.

    Attributes:
        expr (Expr): Right-hand side expression for the rule.
            (index information inherited from IndexedElem)
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, e, **kwargs):
        """Standard constructor for simple descriptor rules.

        Args:
            e (float/int/Expr): An expression to use are right hand side of
                a rule. If a float or int, a LinearExpr is created for the user.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.expr = LinearExpr(offset=e) if type(e) is float or type(e) is int else e
        kwargs = {**self.expr.index_dict, **kwargs}
        DescriptorRule.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_cons(self, var):
        """Method to create a Pyomo constraint from this rule.

        Args:
        var (MaterialDescriptor): The descriptor to be defined by this rule.

        Returns:
        (list<Constraint>) list of Pyomo constraint objects.
        """
        ConIndexes = IndexedElem.fromComb(var, self)
        return [Constraint(*ConIndexes.index_sets, rule=self._pyomo_rule(var))]

    def _pyomo_rule(self, LHS, operator, RHS):
        """Method to create a function for a Pyomo constraint rule.

        Args:
            LHS (MaterialDescriptor/Expr): The left hand side of a simple rule.
            operator (function): The relationship to encode in a rule.
            RHS (MaterialDescriptor/Expr): The right hand side of a simple rule.

        Returns:
            (function) A function interpretable by Pyomo for a 'rule' argument
        """
        ConIndexes = IndexedElem.fromComb(LHS, RHS)

        def rule(m, *args):
            LHS_index = LHS.mask(args, ConIndexes)
            RHS_index = RHS.mask(args, ConIndexes)
            return operator(LHS._pyomo_expr(LHS_index), RHS._pyomo_expr(RHS_index))

        return rule


class LessThan(SimpleDescriptorRule):
    """A class for rules implementing 'less than or equal to' an expression.

    Spelled out: 'the descriptor is less than or equal to a linear expression'

    Attributes:
        expr (Expr): Right-hand side expression for the rule.
            (index information inherited from IndexedElem)

    See DescriptorRule for more information.
    """

    # === STANDARD CONSTRUCTOR
    # --- Inherited from SimpleDescriptorRule ---

    # === PROPERTY EVALUATION METHODS
    def _pyomo_rule(self, desc):
        """Method to create a function for a Pyomo constraint rule.

        Args:
            desc (MaterialDescriptor/Expr): A descriptor to define as 'less than'
                the expression for this rule.

        Returns:
            (function) A function in the format of a Pyomo rule to construct a
            constraint.
        """

        def less_than(LHS, RHS):
            return LHS <= RHS

        return SimpleDescriptorRule._pyomo_rule(self, desc, less_than, self.expr)


class EqualTo(SimpleDescriptorRule):
    """A class for rules implementing 'equal to' an expression.

    Spelled out: 'the descriptor is equal to a linear expression'

    Attributes:
        expr (Expr): Right-hand side expression for the rule.
            (index information inherited from IndexedElem)

    See DescriptorRule for more information.
    """

    # === STANDARD CONSTRUCTOR
    # --- Inherited from SimpleDescriptorRule ---

    # === PROPERTY EVALUATION METHODS
    def _pyomo_rule(self, desc):
        """Method to create a function for a Pyomo constraint rule.

        Args:
        desc (MaterialDescriptor/Expr): A descriptor to define as 'equal to'
            the expression for this rule.

        Returns:
        (function) A function in the format of a Pyomo rule to construct a
            constraint.
        """

        def equal_to(LHS, RHS):
            return LHS == RHS

        return SimpleDescriptorRule._pyomo_rule(self, desc, equal_to, self.expr)


class GreaterThan(SimpleDescriptorRule):
    """A class for rules implementing 'greater than or equal to' an expr.

    Spelled out: 'descriptor is greater than or equal to a linear expression'

    Attributes:
        expr (Expr): Right-hand side expression for the rule.
            (index information inherited from IndexedElem)

    See DescriptorRule for more information.
    """

    # === STANDARD CONSTRUCTOR
    # --- Inherited from SimpleDescriptorRule ---

    # === PROPERTY EVALUATION METHODS
    def _pyomo_rule(self, desc):
        """Method to create a function for a Pyomo constraint rule.

        Args:
            desc (MaterialDescriptor/Expr): A descriptor to define as 'greater
                than' the expression for this rule.

        Returns:
            (function) A function in the format of a Pyomo rule to construct a
            constraint.
        """

        def greater_than(LHS, RHS):
            return LHS >= RHS

        return SimpleDescriptorRule._pyomo_rule(self, desc, greater_than, self.expr)


class FixedTo(DescriptorRule):
    """A class for rules that fix descriptors to required values.

    Spelled out: 'the descriptor is fixed to a scalar value'

    Attributes:
        val (float): the value that the descriptor is fixed to.
            (index information inherited from IndexedElem)

    See DescriptorRule for more information.
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, val, **kwargs):
        """Standard constructor for FixedTo rules.

        Args:
            val (float): The value to fix descriptors to.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.val = val
        DescriptorRule.__init__(self, **kwargs)

    def _pyomo_cons(self, var):
        """Method to create a Pyomo constraint from this rule.

        Args:
        var (MaterialDescriptor): The descriptor to be defined by this rule.

        Returns:
        (list<Constraint>) list of Pyomo constraint objects.
        """
        # NOTE: This method is used to ensure that basic variables that
        #       are fixed get referenced to write basic constraints in
        #       the model. Don't make constraints, but do instantiate
        #       variables
        Comb = IndexedElem.fromComb(var, self)
        for k in Comb.keys():
            var._pyomo_var[k]
        return []


class Disallow(DescriptorRule):
    """A class for rules that disallow a previously-identified design.

    Spelled out: 'the descriptors must attain a different solution than
        a given design'

    Attributes:
        D (Design): the design from which variable values to disallow are infered
            (index information inherited from IndexedElem)

    See DescriptorRule for more information.
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, D, **kwargs):
        """Standard constructor for the Disallow rule.

        Args:
            D (Design): A design object to make infeasible in the resulting model
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.D = D
        DescriptorRule.__init__(self, **kwargs)

    def _pyomo_expr(self, var):
        """Method to create the integer cut for this disallowed design.

        Args:
            var (MaterialDescriptor): The descriptor to be defined by this rule.

        Returns:
            An instance of a Pyomo expression.
        """
        if var.name == "Yi":
            result = 0
            for i in range(len(self.D.Canvas)):
                if self.D.Contents[i] is None:
                    result += var._pyomo_var[i]
                else:
                    result += 1 - var._pyomo_var[i]
        elif var.name == "Yik":
            result = 0
            for i in range(len(self.D.Canvas)):
                for k in self.D.NonVoidElems:
                    if self.D.Contents[i] is not k:
                        result += var._pyomo_var[i, k]
                    else:
                        result += 1 - var._pyomo_var[i, k]
        else:
            # NOTE: This rule was intended to disallow structures
            #       or labelings of structures (i.e., Yi or Yik
            #       variables). It is not clear how we can generally
            #       disallow any general combination of variables.
            raise ValueError("Decide what to do in this case...")
        return result

    def _pyomo_cons(self, var):
        """
        Method to create a Pyomo constraint from this rule.

        Args:
            var (MaterialDescriptor): The descriptor to be defined by this rule.

        Returns:
            (list<Constraint>) list of Pyomo constraint objects.
        """
        return Constraint(expr=(self._pyomo_expr(var) >= 1))


class PiecewiseLinear(DescriptorRule):
    """A class for rules implementing 'equal to a piecewise linear function'.

    Spelled out: 'the descriptor is equal to a piecewise linear expression'

    Note: Innequalities of 'less than' or 'greater than' a piecewie function
    can be achieved by introducing an auxiliary descriptor to be equal to the
    piecewise function. Then, inequalities can be introduced using the
    auxiliary descriptor. Alternatively, users can modify the con_type
    attribute that is interpreted by Pyomo.

    Attributes:
        values (list<float>): values of univariate piecewise linear function at
            each breakpoint.
        breakpoints (list<float>): breakpoints of the piecewise linear function.
        input_desc (MaterialDescriptor): descriptor as the arugment to the
            piecewise linear function
        con_type (string): indicates the bound type of the piecewise function.
            Options:
            “UB” - relevant descriptor is bounded above by piecewise function.
            “LB” - relevant descriptor is bounded below by piecewise function.
            “EQ” - relevant descriptor is equal to the piecewise function. (Default)

    See DescriptorRule for more information.
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, values, breakpoints, input_desc, con_type="EQ", **kwargs):
        """Standard constructor for simple descriptor rules.

        Args:
            values (list<float>): values of the function.
            breakpoints (list<float>): breakpoints of the function.
            input_desc (MaterialDescriptor): arugment to the function
            con_type (string): Optional, indicates the bound type of the
                piecewise function
                Options:
                    “UB” - bounded above by piecewise function.
                    “LB” - bounded below by piecewise function.
                    “EQ” - equal to the piecewise function. (Default)
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.values = values
        self.breakpoints = breakpoints
        self.input_desc = input_desc
        self.con_type = con_type.upper()
        DescriptorRule.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_cons(self, var):
        """Method to create a Pyomo constraint from this rule.

        Args:
            var (MaterialDescriptor): The descriptor to be defined by this rule.

        Returns:
            (list<Block>) list of Pyomo model block objects created by Piecewise
                function.
        """
        Comb = IndexedElem.fromComb(var, self)
        return [
            Piecewise(
                *Comb.index_sets,
                var._pyomo_var,
                self.input_desc._pyomo_var,
                pw_pts=self.breakpoints,
                f_rule=self.values,
                pw_constr_type=self.con_type,
                pw_repn="MC"
            )
        ]


class Implies(DescriptorRule):
    """A class for rules that define simple logical implications.

    Spelled out: 'if this descriptor is true (i.e., equal to one),
        then another set of simple rules also apply'

    Attributes:
        concs (list<tuple<MaterialDescriptor,SimpleDescriptorRule>>):
            list of conclusions to enforce if the logical predicate is true.
            (index information inherited from IndexedElem)

    See DescriptorRule for more information.
    """

    DEFAULT_BIG_M = 9999

    # === STANDARD CONSTRUCTOR
    def __init__(self, concs, **kwargs):
        """Standard constructor for Implies rule.

        Args:
            concs (list<tuple<MaterialDescriptor,SimpleDescriptorRule>>):
                list of conclusions to conditionally enforce. Also, a single
                conclusion can be provided (i.e., a tuple<MaterialDescriptor,
                SimpleDescriptorRule>) and will be placed in a list.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.concs = concs if type(concs) is list else [concs]
        Comb = IndexedElem.fromComb(
            *(desc for desc, conc in self.concs), *(conc for desc, conc in self.concs)
        )
        kwargs = {**Comb.index_dict, **kwargs}
        DescriptorRule.__init__(self, **kwargs)

    def _pyomo_cons(self, var):
        """Method to create a Pyomo constraint from this rule.

        Args:
            var (MaterialDescriptor): The descriptor to be defined by this rule.

        Returns:
            (list<Constraint>) list of Pyomo constraint objects.
        """

        result = []
        for (desc, conc) in self.concs:
            ConIndexes = IndexedElem.fromComb(var, desc, conc)

            def rule_lb(m, *args):
                v = var._pyomo_expr(index=var.mask(args, ConIndexes))
                d = desc._pyomo_expr(index=desc.mask(args, ConIndexes))
                c = conc.expr._pyomo_expr(index=conc.expr.mask(args, ConIndexes))
                body_lb = getLB(d - c)
                MLB = body_lb if body_lb is not None else -Implies.DEFAULT_BIG_M
                return MLB * (1 - v) <= d - c

            def rule_ub(m, *args):
                v = var._pyomo_expr(index=var.mask(args, ConIndexes))
                d = desc._pyomo_expr(index=desc.mask(args, ConIndexes))
                c = conc.expr._pyomo_expr(index=conc.expr.mask(args, ConIndexes))
                body_ub = getUB(d - c)
                MUB = body_ub if body_ub is not None else Implies.DEFAULT_BIG_M
                return d - c <= MUB * (1 - v)

            if isinstance(conc, LessThan) or isinstance(conc, EqualTo):
                result.append(Constraint(*ConIndexes.index_sets, rule=rule_ub))
            if isinstance(conc, GreaterThan) or isinstance(conc, EqualTo):
                result.append(Constraint(*ConIndexes.index_sets, rule=rule_lb))
        return result


class NegImplies(DescriptorRule):
    """A class for rules that define logical implications with negation.

    Spelled out: 'if this descriptor is not true (i.e., is equal to zero),
        then another simple rule also applies'

    Attributes:
        concs (list<tuple<MaterialDescriptor,SimpleDescriptorRule>>):
            list of conclusions to enforce if the logical predicate is false.
            (index information inherited from IndexedElem)

    See DescriptorRule for more information.
    """

    DEFAULT_BIG_M = 9999

    # === STANDARD CONSTRUCTOR
    def __init__(self, concs, **kwargs):
        """Standard constructor for NegImplies rule.

        Args:
            concs (list<tuple<MaterialDescriptor,SimpleDescriptorRule>>):
                list of conclusions to conditionally enforce. Also, a single
                conclusion can be provided (i.e., a tuple<MaterialDescriptor,
                SimpleDescriptorRule>) and will be placed in a list.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.concs = concs if type(concs) is list else [concs]
        Comb = IndexedElem.fromComb(
            *(desc for desc, conc in self.concs), *(conc for desc, conc in self.concs)
        )
        kwargs = {**Comb.index_dict, **kwargs}
        DescriptorRule.__init__(self, **kwargs)

    def _pyomo_cons(self, var):
        """Method to create a Pyomo constraint from this rule.

        Args:
            var (MaterialDescriptor): The descriptor to be defined by this rule.

        Returns:
            (list<Constraint>) list of Pyomo constraint objects.
        """
        result = []
        for (desc, conc) in self.concs:
            ConIndexes = IndexedElem.fromComb(var, desc, conc)

            def rule_lb(m, *args):
                v = var._pyomo_expr(index=var.mask(args, ConIndexes))
                d = desc._pyomo_expr(index=desc.mask(args, ConIndexes))
                c = conc.expr._pyomo_expr(index=conc.expr.mask(args, ConIndexes))
                body_lb = getLB(d - c)
                MLB = body_lb if body_lb is not None else -NegImplies.DEFAULT_BIG_M
                return MLB * (v) <= d - c

            def rule_ub(m, *args):
                v = var._pyomo_expr(index=var.mask(args, ConIndexes))
                d = desc._pyomo_expr(index=desc.mask(args, ConIndexes))
                c = conc.expr._pyomo_expr(index=conc.expr.mask(args, ConIndexes))
                body_ub = getUB(d - c)
                MUB = body_ub if body_ub is not None else NegImplies.DEFAULT_BIG_M
                return d - c <= MUB * v

            if isinstance(conc, LessThan) or isinstance(conc, EqualTo):
                result.append(Constraint(*ConIndexes.index_sets, rule=rule_ub))
            if isinstance(conc, GreaterThan) or isinstance(conc, EqualTo):
                result.append(Constraint(*ConIndexes.index_sets, rule=rule_lb))
        return result


class ImpliesSiteCombination(DescriptorRule):
    """A class for rules that define logical implications between two sites.

    Spelled out: 'if this bond-indexed descriptor is true (i.e., is equal
        to one), then a pair of simple rules hold on the two bonding sites'

    Attributes:
        canv (Canvas): the data structure to identify neighbor connections
            to apply rules over.
        concis (list<tuple<MaterialDescriptor,SimpleDescriptorRule>>):
            list of conclusions to enforce at the first site in the pair if
            the logical predicate is true.
        concjs (list<tuple<MaterialDescriptor,SimpleDescriptorRule>>):
            list of conclusions to enforce at the second site in the pair if
            the logical predicate is true.
        symmetric_bonds (bool): flag to indicate if implications should be
            applied over symmetric bond indices
            (index information inherited from IndexedElem)

    See DescriptorRule for more information.
    """

    DEFAULT_BIG_M = 9999

    # === STANDARD CONSTRUCTOR
    def __init__(self, canv, concis, concjs, symmetric_bonds=False, **kwargs):
        """Standard constructor for ImpliesSiteCombination rules.

        Args:
            canv (Canvas): the data structure to identify neighbor connections
                to apply rules over.
            concis (list<tuple<MaterialDescriptor,SimpleDescriptorRule>>):
                list of conclusions to conditionally enforce at the first
                site in a bond.
                Note: single conclusions can be provided (i.e., a
                tuple<MaterialDescriptor,SimpleDescriptorRule>) and will
                be placed in lists.
            concjs (list<tuple<MaterialDescriptor,SimpleDescriptorRule>>):
                list of conclusions to conditionally enforce at the second
                site in a bond.
                Note: single conclusions can be provided (i.e., a
                tuple<MaterialDescriptor,SimpleDescriptorRule>) and will
                be placed in lists.
            symmetric_bonds (bool): flag to indicate if a symmetric verions
                of bonds should be enumerated or if both directions should
                be included.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """

        self.canv = canv
        self.concis = concis if type(concis) is list else [concis]
        self.concjs = concjs if type(concjs) is list else [concjs]
        self.symmetric_bonds = symmetric_bonds
        Combi = IndexedElem.fromComb(
            *(desc for desc, conc in self.concis), *(conc for desc, conc in self.concis)
        )
        Combj = IndexedElem.fromComb(
            *(desc for desc, conc in self.concjs), *(conc for desc, conc in self.concjs)
        )
        assert Combi.sites is not None and Combj.sites is not None
        if "bonds" not in kwargs:
            kwargs["bonds"] = [
                (i, j)
                for i in Combi.sites
                for j in canv.NeighborhoodIndexes[i]
                if (
                    j is not None
                    and j in Combj.sites
                    and (not symmetric_bonds or j > i)
                )
            ]
        if sum(Combi.dims) > 1:
            raise NotImplementedError(
                "Additional indexes are not supported, please contact MatOpt developer for "
                "possible feature addition"
            )
        DescriptorRule.__init__(self, **kwargs)

    def _pyomo_cons(self, var):
        """Method to create a Pyomo constraint from this rule.

        Args:
            var (MaterialDescriptor): The descriptor to be defined by this rule.

        Returns:
            (list<Constraint>) list of Pyomo constraint objects.
        """
        assert var.bonds is not None
        Comb = IndexedElem.fromComb(var, self)
        result = []
        # NOTE: After much confusion, I found a bug in the line of code
        #       below. Be careful not to use variable names "expr"
        #       because it gets mixed up with the Pyomo module "expr".
        #       No error, but it gives garbage expressions and wasn't
        #       clear to me what was being generated...
        # NOTE: Not clear if this was caused by module "expr" or the
        #       conflict of two local "expr,conc" objects in the two
        #       for loops...
        # for expr,conc in self.concis:
        for expri, conci in self.concis:

            def rule_i_lb(m, i, j):
                e = expri._pyomo_expr(index=(i,))
                c = conci.expr._pyomo_expr(index=(i,))
                body = e - c
                body_LB = getLB(body)
                MLBi = (
                    body_LB
                    if body_LB is not None
                    else -ImpliesSiteCombination.DEFAULT_BIG_M
                )
                return MLBi * (1 - var._pyomo_var[i, j]) <= body

            def rule_i_ub(m, i, j):
                e = expri._pyomo_expr(index=(i,))
                c = conci.expr._pyomo_expr(index=(i,))
                body = e - c
                body_UB = getUB(body)
                MUBi = (
                    body_UB
                    if body_UB is not None
                    else ImpliesSiteCombination.DEFAULT_BIG_M
                )
                return body <= MUBi * (1 - var._pyomo_var[i, j])

            if isinstance(conci, GreaterThan) or isinstance(conci, EqualTo):
                result.append(Constraint(*Comb.index_sets, rule=rule_i_lb))
            if isinstance(conci, LessThan) or isinstance(conci, EqualTo):
                result.append(Constraint(*Comb.index_sets, rule=rule_i_ub))
        # NOTE: See note above for variable name "expr"
        # for expr,conc in self.concjs:
        for exprj, concj in self.concjs:

            def rule_j_lb(m, i, j):
                e = exprj._pyomo_expr(index=(j,))
                c = concj.expr._pyomo_expr(index=(j,))
                body = e - c
                body_LB = getLB(body)
                MLBj = (
                    body_LB
                    if body_LB is not None
                    else -ImpliesSiteCombination.DEFAULT_BIG_M
                )
                return MLBj * (1 - var._pyomo_var[i, j]) <= body

            def rule_j_ub(m, i, j):
                e = exprj._pyomo_expr(index=(j,))
                c = concj.expr._pyomo_expr(index=(j,))
                body = e - c
                body_UB = getUB(body)
                MUBj = (
                    body_UB
                    if body_UB is not None
                    else ImpliesSiteCombination.DEFAULT_BIG_M
                )
                return body <= MUBj * (1 - var._pyomo_var[i, j])

            if isinstance(concj, GreaterThan) or isinstance(concj, EqualTo):
                result.append(Constraint(*Comb.index_sets, rule=rule_j_lb))
            if isinstance(concj, LessThan) or isinstance(concj, EqualTo):
                result.append(Constraint(*Comb.index_sets, rule=rule_j_ub))
        return result


class ImpliesNeighbors(DescriptorRule):
    """A class for rules that define logical implications on neighbor sites.

    Spelled out: 'if this site-indexed descriptor is true (i.e., is equal
        to one), then a set of simple rules hold on each of the neighboring
        sites'

    Attributes:
        concs (list<tuple<MaterialDescriptor,SimpleDescriptorRule>>):
            list of conclusions to enforce if the logical predicate is true.
        neighborhoods (list<list<int>>): neighborhood data structure to use
            if you do not want to use the neighborhoods of the descriptor
            that this rule is attached to.
            (index information inherited from IndexedElem)

    See DescriptorRule for more information on rules and Canvas for more
    information on 'neighborhoods'.
    """

    DEFAULT_BIG_M = 9999

    # === STANDARD CONSTRUCTOR
    def __init__(self, concs, neighborhoods=None, **kwargs):
        """Standard constructor for ImpliesNeighbors rules.

        Args:
            concs (list<tuple<MaterialDescriptor,SimpleDescriptorRule>>):
                list of conclusions to conditionally enforce. Also, a single
                conclusion can be provided (i.e., a tuple<MaterialDescriptor,
                SimpleDescriptorRule>) and will be placed in a list.
            neighborhoods (list<list<int>>) Optional, data structure to use
                as neighborhoods of interest. If not provided, then the
                neighborhoods of the descriptor that this rule is attached to
                is used.
            **kwargs: Optional, index information passed to IndexedElem if
                interested in a subset of indices
                Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.concs = concs if type(concs) is list else [concs]
        self.neighborhoods = neighborhoods
        Comb = IndexedElem.fromComb(
            *(desc for desc, conc in self.concs), *(conc for desc, conc in self.concs)
        )
        assert Comb.sites is not None
        kwargs = {**Comb.index_dict, **kwargs}
        DescriptorRule.__init__(self, **kwargs)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_cons(self, var):
        """Method to create a Pyomo constraint from this rule.

        Args:
            var (MaterialDescriptor): The descriptor to be defined by this rule.

        Returns:
            (list<Constraint>) list of Pyomo constraint objects.
        """
        var_dict_wo_s = var.index_dict
        var_dict_wo_s.pop("sites")  # no need to capture these sites
        neighborhoods = (
            self.neighborhoods
            if self.neighborhoods is not None
            else var.canv.NeighborhoodIndexes
        )
        bonds = [(i, j) for i in var.sites for j in neighborhoods[i] if j is not None]
        result = []
        # NOTE: After much confusion, I found a bug in the line of code
        #       below. Be careful not to use variable names "expr"
        #       because it gets mixed up with the Pyomo module "expr".
        #       No error, but it gives garbage expressions and wasn't
        #       clear to me what was being generated...
        # for expr,conc in self.concs:
        for expr_, conc in self.concs:
            Comb = IndexedElem.fromComb(expr_, conc)
            r_dict_wo_s = Comb.index_dict
            r_dict_wo_s.pop("sites")  # no need to capture these sites
            ConIndexes = IndexedElem.fromComb(
                IndexedElem(bonds=bonds),
                IndexedElem(**var_dict_wo_s),
                IndexedElem(**r_dict_wo_s),
            )

            def rule_lb(m, *args):
                i, j, *args = args
                v = var._pyomo_var[var.mask((i, None, *args), ConIndexes)]
                e = expr_._pyomo_expr(index=expr_.mask((j, None, *args), ConIndexes))
                c = conc.expr._pyomo_expr(
                    index=conc.expr.mask((j, None, *args), ConIndexes)
                )
                body = e - c
                body_LB = getLB(body)
                MLB = (
                    body_LB if body_LB is not None else -ImpliesNeighbors.DEFAULT_BIG_M
                )
                return MLB * (1 - v) <= body

            def rule_ub(m, *args):
                i, j, *args = args
                v = var._pyomo_var[var.mask((i, None, *args), ConIndexes)]
                e = expr_._pyomo_expr(index=expr_.mask((j, None, *args), ConIndexes))
                c = conc.expr._pyomo_expr(
                    index=conc.expr.mask((j, None, *args), ConIndexes)
                )
                body = e - c
                body_UB = getUB(body)
                MUB = body_UB if body_UB is not None else ImpliesNeighbors.DEFAULT_BIG_M
                return body <= MUB * (1 - v)

            if isinstance(conc, GreaterThan) or isinstance(conc, EqualTo):
                result.append(Constraint(*ConIndexes.index_sets, rule=rule_lb))
            if isinstance(conc, LessThan) or isinstance(conc, EqualTo):
                result.append(Constraint(*ConIndexes.index_sets, rule=rule_ub))
        return result


class MaterialDescriptor(IndexedElem):
    """A class to represent material geometric and energetic descriptors.

    This class holds the information to define mathematical optimization
    variables for the properties of materials. Additionally, each descriptor
    has a 'rules' list to which the user can append rules defining the
    descriptor and constraining the design space.

    Attributes:
        name (string): A unique (otherwise Pyomo will complain) name
        canv (``Canvas``): The canvas that the descriptor will be indexed over
        atoms (list<``BBlock``>): The building blocks to index the descriptor over.
        confDs (list<``Design``>): The designs for conformations to index over.
        integer (bool): Flag to indicate if the descriptor takes integer values.
        binary (bool): Flag to indicate if the descriptor takes boolean values.
        rules (list<``DescriptorRules``>): List of rules to define and constrain
            the material descriptor design space.
        bounds (tuple/dict/func): If tuple, the lower and upper bounds on the
            descriptor values across all indices. If dict, the bounds can be
            individually set for each index.

    See ``IndexedElem`` for more information on indexing.
    See ``DescriptorRule`` for information on defining descriptors.
    """

    DBL_TOL = 1e-5

    # === STANDARD CONSTRUCTOR
    def __init__(
        self,
        name,
        canv=None,
        atoms=None,
        confDs=None,
        bounds=(None, None),
        integer=False,
        binary=False,
        rules=[],
        **kwargs
    ):
        """Standard constuctor for material descriptors.

        Note: It is generally not necessary for users to create
              MaterialDescriptors themselves. Instead, use the
              MatOptModel.add____Descriptor() methods for the right
              type of descriptor (i.e., Site, Bond, etc.).

        Args:
        name (string): A unique (otherwise Pyomo will complain) name
        canv (Canvas): The canvas that the descriptor will be indexed over
        atoms (list<BBlock>): Building blocks to index the descriptor over.
        confDs (list<Design>): The designs for conformations to index over.
        bounds (tuple/dict/func): If tuple, the lower and upper bounds on the
            descriptor values across all indices. If dict, the bounds can be
            individually set for each index. Otherwise, advanced users can
            specify a function to be interpreted by Pyomo.
        integer (bool): Flag to indicate if the descriptor is integer.
        binary (bool): Flag to indicate if the descriptor is boolean.
        rules (list<DescriptorRules>): List of rules to define and constrain
            the material descriptor design space.
        **kwargs: Optional, index information passed to IndexedElem if
            interested in a subset of indices.
            Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self._name = name
        self._canv = canv
        self._atoms = atoms
        self._confDs = confDs
        self._integer = integer or binary
        self._binary = binary
        self._rules = rules if type(rules) is list else [rules]
        self._bounds = bounds
        self._pyomo_var = None  # Will be set by MatOptModel._make_pyomo_model
        IndexedElem.__init__(self, **kwargs)

    # === AUXILIARY METHODS
    def _fix_pyomo_var_by_rule(self, r, m):
        if self.name in ("Yik", "Yi", "Xijkl", "Xij", "Cikl", "Ci", "Zic"):
            return self.__fix_basic_pyomo_vars_by_rule(r, m)
        else:
            Comb = IndexedElem.fromComb(self, r)
            for k in Comb.keys():
                self._pyomo_var[k].fix(r.val)

    def __fix_basic_pyomo_vars_by_rule(self, r, m):
        Comb = IndexedElem.fromComb(self, r)
        if self.name == "Yik":
            for i in Comb.sites:
                for k in Comb.site_types:
                    fixYik(m, i, k, r.val)
        elif self.name == "Yi":
            for i in Comb.sites:
                fixYi(m, i, r.val)
        elif self.name == "Xijkl":
            for i, j in Comb.bonds:
                for k, l in Comb.bond_types:
                    fixXijkl(m, i, j, k, l, r.val)
        elif self.name == "Xij":
            for i, j in Comb.bonds:
                fixXij(m, i, j, r.val)
        elif self.name == "Cikl":
            for i in Comb.sites:
                for k, l in Comb.bond_types:
                    fixCikl(m, i, k, l, r.val)
        elif self.name == "Ci":
            for i in Comb.sites:
                fixCi(m, i, r.val)
        elif self.name == "Zic":
            for i in Comb.sites:
                for c in Comb.confs:
                    fixZic(m, i, c, r.val)

    # === PROPERTY EVALUATION METHODS
    def _pyomo_cons(self, m):
        """Create a list of Pyomo constraints related to this descriptor."""
        result = []
        for rule in self.rules:
            if rule is not None:
                result.extend(rule._pyomo_cons(self))
        return result

    @property
    def _pyomo_bounds(self):
        """Creates a bound rule/tuple that can interpreted by Pyomo."""
        if type(self.bounds) is tuple:
            return self.bounds
        elif type(self.bounds) is dict:

            def rule_gen(m, *args):
                if args is not None and len(args) == 1:
                    args = args[0]
                return self.bounds[args]

            return rule_gen
        else:
            # Else, assume that the user knows what they're doing
            # with functions for pyomo bounds
            return self.bounds

    def _pyomo_expr(self, index=None):
        """Interprets a variable as a Pyomo expression.

        Note: This is just necessary so that we can conveniently interpret
            MaterialDescriptor objects in place of Expr objects.
        """
        return self._pyomo_var[index]

    @property
    def values(self):
        """Creates a dictionary of desriptor values after optimization.

        Note: Uses the Pyomo 'value' function and only works after the
            optimization of a model.

        Returns:
            (dict) Dictionary of keys to values after optimization.
        """
        return {index: value(self._pyomo_var[index]) for index in self.keys()}

    # === BASIC QUERY METHODS
    @property
    def name(self):
        return self._name

    @property
    def canv(self):
        return self._canv

    @property
    def atoms(self):
        return self._atoms

    @property
    def confDs(self):
        return self._confDs

    @property
    def bounds(self):
        return self._bounds

    @property
    def integer(self):
        return self._integer

    @property
    def binary(self):
        return self._binary

    @property
    def continuous(self):
        return not self.integer

    @property
    def rules(self):
        return self._rules


class MatOptModel(object):
    """A class for the specification of a materials optimization problem.

    Once all the material information is specified, we use this class to
    specify the material design problem of interest. This class is intended
    to be interpretable without mathematical optimization background while
    the conversion to Pyomo optimization models happens automatically.

    Attributes:
        canv (``Canvas``): The canvas of the material design space
        atoms (list<``BBlock``>): The list of building blocks to consider.
            Note: This list does not need to include a void-atom type. We use
            'None' to represent the absence of any building block at a given site.
        confDs (list<``Design``>): The list of conformations to consider.
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, canv, atoms=None, confDs=None):
        """Standard constructor for materials optimization problems.

        Args:
        canv (``Canvas``): The canvas of the material design space
        atoms (list<``BBlock``>): The list of building blocks to consider.
            Note: This list does not need to include a void-atom type. We use
            'None' to represent the absence of any building block at a given site.
        confDs (list<``Design``>): The list of conformations to consider.
        """
        self._canv = canv
        self._atoms = atoms
        self._confDs = confDs
        self._descriptors = []
        self.addSitesDescriptor("Yi", binary=True, rules=None)
        self.addBondsDescriptor("Xij", binary=True, rules=None)
        self.addNeighborsDescriptor("Ci", integer=True, rules=None)
        self.addSitesTypesDescriptor("Yik", binary=True, rules=None)
        self.addBondsTypesDescriptor("Xijkl", binary=True, rules=None)
        self.addNeighborsTypesDescriptor("Cikl", integer=True, rules=None)
        self.addSitesConfsDescriptor("Zic", binary=True, rules=None)

    # === MANIPULATION METHODS
    def addGlobalDescriptor(
        self, name, bounds=(None, None), integer=False, binary=False, rules=None
    ):
        """
        Method to add scalar descriptor to the model.

        Args:
            name (string): A unique (otherwise Pyomo will complain) name.
            bounds (tuple/dict/func): If tuple, the lower and upper bounds on the
                descriptor values across all indices. If dict, the bounds can be
                individually set for each index. Otherwise, advanced users can
                specify a function to be interpreted by Pyomo.
            integer (bool): Flag to indicate if the descriptor is integer.
            binary (bool): Flag to indicate if the descriptor is boolean.
            rules (list<DescriptorRules>): List of rules to define and constrain
                the material descriptor design space.
        """
        assert not hasattr(self, name)
        Desc = MaterialDescriptor(
            name=name, bounds=bounds, integer=integer, binary=binary, rules=rules
        )
        setattr(self, name, Desc)
        self._descriptors.append(Desc)

    def addSitesDescriptor(
        self,
        name,
        sites=None,
        bounds=(None, None),
        integer=False,
        binary=False,
        rules=None,
    ):
        """Method to add a site-indexed descriptor to the model.

        Args:
            name (string): A unique (otherwise Pyomo will complain) name.
            sites (list<int>): Optional, subset of canvas sites to index the new
                descriptor over.
                Default: None (i.e., all sites in canvas are considered)
            bounds (tuple/dict/func): If tuple, the lower and upper bounds on the
                descriptor values across all indices. If dict, the bounds can be
                individually set for each index. Otherwise, advanced users can
                specify a function to be interpreted by Pyomo.
            integer (bool): Flag to indicate if the descriptor is integer.
            binary (bool): Flag to indicate if the descriptor is boolean.
            rules (list<DescriptorRules>): List of rules to define and constrain
                the material descriptor design space.
        """
        assert not hasattr(self, name)
        sites = sites if sites is not None else list(range(len(self.canv)))
        Desc = MaterialDescriptor(
            name=name,
            canv=self.canv,
            sites=sites,
            bounds=bounds,
            integer=integer,
            binary=binary,
            rules=rules,
        )
        setattr(self, name, Desc)
        self._descriptors.append(Desc)

    def addBondsDescriptor(
        self,
        name,
        bonds=None,
        bounds=(None, None),
        integer=False,
        binary=False,
        rules=None,
        symmetric_bonds=False,
    ):
        """Method to add a bond-indexed descriptor to the model.

        Args:
            name (string): A unique (otherwise Pyomo will complain) name.
            bonds (list<tuple<int,int>>): Optional, subset of canvas neighbor
                pairs to index the new descriptor over.
                Default: None (i.e., all neighbor pairs included)
            bounds (tuple/dict/func): If tuple, the lower and upper bounds on the
                descriptor values across all indices. If dict, the bounds can be
                individually set for each index. Otherwise, advanced users can
                specify a function to be interpreted by Pyomo.
            integer (bool): Flag to indicate if the descriptor is integer.
            binary (bool): Flag to indicate if the descriptor is boolean.
            rules (list<DescriptorRules>): List of rules to define and constrain
                the material descriptor design space.
        """
        assert not hasattr(self, name)
        bonds = (
            bonds
            if bonds is not None
            else [
                (i, j)
                for i in range(len(self.canv))
                for j in self.canv.NeighborhoodIndexes[i]
                if (j is not None and (not symmetric_bonds or j > i))
            ]
        )
        Desc = MaterialDescriptor(
            name=name,
            canv=self.canv,
            bonds=bonds,
            bounds=bounds,
            integer=integer,
            binary=binary,
            rules=rules,
        )
        setattr(self, name, Desc)
        self._descriptors.append(Desc)

    def addNeighborsDescriptor(
        self,
        name,
        sites=None,
        bounds=(None, None),
        integer=False,
        binary=False,
        rules=None,
    ):
        """Method to add a neighborhood-indexed descriptor to the model.

        Args:
            name (string): A unique (otherwise Pyomo will complain) name.
            sites (list<int>): Optional, subset of canvas sites to index the new
                descriptor over.
                Default: None (i.e., all sites in canvas are considered)
            bounds (tuple/dict/func): If tuple, the lower and upper bounds on the
                descriptor values across all indices. If dict, the bounds can be
                individually set for each index. Otherwise, advanced users can
                specify a function to be interpreted by Pyomo.
            integer (bool): Flag to indicate if the descriptor is integer.
            binary (bool): Flag to indicate if the descriptor is boolean.
            rules (list<DescriptorRules>): List of rules to define and constrain
                the material descriptor design space.
        """
        assert not hasattr(self, name)
        sites = sites if sites is not None else list(range(len(self.canv)))
        Desc = MaterialDescriptor(
            name=name,
            canv=self.canv,
            sites=sites,
            bounds=bounds,
            integer=integer,
            binary=binary,
            rules=rules,
        )
        setattr(self, name, Desc)
        self._descriptors.append(Desc)

    def addGlobalTypesDescriptor(
        self,
        name,
        site_types=None,
        bond_types=None,
        bounds=(None, None),
        integer=False,
        binary=False,
        rules=None,
    ):
        """Method to add a type-indexed descriptor to the model.

        Args:
            name (string): A unique (otherwise Pyomo will complain) name.
            site_types (list<BBlock>): Optional, subset of building block types
                to index the new descriptor over.
                Note: If both site_types and bond_types are left to None, then
                we decide to index over building block types by default.
            bond_types (list<tuple<BBlock,BBlock>>): Optional, subset of building
                block pairs to index the new descriptor over.
                Note: If both site_types and bond_types are left to None, then
                we decide to index over building block types by default.
            bounds (tuple/dict/func): If tuple, the lower and upper bounds on the
                descriptor values across all indices. If dict, the bounds can be
                individually set for each index. Otherwise, advanced users can
                specify a function to be interpreted by Pyomo.
            integer (bool): Flag to indicate if the descriptor is integer.
            binary (bool): Flag to indicate if the descriptor is boolean.
            rules (list<DescriptorRules>): List of rules to define and constrain
                the material descriptor design space.
        """
        assert not hasattr(self, name)
        site_types = (
            site_types
            if site_types is not None
            else (self.atoms if bond_types is None else None)
        )
        bond_types = bond_types
        Desc = MaterialDescriptor(
            name=name,
            atoms=self.atoms,
            site_types=site_types,
            bond_types=bond_types,
            bounds=bounds,
            integer=integer,
            binary=binary,
            rules=rules,
        )
        setattr(self, name, Desc)
        self._descriptors.append(Desc)

    def addSitesTypesDescriptor(
        self,
        name,
        sites=None,
        site_types=None,
        bounds=(None, None),
        integer=False,
        binary=False,
        rules=None,
    ):
        """Method to add a site-and-type-indexed descriptor to the model.

        Args:
            name (string): A unique (otherwise Pyomo will complain) name.
            sites (list<int>): Optional, subset of canvas sites to index the new
                descriptor over.
                Default: None (i.e., all sites in canvas are considered)
            site_types (list<BBlock>): Optional, subset of building block types
                to index the new descriptor over.
                Default: None (i.e., all building block types are considered)
            bounds (tuple/dict/func): If tuple, the lower and upper bounds on the
                descriptor values across all indices. If dict, the bounds can be
                individually set for each index. Otherwise, advanced users can
                specify a function to be interpreted by Pyomo.
            integer (bool): Flag to indicate if the descriptor is integer.
            binary (bool): Flag to indicate if the descriptor is boolean.
            rules (list<DescriptorRules>): List of rules to define and constrain
                the material descriptor design space.
        """
        assert not hasattr(self, name)
        sites = sites if sites is not None else list(range(len(self.canv)))
        site_types = site_types if site_types is not None else self.atoms
        Desc = MaterialDescriptor(
            name=name,
            atoms=self.atoms,
            canv=self.canv,
            sites=sites,
            site_types=site_types,
            bounds=bounds,
            integer=integer,
            binary=binary,
            rules=rules,
        )
        setattr(self, name, Desc)
        self._descriptors.append(Desc)

    def addBondsTypesDescriptor(
        self,
        name,
        bonds=None,
        bond_types=None,
        bounds=(None, None),
        integer=False,
        binary=False,
        rules=None,
        symmetric_bonds=False,
    ):
        """Method to add a bond-and-type-indexed descriptor to the model.

        Args:
            name (string): A unique (otherwise Pyomo will complain) name.
            bonds (list<tuple<int,int>>): Optional, subset of canvas neighbor
                pairs to index the new descriptor over.
                Default: None (i.e., all neighbor pairs included)
            bond_types (list<tuple<BBlock,BBlock>>): Optional, subset of
                building block pairs to index the new descriptor over.
                Default: None (i.e., all pairs of building blocks considered)
            bounds (tuple/dict/func): If tuple, the lower and upper bounds on the
                descriptor values across all indices. If dict, the bounds can be
                individually set for each index. Otherwise, advanced users can
                specify a function to be interpreted by Pyomo.
            integer (bool): Flag to indicate if the descriptor is integer.
            binary (bool): Flag to indicate if the descriptor is boolean.
            rules (list<DescriptorRules>): List of rules to define and constrain
                the material descriptor design space.
        """
        assert not hasattr(self, name)
        bonds = (
            bonds
            if bonds is not None
            else [
                (i, j)
                for i in range(len(self.canv))
                for j in self.canv.NeighborhoodIndexes[i]
                if (j is not None and (not symmetric_bonds or j > i))
            ]
        )
        bond_types = (
            bond_types
            if bond_types is not None
            else [(k, l) for k in self.atoms for l in self.atoms]
        )
        Desc = MaterialDescriptor(
            name=name,
            atoms=self.atoms,
            canv=self.canv,
            bonds=bonds,
            bond_types=bond_types,
            bounds=bounds,
            integer=integer,
            binary=binary,
            rules=rules,
        )
        setattr(self, name, Desc)
        self._descriptors.append(Desc)

    def addNeighborsTypesDescriptor(
        self,
        name,
        sites=None,
        bond_types=None,
        bounds=(None, None),
        integer=False,
        binary=False,
        rules=None,
    ):
        """Method to add a neighborhood-bond-type-indexed descriptor.

        Args:
            name (string): A unique (otherwise Pyomo will complain) name.
            sites (list<int>): Optional, subset of canvas sites to index the new
                descriptor over.
                Default: None (i.e., all sites in canvas are considered)
            bond_types (list<tuple<BBlock,BBlock>>): Optional, subset of building
                block pairs to index the new descriptor over.
                Default: None (i.e., all pairs of building blocks considered)
            bounds (tuple/dict/func): If tuple, the lower and upper bounds on the
                descriptor values across all indices. If dict, the bounds can be
                individually set for each index. Otherwise, advanced users can
                specify a function to be interpreted by Pyomo.
            integer (bool): Flag to indicate if the descriptor is integer.
            binary (bool): Flag to indicate if the descriptor is boolean.
            rules (list<DescriptorRules>): List of rules to define and constrain
                the material descriptor design space.
        """
        assert not hasattr(self, name)
        sites = sites if sites is not None else list(range(len(self.canv)))
        bond_types = (
            bond_types
            if bond_types is not None
            else [(k, l) for k in self.atoms for l in self.atoms]
        )
        Desc = MaterialDescriptor(
            name=name,
            atoms=self.atoms,
            canv=self.canv,
            sites=sites,
            bond_types=bond_types,
            bounds=bounds,
            integer=integer,
            binary=binary,
            rules=rules,
        )
        setattr(self, name, Desc)
        self._descriptors.append(Desc)

    def addSitesConfsDescriptor(
        self,
        name,
        sites=None,
        confs=None,
        bounds=(0, 1),
        integer=True,
        binary=True,
        rules=None,
    ):
        """Method to add a site-and-conformation-indexed descriptor.

        Args:
            name (string): A unique (otherwise Pyomo will complain) name.
            sites (list<int>): Optional, subset of canvas sites to index the new
                descriptor over.
                Default: None (i.e., all sites in canvas are considered)
            confs (list<int>): Optional, subset of conformation indices to index
                the new descriptor over.
                Default: None (i.e., all conformations included)
            bounds (tuple/dict/func): If tuple, the lower and upper bounds on the
                descriptor values across all indices. If dict, the bounds can be
                individually set for each index. Otherwise, advanced users can
                specify a function to be interpreted by Pyomo.
            integer (bool): Flag to indicate if the descriptor is integer.
            binary (bool): Flag to indicate if the descriptor is boolean.
            rules (list<DescriptorRules>): List of rules to define and constrain
                the material descriptor design space.
        """
        sites = sites if sites is not None else list(range(len(self.canv)))
        confs = (
            confs
            if confs is not None
            else (list(range(len(self.confDs))) if self.confDs is not None else None)
        )
        Desc = MaterialDescriptor(
            name=name,
            canv=self.canv,
            confDs=self.confDs,
            sites=sites,
            confs=confs,
            bounds=bounds,
            integer=integer,
            binary=binary,
            rules=rules,
        )
        setattr(self, name, Desc)
        self._descriptors.append(Desc)

    # === PROPERTY EVALUATION METHODS
    def maximize(self, func, **kwargs):
        """Method to maximize a target functionality of the material model.

        Args:
            func (``MaterialDescriptor``/``Expr``): Material functionality to optimize.
            **kwargs: Arguments to ``MatOptModel.optimize``

        Returns:
            (``Design``/list<``Design``>) Optimal designs.

        Raises:
            ``pyomo.common.errors.ApplicationError`` if MatOpt can not find
            usable solver (CPLEX or NEOS-CPLEX)

        See ``MatOptModel.optimize`` method for details.
        """
        return self.optimize(func, sense=maximize, **kwargs)

    def minimize(self, func, **kwargs):
        """Method to minimize a target functionality of the material model.

        Args:
            func (``MaterialDescriptor``/``Expr``): Material functionality to optimize.
            **kwargs: Arguments to ``MatOptModel.optimize``

        Returns:
            (``Design``/list<``Design``>) Optimal designs.

        Raises:
            ``pyomo.common.errors.ApplicationError`` if MatOpt can not find usable
            solver (CPLEX or NEOS-CPLEX)

        See ``MatOptModel.optimize`` method for details.
        """
        return self.optimize(func, sense=minimize, **kwargs)

    def optimize(
        self,
        func,
        sense,
        nSolns=1,
        tee=True,
        disp=1,
        keepfiles=False,
        tilim=3600,
        trelim=None,
        solver="cplex",
    ):
        """Method to create and optimize the materials design problem.

        This method automatically creates a new optimization model every
        time it is called. Then, it solves the model via Pyomo with the
        CPLEX solver.

        If multiple solutions (called a 'solution pool') are desired, then
        the nSolns argument can be provided and the populate method will
        be called instead.

        Args:
            func (``MaterialDescriptor``/``Expr``): Material functionality to optimize.
            sense (int): flag to indicate the choice to minimize or maximize the
                functionality of interest.
                Choices: minimize/maximize (Pyomo constants 1,-1 respectively)
            nSolns (int): Optional, number of Design objects to return.
                Default: 1 (See ``MatOptModel.populate`` for more information)
            tee (bool): Optional, flag to turn on solver output.
                Default: True
            disp (int): Optional, flag to control level of MatOpt output.
                Choices: 0: No MatOpt output (other than solver tee) 1: MatOpt
                output for outer level method 2: MatOpt output for solution pool &
                individual solns. Default: 1
            keepfiles (bool): Optional, flag to save temporary pyomo files.
                Default: True
            tilim (float): Optional, solver time limit (in seconds).
                Default: 3600
            trelim (float): Optional, solver tree memeory limit (in MB).
                Default: None (i.e., Pyomo/CPLEX default)
            solver (str): Solver choice. Currently only cplex or neos-cplex are supported
                Default: cplex

        Returns:
            (``Design``/list<``Design``>) Optimal design or designs, depending
            on the number of solutions requested by argument ``nSolns``.

        Raises:
            ``pyomo.common.errors.ApplicationError`` if MatOpt can not find
            usable solver (CPLEX or NEOS-CPLEX)
        """
        if nSolns > 1:
            return self.populate(
                func,
                sense=sense,
                nSolns=nSolns,
                tee=tee,
                disp=disp,
                keepfiles=keepfiles,
                tilim=tilim,
                trelim=trelim,
                solver=solver,
            )
        elif nSolns == 1:
            self._pyomo_m = self._make_pyomo_model(func, sense)
            return self.__solve_pyomo_model(tee, disp, keepfiles, tilim, trelim, solver)

    def populate(
        self,
        func,
        sense,
        nSolns,
        tee=True,
        disp=1,
        keepfiles=False,
        tilim=3600,
        trelim=None,
        solver="cplex",
    ):
        """Method to a pool of solutions that optimize the material model.

        This method automatically creates a new optimization model every
        time it is called. Then, it solves the model via Pyomo with the
        CPLEX solver.

        The populate method iteratively solves the model, interprets the
        solution as a Design object, creates a constraint to disallow that
        design, and resolves to find the next best design. We build a pool
        of Designs that are gauranteed to be the nSolns-best solutions in the
        material design space.

        Args:
            func (``MaterialDescriptor``/``Expr``): Material functionality to optimize.
            sense (int): flag to indicate the choice to minimize or maximize
                the functionality of interest.
                Choices: minimize/maximize (Pyomo constants 1,-1 respectively)
            nSolns (int): Optional, number of Design objects to return.
                Default: 1 (See ``MatOptModel.populate`` for more information)
            tee (bool): Optional, flag to turn on solver output.
                Default: True
            disp (int): Optional, flag to control level of MatOpt output.
                Choices: 0: No MatOpt output (other than solver tee) 1: MatOpt
                output for outer level method 2: MatOpt output for solution
                pool & individual solns. Default: 1
            keepfiles (bool): Optional, flag to save temporary pyomo files.
                Default: True
            tilim (float): Optional, solver time limit (in seconds).
                Default: 3600
            trelim (float): Optional, solver tree memeory limit (in MB).
                Default: None (i.e., Pyomo/CPLEX default)
            solver (str): Solver choice. Currently only cplex or neos-cplex are
                supported Default: cplex

        Returns:
            (list<``Design``>) A list of optimal Designs in order of decreasing
            optimality.

        Raises:
            ``pyomo.common.errors.ApplicationError`` if MatOpt can not find
            usable solver (CPLEX or NEOS-CPLEX)
        """
        self._pyomo_m = self._make_pyomo_model(func, sense)
        self._pyomo_m.iSolns = Set(initialize=list(range(nSolns)))
        self._pyomo_m.IntCuts = Constraint(self._pyomo_m.iSolns)

        def dispPrint(*args):
            if disp > 0:
                print(*args)
            else:
                pass

        Ds = []
        for iSoln in range(nSolns):
            dispPrint("Starting populate for solution #{}... ".format(iSoln))
            D = self.__solve_pyomo_model(
                tee, disp - 1, keepfiles, tilim, trelim, solver
            )
            if D is not None:
                dispPrint(
                    "Found solution with objective: {}".format(value(self._pyomo_m.obj))
                )
                Ds.append(D)
                if len(self._pyomo_m.Yik) > 0:
                    self._pyomo_m.IntCuts.add(
                        index=iSoln, expr=(Disallow(D)._pyomo_expr(self.Yik) >= 1)
                    )
                elif len(self._pyomo_m.Yi) > 0:
                    self._pyomo_m.IntCuts.add(
                        index=iSoln, expr=(Disallow(D)._pyomo_expr(self.Yi) >= 1)
                    )
                else:
                    raise NotImplementedError("Decide what to do " "in this case...")
            else:
                dispPrint("No solution found. Terminating populate.")
                break
        dispPrint("Identified {} solutions via populate.".format(len(Ds)))
        return Ds

    def _make_pyomo_model(self, obj_expr, sense):
        """Method to create a Pyomo concrete model object.

        This method creates a Pyomo model and also modifies several objects
        in the MatOpt framework. It creates Pyomo variable objects and
        attches references to those variables on each of the
        MaterialDescriptors attached to the MatOptModel.

        Args:
            obj_expr (MaterialDescriptor/Expr): Material functionality to
                optimize.
            sense (int): flag to indicate the choice to minimize or maximize the
                functionality of interest.
                Choices: minimize/maximize (Pyomo constants 1,-1 respectively)

        Returns:
            (ConcreteModel) Pyomo model object.
        """
        m = makeMyPyomoBaseModel(self.canv, Atoms=self.atoms, Confs=self.confDs)
        self.Yi._pyomo_var = m.Yi
        self.Xij._pyomo_var = m.Xij
        self.Ci._pyomo_var = m.Ci
        self.Yik._pyomo_var = m.Yik
        self.Xijkl._pyomo_var = m.Xijkl
        self.Cikl._pyomo_var = m.Cikl
        self.Zic._pyomo_var = m.Zic
        for desc in self._descriptors:
            if desc.name not in ("Yik", "Yi", "Xijkl", "Xij", "Cikl", "Ci", "Zic"):
                v = Var(
                    *desc.index_sets,
                    domain=(
                        Binary if desc.binary else (Integers if desc.integer else Reals)
                    ),
                    bounds=desc._pyomo_bounds,
                    dense=False
                )
                setattr(m, desc.name, v)
                setattr(desc, "_pyomo_var", v)
        for desc in self._descriptors:
            for c, pyomo_con in enumerate(desc._pyomo_cons(m)):
                setattr(m, "Assign{}_{}".format(desc.name, c), pyomo_con)
        if sum(obj_expr.dims) == 0:
            m.obj = Objective(expr=obj_expr._pyomo_expr(index=(None,)), sense=sense)
        else:
            raise TypeError(
                "The MaterialDescriptor chosen is not supported to be an objective, please contact MatOpt "
                "developer for potential fix"
            )
        # NOTE: The timing of the call to addConsForGeneralVars is important
        #       We need to call it after all user-defined descriptors are
        #       encoded.
        #       Else, lots of constraints for basic variables that are not
        #       necessary will be written.
        addConsForGeneralVars(m)
        for desc in self._descriptors:
            for r in desc.rules:
                if isinstance(r, FixedTo):
                    desc._fix_pyomo_var_by_rule(r, m)
        return m

    def __solve_pyomo_model(self, tee, disp, keepfiles, tilim, trelim, solver):
        """Method to solve the formulated Pyomo optimization model.

        This function is intended to standardize the printout and
        approach for converting results into Designs that is used
        by the optimize and populate methods.

        Args:
            tee (bool): Flag to turn on solver output.
            disp (int): Flag to control level of MatOpt output.
            keepfiles (bool): Flag to save temporary pyomo files.
            tilim (float): Solver time limit (in seconds).
            trelim (float): Solver tree memeory limit (in MB).
            solver (str): Solver choice. Currently only cplex or
                neos-cplex are supported

        Returns:
            (Design) The best design identified by the solver, if any.
                The quality of the solution (optimal vs best found at
                termination) can be found by reading the output display.
                In the case that the model was infeasible or no solution
                could be identified, the method returns 'None'.
        """
        if solver == "cplex":
            opt = SolverFactory("cplex")
            opt.options["mip_tolerances_absmipgap"] = 0.0
            opt.options["mip_tolerances_mipgap"] = 0.0
            if tilim is not None:
                opt.options["timelimit"] = tilim
            if trelim is not None:
                opt.options["mip_limits_treememory"] = trelim
            res = opt.solve(
                self._pyomo_m, tee=tee, symbolic_solver_labels=True, keepfiles=keepfiles
            )
        elif solver == "neos-cplex":
            with SolverManagerFactory("neos") as manager:
                opt = SolverFactory("cplex")
                opt.options["absmipgap"] = 0.0  # NOTE: different option names
                opt.options["mipgap"] = 0.0
                if tilim is not None:
                    opt.options["timelimit"] = tilim
                if trelim is not None:
                    opt.options["treememory"] = trelim
                res = manager.solve(self._pyomo_m, opt=opt)
        else:
            raise NotImplementedError(
                "MatOpt is tailored to perform best with CPLEX (locally or through NEOS), "
                "please contact MatOpt developer for additional solver support "
            )
        solver_status = res.solver.status
        solver_term = res.solver.termination_condition
        soln_status = res.solution.status
        has_solution = (
            soln_status == SolutionStatus.optimal
            or soln_status == SolutionStatus.feasible
            or soln_status == SolutionStatus.bestSoFar
            or soln_status == SolutionStatus.globallyOptimal
            or soln_status == SolutionStatus.locallyOptimal
        )
        # NOTE: The block below is a hack to get around the fact that Pyomo
        #       solution statuses are not flagged correctly all the time. If
        #       solution status was unknown (but actually optimal or feasible)
        #       then this should (hopefully) flag the solution as available
        if soln_status == SolutionStatus.unknown:
            value(self._pyomo_m.obj)
            has_solution = True

        def dispPrint(*args):
            if disp > 0:
                print(*args)
            else:
                pass

        if solver_status == SolverStatus.ok:
            dispPrint("The solver exited normally.")
            if (
                solver_term == TerminationCondition.optimal
                or solver_term == TerminationCondition.locallyOptimal
                or solver_term == TerminationCondition.globallyOptimal
            ):
                #  NOTE: This assertion should be re-enabled when Pyomo bug
                #        described above is fixed.
                # assert(soln_status==SolutionStatus.optimal)
                dispPrint("A feasible and provably optimal solution " "is available.")
            else:
                dispPrint(
                    "The solver exited due to termination criteria: {}".format(
                        solver_term
                    )
                )
                if has_solution:
                    dispPrint(
                        "A feasible (but not provably optimal) "
                        "solution is available."
                    )
                else:
                    dispPrint("No solution available.")
        else:
            dispPrint(
                "The solver did not exit normally. Status: {}".format(solver_status)
            )
        if has_solution:
            dispPrint("The Design has objective: {}".format(value(self._pyomo_m.obj)))
            result = Design(self.canv)
            setDesignFromModel(result, self._pyomo_m)
            return result
        else:
            return None

    # === BASIC QUERY METHODS
    @property
    def canv(self):
        return self._canv

    @property
    def atoms(self):
        return self._atoms

    @property
    def confDs(self):
        return self._confDs

    @property
    def descriptors(self):
        return self._descriptors
