##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
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
    def __init__(self, sites=None, bonds=None, site_types=None, bond_types=None,
                 confs=None):
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
        """Constructor of IndexedElem from other IndexedElem objects.

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
            sites = (LHS.sites if LHS.sites is not None else RHS.sites)
        if LHS.bonds is not None and RHS.bonds is not None:
            bonds = list(set(LHS.bonds) & set(RHS.bonds))
        else:
            bonds = (LHS.bonds if LHS.bonds is not None else RHS.bonds)
        if LHS.site_types is not None and RHS.site_types is not None:
            site_types = list(set(LHS.site_types) & set(RHS.site_types))
        else:
            site_types = (LHS.site_types if LHS.site_types is not None
                          else RHS.site_types)
        if LHS.bond_types is not None and RHS.bond_types is not None:
            bond_types = list(set(LHS.bond_types) & set(RHS.bond_types))
        else:
            bond_types = (LHS.bond_types if LHS.bond_types is not None
                          else RHS.bond_types)
        if LHS.confs is not None and RHS.confs is not None:
            confs = list(set(LHS.confs) & set(RHS.confs))
        else:
            confs = (LHS.confs if LHS.confs is not None else RHS.confs)
        return cls(sites=sites, bonds=bonds,
                   site_types=site_types, bond_types=bond_types,
                   confs=confs)

    def mask(self, index, Comb):
        """Method to identify the indexes relevant to this object.

        Given an instance of index that was generated by another IndexedElem 
        object (Comb), we identify which parts of that index were relevant
        to this object. 

        Example: 
        VarIndexes = IndexedElem(sites=[1,2])
        CoefIndexes = IndexedElem(site_types=['A','B'])
        Comb = IndexedElem.fromComb(VarIndexes,CoefIndexes)
        for k in Comb.keys():
            site = VarIndexes.mask(k,Comb)
            site_type = CoefIndexes.mask(k,Comb)

        Args:
            index (tuple<int/BBlock>): index from which to identify relevant 
                parts
            Comb (IndexedElem): object from which the index was generated

        Returns:
            (tuple<int/BBlock>) index with indices relevant to this object 
                remaining
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

        return (self.sites is not None,
                self.bonds is not None,
                self.site_types is not None,
                self.bond_types is not None,
                self.confs is not None)

    @property
    def index_sets(self):
        """Sets (actually lists) of indices.

        Note that in the cases that there are no relevant indices, a
        dummy list [[None]] is returned to allow the [None] key to be 
        included.
        
        Returns:
            (list<list<int/BBlock>>) lists of indices of each relevant type.
        """

        result = [s for s in (self.sites,
                              self.bonds,
                              self.site_types,
                              self.bond_types,
                              self.confs) if s is not None]
        if not result:
            result = [[None]]
        return result

    @property
    def index_dict(self):
        """Dictionary of relevant attributes.

        Returns:
            (dict<string:list<int/BBlock>>) attributes of this IndexedElem.
        """

        return {'sites': self.sites,
                'bonds': self.bonds,
                'site_types': self.site_types,
                'bond_types': self.bond_types,
                'confs': self.confs}

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
            (generator<tuple<int/BBlock>>) keys generator 
                (similar to dict.keys()).
        """

        index_sets = self.index_sets
        if len(index_sets) > 1:
            return product(*self.index_sets)
        elif len(index_sets) == 1:
            return (k for k in index_sets[0])
        else:
            raise NotImplementedError('There should always be at least '
                                      'a [None] key')


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
        vals (list<float>/dict/other) Any data structure that supports the
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
            each descriptor by.
            Default: 1.0
        descs (Descriptor/list<Descriptor>): descriptors to add 
        offset (float/int): Optional, scalar value to add to the expression
            Default: 0.0
        **kwargs: Optional, index information passed to IndexedElem if 
            interested in a subset of indices
            Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.coefs = (coefs if type(coefs) is list else [coefs])
        self.descs = (descs if type(descs) is list else [descs])
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
        many terms)
    (index information inherited from IndexedElem)    
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, coefi, desci, coefj=None, descj=None, offset=0.0,
                 symmetric_bonds=False, **kwargs):
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
        self.coefj = (coefj if coefj is not None else coefi)
        self.descj = (descj if descj is not None else desci)
        self.offset = offset
        self.symmetric_bonds = symmetric_bonds
        if 'bonds' not in kwargs:
            kwargs['bonds'] = [(i, j)
                               for i in desci.sites
                               for j in desci.canv.NeighborhoodIndexes[i]
                               if (j is not None and (not symmetric_bonds
                                                      or j > i))]
        if 'bond_types' in kwargs:
            pass  # use the kwargs bond_types
        elif coefi.bond_types is not None:
            kwargs['bond_types'] = coefi.bond_types
        elif desci.site_types is not None:
            kwargs['bond_types'] = [(k, l)
                                    for k in desci.site_types
                                    for l in desci.site_types]
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
            raise NotImplementedError('Decide how to split the extra '
                                      'indices in this case...')
        if (type(self.coefi) is float or
                type(self.coefi) is int or
                type(self.coefi) is SimpleParam):
            ci = self.coefi
        else:
            coefi_index = self.coefi.mask((i, i, j, k, k, l),
                                          IndexedElem(sites=[i],
                                                      bonds=[(i, j)],
                                                      site_types=[k],
                                                      bond_types=[(k, l)]))
            ci = self.coefi[coefi_index]

        if (type(self.coefj) is float or
                type(self.coefj) is int or
                type(self.coefj) is SimpleParam):
            cj = self.coefj
        else:
            coefj_index = self.coefj.mask((j, j, i, l, l, k),
                                          IndexedElem(sites=[j],
                                                      bonds=[(j, i)],
                                                      site_types=[l],
                                                      bond_types=[(l, k)]))
            cj = self.coefj[coefj_index]
        desci_index = self.desci.mask((i, i, j, k, k, l),
                                      IndexedElem(sites=[i],
                                                  bonds=[(i, j)],
                                                  site_types=[k],
                                                  bond_types=[(k, l)]))
        descj_index = self.descj.mask((j, j, i, l, l, k),
                                      IndexedElem(sites=[j],
                                                  bonds=[(j, i)],
                                                  site_types=[l],
                                                  bond_types=[(l, k)]))
        di = self.desci._pyomo_expr(index=desci_index)
        dj = self.descj._pyomo_expr(index=descj_index)
        return self.offset + ci * di + cj * dj


class SumNeighborSites(Expr):
    """A class for expressions for summation across neighbor sites.

    Attributes:
    desc (Descriptor): descriptors to sum around a site
    coefs (float/list<float>): Optional, coefficients to multiple each 
        neighbor term by. 
        Default=1.0
    offset: Optional, term to add to the expression. 
        Default=0.0
    (index information inherited from IndexedElem)    
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, desc, coefs=1.0, offset=0.0, **kwargs):
        """Standard constructor for expressions of neighbor summations.

        Args: 
        desc (Descriptor): descriptors to sum around a site
        coefs (float/list<float>): Optional, coefficients to multiple each 
            neighbor term by. 
            Default=1.0
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
                result += (self.coefs if
                           (type(self.coefs) is float or
                            type(self.coefs) is int or
                            type(self.coefs) is SimpleParam)
                           else self.coefs[n]) * self.desc._pyomo_var[(j, *index)]
        return result


class SumNeighborBonds(Expr):
    """A class for expressions from summation of neighbor bond descriptors.

    Attributes:
    desc (Descriptor/Expr): descriptors to sum over
    coefs (float/list<float>): coefficients to multiply bonds to 
        neighbor sites 
    offset (float): coefficient to add to the expression
    symmetric_bonds (bool): flag to indicate if bond variables should be
        added in a symmetric way
    (index information inherited from IndexedElem)    
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, desc, coefs=1.0, offset=0.0,
                 symmetric_bonds=False, **kwargs):
        """Standard constructor for summation of neighboring bonds

        Args:
        desc (Descriptor): descriptors to sum around a site
        coefs (float/list<float>): Optional, coefficients to multiple each 
            neighbor term by. 
            Default=1.0
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
                result += (self.coefs if
                           (type(self.coefs) is float or
                            type(self.coefs) is int or
                            type(self.coefs) is SimpleParam)
                           else
                           self.coefs[n]) * self.desc._pyomo_expr(index=(i, j, *index))
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
    def __init__(self, desc, coefs=1.0, offset=0.0,
                 sites_to_sum=None, **kwargs):
        """Standard constructor for summation of site contributions.

        Args:
        desc (Descriptor): descriptors or expressions to sum across all sites
        coefs (float/list<float>): Optional, coefficients to multiple each 
            site term by. 
            Default=1.0
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
        self.sites_to_sum = (sites_to_sum if sites_to_sum is not None
                             else desc.sites)
        kwargs = {**desc.index_dict, **kwargs}
        kwargs.pop('sites')
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
            result += (self.coefs if
                       (type(self.coefs) is float or
                        type(self.coefs) is int or
                        type(self.coefs) is SimpleParam)
                       else
                       self.coefs[(i, *index)]) * self.desc._pyomo_var[(i, *index)]
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
    def __init__(self, desc, coefs=1.0, offset=0.0,
                 bonds_to_sum=None, **kwargs):
        """Standard constructor for summation of bond contributions.

        Args:
        desc (Descriptor): descriptors or expressions to sum across all bonds
        coefs (float/list<float>): Optional, coefficients to multiple each 
            bond term by. 
            Default=1.0
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
        self.bonds_to_sum = (bonds_to_sum if bonds_to_sum is not None
                             else desc.bonds)
        kwargs = {**desc.index_dict, **kwargs}
        kwargs.pop('bonds')
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
            result += ((self.coefs if
                        (type(self.coefs) is float or
                         type(self.coefs) is int or
                         type(self.coefs) is SimpleParam)
                        else self.coefs[(i, j, *index)])
                       * self.desc._pyomo_var[(i, j, *index)])
        return result


class SumSiteTypes(Expr):
    """A class for expressions formed by summation over building block types.

    Attributes:
    desc (Descriptor/Expr): descriptors or expressions to sum over
    coefs (float/list<float>): coefficients to multiply contributions 
        from each building block type
    offset (float): coefficient to add to the expression
    site_types_to_sum (list<BBlock>): building block types to consider in 
        the summation
    (index information inherited from IndexedElem)    
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, desc, coefs=1.0, offset=0.0,
                 site_types_to_sum=None, **kwargs):
        """Standard constructor for summation of contributions by site-type.

        Args:
        desc (Descriptor): descriptors or expressions to sum across site types
        coefs (float/list<float>): Optional, coefficients to multiple each 
            site-type term by. 
            Default=1.0
        offset (float): Optional, coefficient to add to the expression. 
            Default=0.0
        bonds_types_to_sum (list<int>): Optional, subset of site types
            to sum.
            Default=None, meaning all site-types in the desc object are 
                considered.
        **kwargs: Optional, index information passed to IndexedElem if 
            interested in a subset of indices
            Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.desc = desc
        self.coefs = coefs
        self.offset = offset
        self.site_types_to_sum = (site_types_to_sum
                                  if site_types_to_sum is not None
                                  else desc.site_types)
        kwargs = {**desc.index_dict, **kwargs}
        kwargs.pop('site_types')
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
        assert (index is not None)
        i, *index = index
        if index == (None,):
            index = ()
        result = self.offset
        for k in self.site_types_to_sum:
            result += ((self.coefs if
                        (type(self.coefs) is float or
                         type(self.coefs) is int or
                         type(self.coefs) is SimpleParam)
                        else self.coefs[(i, k, *index)])
                       * self.desc._pyomo_var[(i, k, *index)])
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
    def __init__(self, desc, coefs=1.0, offset=0.0,
                 bond_types_to_sum=None, **kwargs):
        """Standard constructor for summation of contributions by bond-type.

        Args:
        desc (Descriptor): descriptors or expressions to sum across bond types
        coefs (float/list<float>): Optional, coefficients to multiple each 
            bond-type term by. 
            Default=1.0
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
        self.bond_types_to_sum = (bond_types_to_sum
                                  if bond_types_to_sum is not None
                                  else desc.bond_types)
        kwargs = {**desc.index_dict, **kwargs}
        kwargs.pop('bond_types')
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
        assert (index is not None)
        i, j, *index = index
        if index == (None,):
            index = ()
        result = self.offset
        for k, l in self.bond_types_to_sum:
            result += ((self.coefs if
                        (type(self.coefs) is float or
                         type(self.coefs) is int or
                         type(self.coefs) is SimpleParam)
                        else self.coefs[(i, j, k, l, *index)])
                       * self.desc._pyomo_var[(i, j, k, l, *index)])
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
        the summation
    (index information inherited from IndexedElem)    
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, desc, coefs=1.0, offset=0.0,
                 sites_to_sum=None,
                 site_types_to_sum=None, **kwargs):
        """Standard constructor for summation of site contributions.

        Args:
        desc (Descriptor): descriptors or expressions to sum across all 
            sites and site types
        coefs (float/list<float>): Optional, coefficients to multiple each 
            site and site-type term by. 
            Default=1.0
        offset (float): Optional, coefficient to add to the expression. 
            Default=0.0
        sites_to_sum (list<int>): Optional, subset of canvas sites to sum.
            Default=None, meaning all sites in the desc object are considered.
        site_types_to_sum (list<BBlock>): Optional, subset of site types to 
            sum.
            Default=None, meaning all site types in the desc object are 
                considered. 
        **kwargs: Optional, index information passed to IndexedElem if 
            interested in a subset of indices
            Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.desc = desc
        self.coefs = coefs
        self.offset = offset
        self.sites_to_sum = (sites_to_sum if sites_to_sum is not None
                             else desc.sites)
        self.site_types_to_sum = (site_types_to_sum
                                  if site_types_to_sum is not None
                                  else desc.site_types)
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
                result += ((self.coefs if
                            (type(self.coefs) is float or
                             type(self.coefs) is int or
                             type(self.coefs) is SimpleParam)
                            else self.coefs[(i, k, *index)])
                           * self.desc._pyomo_var[(i, k, *index)])
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
    def __init__(self, desc, coefs=1.0, offset=0.0,
                 bonds_to_sum=None, bond_types_to_sum=None, **kwargs):
        """Standard constructor for summation of contributions by bond-type.

        Args:
        desc (Descriptor): descriptors or expressions to sum across bonds 
            and bond types
        coefs (float/list<float>): Optional, coefficients to multiple each 
            term by. 
            Default=1.0
        offset (float): Optional, coefficient to add to the expression. 
            Default=0.0
        bonds_to_sum (list<int>): Optional, subset of bonds to sum. 
            Default=None, meaning all bonds in the desc object are considered.
        bonds_types_to_sum (list<int>): Optional, subset of bond types
            to sum.
            Default=None, meaning all bond-types in the desc object are 
                considered.
        **kwargs: Optional, index information passed to IndexedElem if 
            interested in a subset of indices
            Possible choices: sites, bonds, site_types, bond_types, confs.
        """
        self.desc = desc
        self.coefs = coefs
        self.offset = offset
        self.bonds_to_sum = (bonds_to_sum if bonds_to_sum is not None
                             else desc.bonds)
        self.bond_types_to_sum = (bond_types_to_sum
                                  if bond_types_to_sum is not None
                                  else desc.bond_types)
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
                result += ((self.coefs if
                            (type(self.coefs) is float or
                             type(self.coefs) is int or
                             type(self.coefs) is SimpleParam)
                            else self.coefs[i, j, k, l])
                           * self.desc._pyomo_var[i, j, k, l])
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
    def __init__(self, Zic, coefs=1.0, offset=0.0,
                 confs_to_sum=None, **kwargs):
        """Standard constructor for summation of bond contributions.

        Args:
        Zic (Descriptor): descriptors or expressions to sum across 
            conformations
        coefs (float/list<float>): Optional, coefficients to multiple each 
            conformation term by. 
            Default=1.0
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
        self.confs_to_sum = (confs_to_sum if confs_to_sum is not None
                             else Zic.confs)
        kwargs = {**Zic.index_dict, **kwargs}
        kwargs.pop('confs')
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
            result += ((self.coefs if
                        (type(self.coefs) is float or
                         type(self.coefs) is int or
                         type(self.coefs) is SimpleParam)
                        else self.coefs[(i, c, *index)])
                       * self.Zic._pyomo_var[(i, c, *index)])
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
    def __init__(self, Zic, coefs=1.0, offset=0.0,
                 sites_to_sum=None,
                 confs_to_sum=None, **kwargs):
        """Standard constructor for summation of site and conformation 
           contributions.

        Args:
        Zic (Descriptor): descriptors or expressions to sum across all sites 
            and conformations. 
        coefs (float/list<float>): Optional, coefficients to multiple each 
            site and conformation term by. 
            Default=1.0
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
        self.sites_to_sum = (sites_to_sum if sites_to_sum is not None
                             else Zic.sites)
        self.confs_to_sum = (confs_to_sum if confs_to_sum is not None
                             else Zic.confs)
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
                result += ((self.coefs if
                            (type(self.coefs) is float or
                             type(self.coefs) is int or
                             type(self.coefs) is SimpleParam)
                            else self.coefs[(i, c, *index)])
                           * self.Zic._pyomo_var[(i, c, *index)])
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
        self.expr = (LinearExpr(offset=e)
                     if type(e) is float or type(e) is int else e)
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
            return operator(LHS._pyomo_expr(LHS_index),
                            RHS._pyomo_expr(RHS_index))

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

        return SimpleDescriptorRule._pyomo_rule(self, desc,
                                                less_than,
                                                self.expr)


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

        return SimpleDescriptorRule._pyomo_rule(self, desc,
                                                equal_to,
                                                self.expr)


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

        return SimpleDescriptorRule._pyomo_rule(self, desc,
                                                greater_than,
                                                self.expr)


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


