##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Example property package for adsorption of N2 onto a Zeolite.
"""

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           Param,
                           Set,
                           Var,
                           exp,
                           log)
from pyomo.common.config import ConfigBlock, ConfigValue

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        ReactionBlockBase,
                        ProcessBlockData,
                        property_meta)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import (BurntToast,
                                        PropertyNotSupportedError,
                                        PropertyPackageError)

# Some more inforation about this module
__author__ = "Andrew Lee"


# Set up logger
_log = logging.getLogger(__name__)


# TODO: Current reaction blocks assume there is is a config argument named
# property_package which it tries to validate. This does not work so well
# with a heterogeneous system where we have 2 property blocks
@declare_process_block_class("AdsorptionParameterBlock")
class AdsorptionParameterData(ProcessBlockData,
                              property_meta.HasPropertyClassMetadata):
    """
    Property Parameter Block Class
    """
    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare("gas_property_package", ConfigValue(
            description="Reference to associated PropertyPackageParameter "
                        "object for the gas phase.",
            domain=is_physical_parameter_block))
    CONFIG.declare("solid_property_package", ConfigValue(
            description="Reference to associated PropertyPackageParameter "
                        "object for the solid phase.",
            domain=is_physical_parameter_block))
    CONFIG.declare("default_arguments", ConfigBlock(
            description="Default arguments to use with Property Package",
            implicit=True))

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(AdsorptionParameterData, self).build()

        self.reaction_block_class = AdsorptionBlock

        # List of valid phases in property package
        self.phase_list = Set(initialize=['Sol'])

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['Zeo', 'N2'])

        # Reaction Index
        self.rate_reaction_idx = Set(initialize=["R1"])

        # TODO: This is non-standard in form, and uses only the gas phase component list...
        # Reaction Stoichiometry
        self.rate_reaction_stoichiometry = {("R1", "Vap", "O2"): 0,
                                            ("R1", "Vap", "N2"): -1,
                                            ("R1", "Sol", "Zeo"): 0,
                                            ("R1", "Sol", "N2"): 1}

        # Maximum Loading
        self.loading_max = Param(
                self.component_list,
                initialize={"Zeo": 1, "N2": 1},
                doc="Loading of species at complete surface coverage [mol/kg]")

        # Standard Gibbs energy of adsorption
        self.gibbs_mol_rxn = Param(
                initialize=4000,
                doc="Standard molar Gibbs energy of adsorption [J/mol]")

        # Arrhenius constant
        self.arrhenius = Param(initialize=1,
                               doc="Arrhenius constant [?]")

        # Activation energy
        self.energy_activation = Param(initialize=50000,
                                       doc="Activation energy [J/mol]")

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
                'pressure_comp_equil': {'method': '_p_eq', 'units': 'Pa'},
                'k_eq': {'method': '_k_eq', 'units': '?'},
                'k_rxn': {'method': '_k_rxn', 'units': '?'},
                'reaction_rate': {'method': '_rxn_rate', 'units': 'mol/m^3.s'}
                })
        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'kg',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'J',
                               'holdup': 'mol'})


class _ReactionBlock(ReactionBlockBase):
    """
    This Class contains methods which should be applied to Reaction Blocks as a
    whole, rather than individual elements of indexed Reaction Blocks.
    """
    def initialize(blk, outlvl=0, **kwargs):
        '''
        Initialization routine for reaction package.

        Keyword Arguments:
            outlvl : sets output level of initialization routine

                     * 0 = no output (default)
                     * 1 = report after each step

        Returns:
            None
        '''
        if outlvl > 0:
            _log.info('{} Initialization Complete.'.format(blk.name))


@declare_process_block_class("AdsorptionBlock",
                             block_class=_ReactionBlock)
class AdsorptionBlockData(ProcessBlockData):
    """
    An example reaction package for adsorption of N2
    """
    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare("parameters", ConfigValue(
#            domain=is_reaction_parameter_block,
            description="""A reference to an instance of the Reaction Parameter
Block associated with this property package."""))
    CONFIG.declare("solid_state_block", ConfigValue(
#            domain=is_state_block,
            description="""A reference to an instance of a StateBlock for the
solid phase with which this reaction block should be associated."""))
    CONFIG.declare("gas_state_block", ConfigValue(
#            domain=is_state_block,
            description="""A reference to an instance of a StateBlock for the
gas phase with which this reaction block should be associated."""))

    def build(self):
        """
        Callable method for Block construction
        """
        super(AdsorptionBlockData, self).build()

        add_object_reference(self, "_params", self.config.parameters)

        add_object_reference(self,
                             "solid_state_ref",
                             self.config.solid_state_block[self.index()])
        add_object_reference(self,
                             "gas_state_ref",
                             self.config.gas_state_block[self.index()])

    def _p_eq(self):
        self.pressure_comp_equil = Var(
                self._params.component_list,
                initialize=1e3,
                doc="Component partial pressures at equilibrium [Pa]")

        try:
            # Assume Langmuir isotherm for equilibrium
            def langmuir_isotherm(b, j):
                if j == "Zeo":
                    return Constraint.Skip
                else:
                    return (b.solid_state_ref.loading[j] *
                            (1+b.k_eq*b.pressure_comp_equil[j]) ==
                            b._params.loading_max[j]*b.k_eq *
                            b.pressure_comp_equil[j])
            self.langmuir_isotherm = Constraint(
                    self._params.component_list,
                    rule=langmuir_isotherm,
                    doc="Langmuir Isotherm")
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.pressure_comp_equil)
            self.del_component(self.langmuir_isotherm)
            raise

    def _k_eq(self):
        self.k_eq = Var(initialize=1e-3,
                        doc="Langmuir equilibrium coefficient [?]")

        try:
            # Temperature dependence of k_eq: deltaG = -RT*Ln(K_eq)
            def gibbs_relationship(b):
                return (b._params.gibbs_mol_rxn ==
                        -8.314*b.solid_state_ref.temperature*log(b.k_eq))
            self.gibbs_relationship = Constraint(
                    rule=gibbs_relationship,
                    doc="Relationship between Gibbs energy and "
                    "equilibrium coefficient")
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.k_eq)
            self.del_component(self.gibbs_relationship)
            raise

    def _rxn_rate(self):
        self.reaction_rate = Var(self._params.rate_reaction_idx,
                                 initialize=0,
                                 doc="Rate of reaction [?]")

        try:
            # Adsorption rate constraint
            def reaction_rate_eqn(b):
                return (b.reaction_rate["R1"] ==
                        b.k_rxn*(b.gas_state_ref.pressure *
                                 b.gas_state_ref.mole_frac_comp["N2"] -
                                 b.pressure_comp_equil["N2"]))
            self.reaction_rate_eqn = Constraint(
                    rule=reaction_rate_eqn,
                    doc="Reaction rate constraint")
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.reaction_rate)
            self.del_component(self.reaction_rate_eqn)
            raise

    def _k_rxn(self):
        self.k_rxn = Var(initialize=1e-3,
                         doc="Adsorption rate constant [?]")

        try:
            # Arrhenius equation
            def arrhenius_eqn(b):
                return (b.k_rxn == b._params.arrhenius *
                        exp(-b._params.energy_activation /
                            (8.314*b.solid_state_ref.temperature)))
            self.arrhenius_eqn = Constraint(
                    rule=arrhenius_eqn,
                    doc="Arrhenius equation constraint")
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.k_rxn)
            self.del_component(self.arrhenius_eqn)
            raise

    def get_reaction_rate_basis(b):
        return MaterialFlowBasis.mass

    def __getattr__(self, attr):
        """
        This method is used to avoid generating unnecessary property
        calculations in reaction blocks. __getattr__ is called whenever a
        property is called for, and if a propery does not exist, it looks for
        a method to create the required property, and any associated
        components.

        Create a property calculation if needed. Return an attrbute error if
        attr == 'domain' or starts with a _ . The error for _ prevents a
        recursion error if trying to get a function to create a property and
        that function doesn't exist.  Pyomo also ocasionally looks for things
        that start with _ and may not exist.  Pyomo also looks for the domain
        attribute, and it may not exist.
        This works by creating a property calculation by calling the "_"+attr
        function.

        A list of __getattr__ calls is maintained in self.__getattrcalls to
        check for recursive loops which maybe useful for debugging. This list
        is cleared after __getattr__ completes successfully.

        Args:
            attr: an attribute to create and return. Should be a property
                  component.
        """
        def clear_call_list(self, attr):
            """Local method for cleaning up call list when a call is handled.

                Args:
                    attr: attribute currently being handled
            """
            if self.__getattrcalls[-1] == attr:
                if len(self.__getattrcalls) <= 1:
                    del self.__getattrcalls
                else:
                    del self.__getattrcalls[-1]
            else:
                raise PropertyPackageError(
                        "{} Trying to remove call {} from __getattr__"
                        " call list, however this is not the most "
                        "recent call in the list ({}). This indicates"
                        " a bug in the __getattr__ calls. Please "
                        "contact the IDAES developers with this bug."
                        .format(self.name, attr, self.__getattrcalls[-1]))

        # Check that attr is not something we shouldn't touch
        if attr == "domain" or attr.startswith("_"):
            # Don't interfere with anything by getting attributes that are
            # none of my business
            raise PropertyPackageError(
                    '{} {} does not exist, but is a protected '
                    'attribute. Check the naming of your '
                    'components to avoid any reserved names'
                    .format(self.name, attr))

        if attr == "config":
            try:
                self._get_config_args()
                return self.config
            except:
                raise BurntToast("{} getattr method was triggered by a call "
                                 "to the config block, but _get_config_args "
                                 "failed. This should never happen.")

        # Check for recursive calls
        try:
            # Check to see if attr already appears in call list
            if attr in self.__getattrcalls:
                # If it does, indicates a recursive loop.
                if attr == self.__getattrcalls[-1]:
                    # attr method is calling itself
                    self.__getattrcalls.append(attr)
                    raise PropertyPackageError(
                                    '{} _{} made a recursive call to '
                                    'itself, indicating a potential '
                                    'recursive loop. This is generally '
                                    'caused by the {} method failing to '
                                    'create the {} component.'
                                    .format(self.name, attr, attr, attr))
                else:
                    self.__getattrcalls.append(attr)
                    raise PropertyPackageError(
                                    '{} a potential recursive loop has been '
                                    'detected whilst trying to construct {}. '
                                    'A method was called, but resulted in a '
                                    'subsequent call to itself, indicating a '
                                    'recursive loop. This may be caused by a '
                                    'method trying to access a component out '
                                    'of order for some reason (e.g. it is '
                                    'declared later in the same method). See '
                                    'the __getattrcalls object for a list of '
                                    'components called in the __getattr__ '
                                    'sequence.'
                                    .format(self.name, attr))
            # If not, add call to list
            self.__getattrcalls.append(attr)
        except AttributeError:
            # A list of calls if one does not exist, so create one
            self.__getattrcalls = [attr]

        # Get property information from get_supported_properties
        try:
            m = self.config.parameters.get_metadata().properties

            if m is None:
                raise PropertyPackageError(
                        '{} reaction package get_supported_properties'
                        ' method returned None when trying to create '
                        '{}. Please contact the developer of the '
                        'property package'.format(self.name, attr))
        except KeyError:
            # If attr not in get_supported_properties, assume package does not
            # support property
            clear_call_list(self, attr)
            raise PropertyNotSupportedError(
                    '{} {} is not supported by reaction package (property is '
                    'not listed in get_supported_properties).'
                    .format(self.name, attr, attr))

        # Get method name from get_supported_properties
        try:
            if m[attr]['method'] is None:
                # If method is none, property should be constructed
                # by property package, so raise PropertyPackageError
                clear_call_list(self, attr)
                raise PropertyPackageError(
                        '{} {} should be constructed automatically '
                        'by reaction package, but is not present. '
                        'This can be caused by methods being called '
                        'out of order.'.format(self.name, attr))
            elif m[attr]['method'] is False:
                # If method is False, package does not support property
                # Raise NotImplementedError
                clear_call_list(self, attr)
                raise PropertyNotSupportedError(
                        '{} {} is not supported by reaction package '
                        '(property method is listed as False in '
                        'get_supported_properties).'
                        .format(self.name, attr))
            elif isinstance(m[attr]['method'], str):
                # Try to get method name in from PropertyBlock object
                try:
                    f = getattr(self, m[attr]['method'])
                except AttributeError:
                    # If fails, method does not exist
                    clear_call_list(self, attr)
                    raise PropertyPackageError(
                            '{} {} get_supported_properties method '
                            'returned a name that does not correspond'
                            ' to any method in the reaction package. '
                            'Please contact the developer of the '
                            'reaction package.'.format(self.name, attr))
            else:
                # Otherwise method name is invalid
                clear_call_list(self, attr)
                raise PropertyPackageError(
                             '{} {} get_supported_properties method '
                             'returned invalid value for method name. '
                             'Please contact the developer of the '
                             'reaction package.'
                             .format(self.name, attr))
        except KeyError:
            # No method key - raise Exception
            # Need to use an AttributeError so Pyomo.DAE will handle this
            clear_call_list(self, attr)
            raise PropertyNotSupportedError(
                    '{} get_supported_properties method '
                    'does not contain a method for {}. '
                    'Please select a package which supports '
                    'the necessary properties for your process.'
                    .format(self.name, attr))

        # Call attribute if it is callable
        # If this fails, it should return a meaningful error.
        if callable(f):
            try:
                f()
            except Exception:
                # Clear call list and reraise error
                clear_call_list(self, attr)
                raise
        else:
            # If f is not callable, inform the user and clear call list
            clear_call_list(self, attr)
            raise PropertyPackageError(
                    '{} tried calling attribute {} in order to create '
                    'component {}. However the method is not callable.'
                    .format(self.name, f, attr))

        # Clear call list, and return
        comp = getattr(self, attr)
        clear_call_list(self, attr)
        return comp
