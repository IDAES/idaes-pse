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
Base class for control volumes
"""

# Import Python libraries
from enum import Enum

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
from idaes.core import (
    ProcessBlockData,
    MaterialFlowBasis,
    useDefault,
    declare_process_block_class,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
    DefaultBool,
)
from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    PropertyNotSupportedError,
)
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

__author__ = "Andrew Lee"


# Enumerate options for material balances
class MaterialBalanceType(Enum):
    useDefault = -1
    none = 0
    componentPhase = 1
    componentTotal = 2
    elementTotal = 3
    total = 4


# Enumerate options for energy balances
class EnergyBalanceType(Enum):
    useDefault = -1
    none = 0
    enthalpyPhase = 1
    enthalpyTotal = 2
    energyPhase = 3
    energyTotal = 4


# Enumerate options for momentum balances
class MomentumBalanceType(Enum):
    none = 0
    pressureTotal = 1
    pressurePhase = 2
    momentumTotal = 3
    momentumPhase = 4


# Enumerate options for flow direction
class FlowDirection(Enum):
    forward = 1
    backward = 2


# Set up example ConfigBlock that will work with ControlVolume autobuild method
CONFIG_Template = ProcessBlockData.CONFIG()
CONFIG_Template.declare(
    "dynamic",
    ConfigValue(
        default=useDefault,
        domain=DefaultBool,
        description="Dynamic model flag",
        doc="""Indicates whether this model will be dynamic,
**default** - useDefault.
**Valid values:** {
**useDefault** - get flag from parent,
**True** - set as a dynamic model,
**False** - set as a steady-state model}""",
    ),
)
CONFIG_Template.declare(
    "has_holdup",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Holdup construction flag",
        doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}""",
    ),
)
CONFIG_Template.declare(
    "material_balance_type",
    ConfigValue(
        default=MaterialBalanceType.componentPhase,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.componentPhase.
**Valid values:** {
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
    ),
)
CONFIG_Template.declare(
    "energy_balance_type",
    ConfigValue(
        default=EnergyBalanceType.enthalpyTotal,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single ethalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - ethalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
    ),
)
CONFIG_Template.declare(
    "momentum_balance_type",
    ConfigValue(
        default=MomentumBalanceType.pressureTotal,
        domain=In(MomentumBalanceType),
        description="Momentum balance construction flag",
        doc="""Indicates what type of momentum balance should be constructed,
**default** - MomentumBalanceType.pressureTotal.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
    ),
)
CONFIG_Template.declare(
    "has_rate_reactions",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Rate reaction construction flag",
        doc="""Indicates whether terms for rate controlled reactions should be
constructed,
**default** - False.
**Valid values:** {
**True** - include kinetic reaction terms,
**False** - exclude kinetic reaction terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_equilibrium_reactions",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Equilibrium reaction construction flag",
        doc="""Indicates whether terms for equilibrium controlled reactions
should be constructed,
**default** - False.
**Valid values:** {
**True** - include equilibrium reaction terms,
**False** - exclude equilibrium reaction terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_phase_equilibrium",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Phase equilibrium construction flag",
        doc="""Indicates whether terms for phase equilibrium should be
constructed,
**default** = False.
**Valid values:** {
**True** - include phase equilibrium terms
**False** - exclude phase equilibrium terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_mass_transfer",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Mass transfer term construction flag",
        doc="""Indicates whether terms for mass transfer should be constructed,
**default** - False.
**Valid values:** {
**True** - include mass transfer terms,
**False** - exclude mass transfer terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_heat_of_reaction",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Heat of reaction term construction flag",
        doc="""Indicates whether terms for heat of reaction should be constructed,
**default** - False.
**Valid values** {
**True** - include heat of reaction terms,
**False** - exclude heat of reaction terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_heat_transfer",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Heat transfer term construction flag",
        doc="""Indicates whether terms for heat transfer should be constructed,
**default** - False.
**Valid values:** {
**True** - include heat transfer terms,
**False** - exclude heat transfer terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_work_transfer",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Work transfer term construction flag",
        doc="""Indicates whether terms for work transfer should be constructed,
**default** - False.
**Valid values** {
**True** - include work transfer terms,
**False** - exclude work transfer terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_enthalpy_transfer",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Enthalpy transfer term construction flag",
        doc="""Indicates whether terms for enthalpy transfer due to mass trasnfer
should be constructed, **default** - False.
**Valid values** {
**True** - include enthalpy transfer terms,
**False** - exclude enthalpy transfer terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_pressure_change",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
    ),
)
CONFIG_Template.declare(
    "property_package",
    ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
    ),
)
CONFIG_Template.declare(
    "property_package_args",
    ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
    ),
)
CONFIG_Template.declare(
    "reaction_package",
    ConfigValue(
        default=None,
        domain=is_reaction_parameter_block,
        description="Reaction package to use for control volume",
        doc="""Reaction parameter object used to define reaction calculations,
**default** - None.
**Valid values:** {
**None** - no reaction package,
**ReactionParameterBlock** - a ReactionParameterBlock object.}""",
    ),
)
CONFIG_Template.declare(
    "reaction_package_args",
    ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing reaction packages",
        doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see reaction package for documentation.}""",
    ),
)


@declare_process_block_class(
    "ControlVolume",
    doc="This class is not usually used directly. "
    "Use ControlVolume0DBlock or ControlVolume1DBlock"
    " instead.",
)
class ControlVolumeBlockData(ProcessBlockData):
    """
    The ControlVolumeBlockData Class forms the base class for all IDAES
    ControlVolume models. The purpose of this class is to automate the tasks
    common to all control volume blockss and ensure that the necessary
    attributes of a control volume block are present.

    The most signfiicant role of the ControlVolumeBlockData class is to set up
    the construction arguments for the control volume block, automatically link
    to the time domain of the parent block, and to get the information about
    the property and reaction packages.
    """

    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=DefaultBool,
            default=useDefault,
            description="Dynamic model flag",
            doc="""Indicates whether this model will be dynamic,
**default** - useDefault.
**Valid values:** {
**useDefault** - get flag from parent,
**True** - set as a dynamic model,
**False** - set as a steady-state model}""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=useDefault,
            domain=DefaultBool,
            description="Holdup construction flag",
            doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these, **default** - None. **Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "reaction_package",
        ConfigValue(
            default=None,
            domain=is_reaction_parameter_block,
            description="Reaction package to use for control volume",
            doc="""Reaction parameter object used to define reaction calculations,
**default** - None.
**Valid values:** {
**None** - no reaction package,
**ReactionParameterBlock** - a ReactionParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "reaction_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing reaction packages",
            doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see reaction package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "auto_construct",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Argument indicating whether ControlVolume should "
            "automatically construct balance equations",
            doc="""If set to True, this argument will trigger the auto_construct
method which will attempt to construct a set of material, energy and momentum
balance equations based on the parent unit's config block. The parent unit must
have a config block which derives from CONFIG_Base,
**default** - False.
**Valid values:** {
**True** - use automatic construction,
**False** - do not use automatic construciton.}""",
        ),
    )

    def build(self):
        """
        General build method for Control Volumes blocks. This method calls a
        number of sub-methods which automate the construction of expected
        attributes of all ControlVolume blocks.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        super(ControlVolumeBlockData, self).build()

        # Setup dynamics flag and time domain
        self._setup_dynamics()

        # Get property package details
        self._get_property_package()

        # Get indexing sets
        self._get_indexing_sets()

        # Get reaction package details (as necessary)
        self._get_reaction_package()

        if self.config.auto_construct is True:
            self._auto_construct()

    def add_geometry(self, *args, **kwargs):
        """
        Method for defining the geometry of the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_geometry. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_material_balances(
        self, balance_type=MaterialBalanceType.useDefault, **kwargs
    ):
        """
        General method for adding material balances to a control volume.
        This method makes calls to specialised sub-methods for each type of
        material balance.

        Args:
            balance_type - MaterialBalanceType Enum indicating which type of
                    material balance should be constructed.
            has_rate_reactions - whether default generation terms for rate
                    reactions should be included in material balances
            has_equilibrium_reactions - whether generation terms should for
                    chemical equilibrium reactions should be included in
                    material balances
            has_phase_equilibrium - whether generation terms should for phase
                    equilibrium behaviour should be included in material
                    balances
            has_mass_transfer - whether generic mass transfer terms should be
                    included in material balances
            custom_molar_term - a Pyomo Expression representing custom terms to
                    be included in material balances on a molar basis.
            custom_mass_term - a Pyomo Expression representing custom terms to
                    be included in material balances on a mass basis.

        Returns:
            Constraint objects constructed by sub-method
        """
        # Check if balance_type is useDefault, and get default if necessary
        if balance_type == MaterialBalanceType.useDefault:
            try:
                blk = self._get_representative_property_block()
                balance_type = blk.default_material_balance_type()
            except NotImplementedError:
                raise ConfigurationError(
                    "{} property package has not implemented a "
                    "default_material_balance_type, thus cannot use "
                    "MaterialBalanceType.useDefault when constructing "
                    "material balances. Please contact the developer of "
                    "your property package to implement the necessary "
                    "default attributes.".format(self.name)
                )

        self._constructed_material_balance_type = balance_type
        if balance_type == MaterialBalanceType.none:
            mb = None
        elif balance_type == MaterialBalanceType.componentPhase:
            mb = self.add_phase_component_balances(**kwargs)
        elif balance_type == MaterialBalanceType.componentTotal:
            mb = self.add_total_component_balances(**kwargs)
        elif balance_type == MaterialBalanceType.elementTotal:
            mb = self.add_total_element_balances(**kwargs)
        elif balance_type == MaterialBalanceType.total:
            mb = self.add_total_material_balances(**kwargs)
        else:
            raise ConfigurationError(
                "{} invalid balance_type for add_material_balances."
                "Please contact the unit model developer with this bug.".format(
                    self.name
                )
            )

        return mb

    def add_energy_balances(self, balance_type=EnergyBalanceType.useDefault, **kwargs):
        """
        General method for adding energy balances to a control volume.
        This method makes calls to specialised sub-methods for each type of
        energy balance.

        Args:
            balance_type (EnergyBalanceType): Enum indicating which type of
                energy balance should be constructed.
            has_heat_of_reaction (bool): whether terms for heat of reaction
                should be included in energy balance
            has_heat_transfer (bool): whether generic heat transfer terms
                should be included in energy balances
            has_work_transfer (bool): whether generic mass transfer terms
                should be included in energy balances
            has_enthalpy_transfer (bool): whether generic enthalpy transfer
                terms should be included in energy balances
            custom_term (Expression): a Pyomo Expression representing custom
                terms to be included in energy balances

        Returns:
            Constraint objects constructed by sub-method
        """
        # Check if balance_type is useDefault, and get default if necessary
        if balance_type == EnergyBalanceType.useDefault:
            try:
                blk = self._get_representative_property_block()
                balance_type = blk.default_energy_balance_type()
            except NotImplementedError:
                raise ConfigurationError(
                    "{} property package has not implemented a "
                    "default_energy_balance_type, thus cannot use "
                    "EnergyBalanceType.useDefault when constructing "
                    "energy balances. Please contact the developer of "
                    "your property package to implement the necessary "
                    "default attributes.".format(self.name)
                )

        self._constructed_energy_balance_type = balance_type
        if balance_type == EnergyBalanceType.none:
            eb = None
        elif balance_type == EnergyBalanceType.enthalpyTotal:
            eb = self.add_total_enthalpy_balances(**kwargs)
        elif balance_type == EnergyBalanceType.enthalpyPhase:
            eb = self.add_phase_enthalpy_balances(**kwargs)
        elif balance_type == EnergyBalanceType.energyTotal:
            eb = self.add_total_energy_balances(**kwargs)
        elif balance_type == EnergyBalanceType.energyPhase:
            eb = self.add_phase_energy_balances(**kwargs)
        else:
            raise ConfigurationError(
                "{} invalid balance_type for add_energy_balances."
                "Please contact the unit model developer with this bug.".format(
                    self.name
                )
            )

        return eb

    def add_momentum_balances(
        self, balance_type=MomentumBalanceType.pressureTotal, **kwargs
    ):
        """
        General method for adding momentum balances to a control volume.
        This method makes calls to specialised sub-methods for each type of
        momentum balance.

        Args:
            balance_type (MomentumBalanceType): Enum indicating which type of
                momentum balance should be constructed. Default =
                MomentumBalanceType.pressureTotal.
            has_pressure_change (bool): whether default generation terms for
                pressure change should be included in momentum balances
            custom_term (Expression): a Pyomo Expression representing custom
                terms to be included in momentum balances

        Returns:
            Constraint objects constructed by sub-method
        """
        self._constructed_momentum_balance_type = balance_type
        if balance_type == MomentumBalanceType.none:
            mb = None
        elif balance_type == MomentumBalanceType.pressureTotal:
            mb = self.add_total_pressure_balances(**kwargs)
        elif balance_type == MomentumBalanceType.pressurePhase:
            mb = self.add_phase_pressure_balances(**kwargs)
        elif balance_type == MomentumBalanceType.momentumTotal:
            mb = self.add_total_momentum_balances(**kwargs)
        elif balance_type == MomentumBalanceType.momentumPhase:
            mb = self.add_phase_momentum_balances(**kwargs)
        else:
            raise ConfigurationError(
                "{} invalid balance_type for add_momentum_balances."
                "Please contact the unit model developer with this bug.".format(
                    self.name
                )
            )

        return mb

    def _auto_construct(self):
        """
        Placeholder _auto_construct method to ensure a useful exception is
        returned if auto_build is set to True but something breaks in the
        process. Derived ControlVolume classes should overload this.

        Args:
            None

        Returns:
            None
        """
        parent = self.parent_block()

        self.add_geometry()
        self.add_state_blocks()
        self.add_reaction_blocks()

        self.add_material_balances(
            material_balance_type=parent.config.material_balance_type,
            has_rate_reactions=parent.config.has_rate_reactions,
            has_equilibrium_reactions=parent.config.has_equilibrium_reactions,
            has_phase_equilibrium=parent.config.has_phase_equilibrium,
            has_mass_transfer=parent.config.has_mass_transfer,
        )

        self.add_energy_balances(
            energy_balance_type=parent.config.energy_balance_type,
            has_heat_of_reaction=parent.config.has_heat_of_reaction,
            has_heat_transfer=parent.config.has_heat_transfer,
            has_work_transfer=parent.config.has_work_transfer,
            has_enthalpy_transfer=parent.config.has_enthalpy_transfer,
        )

        self.add_momentum_balances(
            has_pressure_change=parent.config.has_pressure_change
        )

        try:
            self.apply_transformation()
        except AttributeError:
            pass

    # Add placeholder methods for adding property and reaction packages
    def add_state_blocks(self, *args, **kwargs):
        """
        Method for adding StateBlocks to the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_state_blocks. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_reaction_blocks(self, *args, **kwargs):
        """
        Method for adding ReactionBlocks to the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_reaction_blocks. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    # Add placeholder methods for all types of material, energy and momentum
    # balance equations which return NotImplementedErrors
    def add_phase_component_balances(self, *args, **kwargs):
        """
        Method for adding material balances indexed by phase and component to
        the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_phase_component_material_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_total_component_balances(self, *args, **kwargs):
        """
        Method for adding material balances indexed by component to
        the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_total_component_material_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_total_element_balances(self, *args, **kwargs):
        """
        Method for adding total elemental material balances indexed to
        the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_total_element_material_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_total_material_balances(self, *args, **kwargs):
        """
        Method for adding a total material balance to
        the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_total_material_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_phase_enthalpy_balances(self, *args, **kwargs):
        """
        Method for adding enthalpy balances indexed by phase to
        the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_phase_enthalpy_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_total_enthalpy_balances(self, *args, **kwargs):
        """
        Method for adding a total enthalpy balance to
        the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_total_enthalpy_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_phase_energy_balances(self, *args, **kwargs):
        """
        Method for adding energy balances (including kinetic energy) indexed by
        phase to the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_phase_energy_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_total_energy_balances(self, *args, **kwargs):
        """
        Method for adding a total energy balance (including kinetic energy)
        to the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_total_energy_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_phase_pressure_balances(self, *args, **kwargs):
        """
        Method for adding pressure balances indexed by
        phase to the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_phase_pressure_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_total_pressure_balances(self, *args, **kwargs):
        """
        Method for adding a total pressure balance to the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_total_pressure_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_phase_momentum_balances(self, *args, **kwargs):
        """
        Method for adding momentum balances indexed by phase to the control
        volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_phase_momentum_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_total_momentum_balances(self, *args, **kwargs):
        """
        Method for adding a total momentum balance to the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_total_momentum_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def _rxn_rate_conv(b, t, x, j, has_rate_reactions):
        """
        Method to determine conversion term for reaction rate terms in material
        balance equations. This method gets the basis of the material flow
        and reaction rate terms and determines the correct conversion factor.
        """
        # If rate reactions are not required, skip the rest and return 1
        if not has_rate_reactions:
            return 1

        if x is None:
            # 0D control volume
            flow_basis = b.properties_out[t].get_material_flow_basis()
            prop = b.properties_out[t]
            rxn_basis = b.reactions[t].get_reaction_rate_basis()
        else:
            # 1D control volume
            flow_basis = b.properties[t, x].get_material_flow_basis()
            prop = b.properties[t, x]
            rxn_basis = b.reactions[t, x].get_reaction_rate_basis()

        # Check for undefined basis
        if flow_basis == MaterialFlowBasis.other:
            raise ConfigurationError(
                "{} contains reaction terms, but the property package "
                "used an undefined basis (MaterialFlowBasis.other). "
                "Rate based reaction terms require the property "
                "package to define the basis of the material flow "
                "terms.".format(b.name)
            )
        if rxn_basis == MaterialFlowBasis.other:
            raise ConfigurationError(
                "{} contains reaction terms, but the reaction package "
                "used an undefined basis (MaterialFlowBasis.other). "
                "Rate based reaction terms require the reaction "
                "package to define the basis of the reaction rate "
                "terms.".format(b.name)
            )

        try:
            if flow_basis == rxn_basis:
                return 1
            elif (
                flow_basis == MaterialFlowBasis.mass
                and rxn_basis == MaterialFlowBasis.molar
            ):
                return prop.mw_comp[j]
            elif (
                flow_basis == MaterialFlowBasis.molar
                and rxn_basis == MaterialFlowBasis.mass
            ):
                return 1 / prop.mw_comp[j]
            else:
                raise BurntToast(
                    "{} encountered unrecognsied combination of bases "
                    "for reaction rate terms. Please contact the IDAES"
                    " developers with this bug.".format(b.name)
                )
        except AttributeError:
            raise PropertyNotSupportedError(
                "{} property package does not support "
                "molecular weight (mw), which is required for "
                "using property and reaction packages with "
                "different bases.".format(b.name)
            )

    def _get_representative_property_block(self):
        try:
            t_ref = self.flowsheet().time.first()
        except AttributeError:
            raise ConfigurationError(
                "{} control volume does not appear to be part of a "
                "flowsheet (could not find a time attribute).".format(self.name)
            )

        try:
            rep_blk = self.properties_out[t_ref]
        except AttributeError:
            try:
                d_ref = self.length_domain.first()
                rep_blk = self.properties[t_ref, d_ref]
            except AttributeError:
                raise BurntToast(
                    "{} Something went wrong when trying to find "
                    "a representative StateBlock. Please contact "
                    "the IDAES developers with this bug.".format(self.name)
                )
        return rep_blk
