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
Tool to interrogate IDAES flowsheets and list the physical properties
required to simualte it.
"""

# Import Pyomo libraries
from pyomo.environ import Set, Var

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock,
                        MaterialBalanceType,
                        EnergyBalanceType)
import idaes.logger as idaeslog

# Some more inforation about this module
__author__ = "Andrew Lee"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("InterrogatorParameterBlock")
class InterrogatorParameterData(PhysicalParameterBlock):
    """
    Interrogator Parameter Block Class

    This class contains the methods and attributes for recording and displaying
    the properties requried by the flowsheet.
    """
    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(InterrogatorParameterData, self).build()

        self.state_block_class = InterrogatorStateBlock

        # List of valid phases in property package
        # TODO : Allow users to define a phase list
        self.phase_list = Set(initialize=['Liq', 'Vap'])

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['A', 'B'])

        # Set up dict to record property calls
        self.required_properties = {}

        # Create dummy Vars to be used in constructing flowsheet
        self.dummy_var_unindexed = Var(initialize=42)
        self.dummy_var_phase = Var(self.phase_list, initialize=42)
        self.dummy_var_comp = Var(self.component_list, initialize=42)
        self.dummy_var_phase_comp = Var(self.phase_list,
                                        self.component_list,
                                        initialize=42)


class _InterrogatorStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """
    def initialize(blk, *args, **kwargs):
        '''
        Dummy initialization routine, This will raise an Exception if a user
        tries to initialize a model using the Interrogator Property Package
        and tell them that the model cannot be solved.
        '''
        raise Exception(
                "Models constructed using the Property Interrogator package "
                "cannot be used to solve a flowsheet. Please rebuild your "
                "flowsheet using a valid property package.")


@declare_process_block_class("InterrogatorStateBlock",
                             block_class=_InterrogatorStateBlock)
class InterrogatorStateBlockData(StateBlockData):
    """
    A dummy state block for interrogating flowsheets and recording physical
    properties called for during construction.
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(InterrogatorStateBlockData, self).build()

    # Define standard methods and log calls before returning dummy variable
    def get_material_flow_terms(self, p, j):
        self._log_call("material flow terms")
        return self._params.dummy_var_unindexed

    def get_enthalpy_flow_terms(self, p):
        self._log_call("enthalpy flow terms")
        return self._params.dummy_var_unindexed

    def get_material_density_terms(self, p, j):
        self._log_call("material density terms")
        return self._params.dummy_var_unindexed

    def get_energy_density_terms(self, p):
        self._log_call("energy density terms")
        return self._params.dummy_var_unindexed

    # Set default values for required attributes so construction doesn't fail
    def default_material_balance_type(self):
        return MaterialBalanceType.componentPhase

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def define_state_vars(b):
        return {"flow_vol": b.flow_vol,
                "conc_mol_comp": b.conc_mol_comp,
                "temperature": b.temperature,
                "pressure": b.pressure}

    def define_display_vars(b):
        raise Exception(
                "Models constructed using the Property Interrogator package "
                "should not be used for report methods.")

    def get_material_flow_basis(b):
        return MaterialFlowBasis.molar

    def _log_call(self, prop):
        # Log call for prop in required_properties
        prop_dict = self._params.required_properties
        try:
            # Check if current block is already listed for property
            if self.name not in prop_dict[prop]:
                prop_dict[prop].append(self.name)
        except KeyError:
            prop_dict[prop] = [self.name]

    def __getattr__(self, prop):
        """
        Overload getattr to log each call for an unknown attribute, assuming
        these are all properties.
        Then, return a dummy variable with the correct indexing set.
        """
        # Log call
        self._log_call(prop)

        # Return dummy var
