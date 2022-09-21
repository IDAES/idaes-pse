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
from pyomo.common.config import ConfigValue
from pyomo.core.base.constraint import Constraint
from pyomo.core.base.var import Var
from pyomo.core.base.set import Set
from pyomo.core.base.block import Block
from pyomo.core.base.units_container import units as pyunits
from pyomo.core.base.reference import Reference
from pyomo.network.port import Port
from pyomo.network.arc import Arc

from idaes.core import declare_process_block_class, StateBlock, UnitModelBlockData
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.config import (
    is_physical_parameter_block,
)

"""
A unit model for a "node" that connects several pipelines,
possibly with its own supply and demand streams.
"""


@declare_process_block_class("PipelineNode")
class PipelineNodeData(UnitModelBlockData):

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "property_package",
        ConfigValue(default=None, domain=is_physical_parameter_block),
    )
    CONFIG.declare(
        "n_inlet_pipelines",
        ConfigValue(default=1, domain=int),
    )
    CONFIG.declare(
        "n_outlet_pipelines",
        ConfigValue(default=1, domain=int),
    )
    CONFIG.declare(
        "n_supplies",
        ConfigValue(default=0, domain=int),
    )
    CONFIG.declare(
        "n_demands",
        ConfigValue(default=0, domain=int),
    )

    def build(self):
        super(PipelineNodeData, self).build()

        # self.config is the ConfigBlock "instantiated" from self.CONFIG
        # in ProcessBlockData.
        config = self.config
        time = self.flowsheet().time
        property_package = config.property_package

        if len(property_package.phase_list) != 1:
            raise ValueError(
                "%s can only be constructed with a "
                "single-phase property package.\n"
                "Got phases %s."
                % (self.__class__, [p for p in property_package.phase_list])
            )
        self.phase = next(iter(property_package.phase_list))
        if self.phase != "Vap":
            raise ValueError(
                '%s can only be constructed with a single phase, "Vap".'
                "Got phase %s." % (self.__class__, self.phase)
            )

        #
        # Make sure property package has the properties that this unit
        # model requires.
        #
        property_dict = property_package.get_metadata().properties
        if "pressure" not in property_dict:
            raise ValueError(
                "Property package supplied to pipeline must have a property "
                "for 'pressure', which was not found in %s." % type(property_package)
            )

        # Build a state block for this node. Note that the flow_mol attribute
        # of this state block should never be used. Should we delete it?
        # We should probably add a constraint that makes sure it is equal
        # to flow "through" this node.
        state_config = {"defined_state": True}
        self.state = property_package.build_state_block(time, **state_config)

        self.add_inlets()
        self.add_outlets()
        self.add_supplies()
        self.add_demands()
        self.add_flow_balance_con()
        self.add_total_flow_con()
        self.add_component_mixing_con()
        self.add_enthalpy_mixing_con()

        self.n_inlet_pipelines = 0
        self.n_outlet_pipelines = 0
        self._inlet_pipelines = {}
        self._outlet_pipelines = {}

    def inlet_pipelines(self):
        return self._inlet_pipelines

    def outlet_pipelines(self):
        return self._outlet_pipelines

    def get_inlet_pipeline(self, i):
        return self._inlet_pipelines[i]

    def get_outlet_pipeline(self, i):
        return self._outlet_pipelines[i]

    def get_port_name(self):
        return "port"

    def _get_port_and_references(self, name, block, doc=None):
        """
        This is a copy of UnitModel.add_port, except it returns components
        instead of adding them to self.
        """
        # Validate block object
        if not isinstance(block, StateBlock):
            raise ConfigurationError(
                "{} block object provided to add_port "
                "method is not an instance of a "
                "StateBlock object. IDAES port objects "
                "should only be associated with "
                "StateBlocks.".format(self.name)
            )

        # Create empty Port
        p = Port(doc=doc)
        # setattr(blk, name, p)

        # Get dict of Port members and names
        member_list = block[self.flowsheet().time.first()].define_port_members()

        # Create References for port members
        references = []
        for s in member_list:
            if not member_list[s].is_indexed():
                slicer = block[:].component(member_list[s].local_name)
            else:
                slicer = block[:].component(member_list[s].local_name)[...]

            r = Reference(slicer)
            # setattr(blk, "_"+s+"_"+name+"_ref", r)
            ref_name = "_" + s + "_" + name + "_ref"
            references.append((ref_name, r))

            # Add Reference to Port
            # Should this be deferred as well?
            p.add(r, s)

        return p, references

    def get_temperature_eq_con(self, state1, state2):
        time = self.flowsheet().time

        def temperature_eq_rule(b, t):
            return state1[t].temperature == state2[t].temperature

        return Constraint(time, rule=temperature_eq_rule)

    def get_pressure_eq_con(self, state1, state2):
        time = self.flowsheet().time

        def pressure_eq_rule(b, t):
            return state1[t].pressure == state2[t].pressure

        return Constraint(time, rule=pressure_eq_rule)

    def get_mole_frac_comp_eq_con(self, state1, state2):
        time = self.flowsheet().time
        comp_list = self.config.property_package.component_list

        def mole_frac_comp_eq_rule(b, t, j):
            return state1[t].mole_frac_comp[j] == state2[t].mole_frac_comp[j]

        return Constraint(time, comp_list, rule=mole_frac_comp_eq_rule)

    def get_port_block_rule(self, outlet=False):
        properties = self.config.property_package
        time = self.flowsheet().time
        state_config = {"defined_state": True}
        port_name = self.get_port_name()

        def block_rule(b, i):
            # Each inlet/outlet needs its own state block so it can have its
            # own flow rate, at the very least.
            b.state = properties.build_state_block(time, **state_config)
            port, refs = self._get_port_and_references(port_name, b.state)
            for ref_name, ref in refs:
                b.add_component(ref_name, ref)
            b.port = port
            b.has_pipeline = False

            # Add equality constraints between the node states and inlet/outlet
            # states for pressure and temperature.
            b.pressure_eq = self.get_pressure_eq_con(self.state, b.state)

            if outlet:
                # For inlets, mole fractions and temperatures are already
                # defined by mixing rules.
                b.mole_frac_comp_eq = self.get_mole_frac_comp_eq_con(
                    self.state, b.state
                )
                b.temperature_eq = self.get_temperature_eq_con(self.state, b.state)

        return block_rule

    def add_inlets(self):
        n_inlets = self.config.n_inlet_pipelines
        # TODO: Should these set indices be something other than integers?
        # For now, the same indices could be valid in inlets, outlets,
        # supplies, and demands, which could lead to hard-to-debug mistakes.
        self.inlet_set = Set(initialize=list(range(n_inlets)), dimen=1)
        rule = self.get_port_block_rule()
        self.inlets = Block(self.inlet_set, rule=rule)

    def add_outlets(self):
        n_outlets = self.config.n_outlet_pipelines
        self.outlet_set = Set(initialize=list(range(n_outlets)), dimen=1)
        rule = self.get_port_block_rule(outlet=True)
        self.outlets = Block(self.outlet_set, rule=rule)

    def get_supply_block_rule(self):
        """
        Add a block for each supply, which contains a state block,
        allowing us to specify different conditions for each supply.

        """
        time = self.flowsheet().time
        properties = self.config.property_package
        state_config = {"defined_state": True}

        def block_rule(b, i):
            # NOTE: Adding a state block here leaves us with temperature
            # and pressure variables that are unused. Strictly, pressure
            # should be set equal to the node's pressure, and temperature
            # should be used to calculate the node's temperature via an
            # energy balance.
            # We don't link pressure to avoid adding extra variables
            # (the node pressure may be used instead) and we don't have an
            # energy balance, so we leave temperature unconnected.
            b.state = properties.build_state_block(time, **state_config)
            # Add a reference here so we can access flow in the same way
            # as on demand blocks.
            b.flow_mol = Reference(b.state[:].flow_mol)

            # Add equation to link pressure of this supply
            # to that of the node.
            b.isobaric_eq = self.get_pressure_eq_con(self.state, b.state)

        return block_rule

    def get_demand_block_rule(self):
        """
        Add a block for each demand, which only contains a variable
        for flow rate. Other (intensive) variables have the same value
        as this node's state block.

        """
        # Should these demand blocks have references to the node's intensive
        # state variables? This would probably be convenient.
        time = self.flowsheet().time

        def block_rule(b, i):
            b.flow_mol = Var(time, initialize=100.0, units=pyunits.kmol / pyunits.hr)

        return block_rule

    def add_supplies(self):
        """
        Add set and indexed block for supplies.
        """
        n_supplies = self.config.n_supplies
        self.supply_set = Set(initialize=list(range(n_supplies)), dimen=1)
        rule = self.get_supply_block_rule()
        self.supplies = Block(self.supply_set, rule=rule)

    def add_demands(self):
        """
        Add set and indexed block for demands.
        """
        n_demands = self.config.n_demands
        self.demand_set = Set(initialize=list(range(n_demands)), dimen=1)
        # Should I give demands a state block? This seems redundant.
        rule = self.get_demand_block_rule()
        self.demands = Block(self.demand_set, rule=rule)

    def add_flow_balance_con(self):
        """
        Adds a total flow rate balance stating that inlet flow must equal
        outlet flow.

        """
        # Should this be a balance on component flow rates instead?
        time = self.flowsheet().time

        def flow_balance_rule(b, t):
            return sum(self.supplies[:].flow_mol[t]) + sum(
                self.inlets[:].state[t].flow_mol
            ) == sum(self.demands[:].flow_mol[t]) + sum(
                self.outlets[:].state[t].flow_mol
            )

        self.flow_balance = Constraint(time, rule=flow_balance_rule)

    def add_total_flow_con(self):
        """
        Calculates the flow rate through this node as the sum of flow rates
        from supplies and inlet pipelines.

        """
        time = self.flowsheet().time

        def total_flow_rule(b, t):
            return self.state[t].flow_mol == sum(self.supplies[:].flow_mol[t]) + sum(
                self.inlets[:].state[t].flow_mol
            )

        self.total_flow_eq = Constraint(time, rule=total_flow_rule)

    def add_component_mixing_con(self):
        """
        Adds a component mixing equation that calculates mole fractions
        of the node from mole fractions and flow rates of the supplies
        and inlet pipelines.

        """
        time = self.flowsheet().time
        component_list = self.config.property_package.component_list

        def component_mixing_rule(b, t, j):
            return (
                sum(self.supplies[i].state[t].flow_mol_comp[j] for i in self.supply_set)
                + sum(self.inlets[i].state[t].flow_mol_comp[j] for i in self.inlet_set)
                == self.state[t].flow_mol_comp[j]
            )

        self.component_mixing_eq = Constraint(
            time, component_list, rule=component_mixing_rule
        )

    def add_enthalpy_mixing_con(self):
        """
        Adds an equation that calculates temperature of the node from
        temperatures, heat capacities, and flow rates of the supplies
        and inlet pipelines.

        """
        time = self.flowsheet().time
        # We enforce that phase_list contains only "Vap"
        p = next(iter(self.config.property_package.phase_list))

        def enthalpy_mixing_rule(b, t):
            supplies = [self.supplies[i].state[t] for i in self.supply_set]
            inlets = [self.inlets[i].state[t] for i in self.inlet_set]
            return sum(supply.get_enthalpy_flow_terms(p) for supply in supplies) + sum(
                inlet.get_enthalpy_flow_terms(p) for inlet in inlets
            ) == self.state[t].get_enthalpy_flow_terms(p)

        self.enthalpy_mixing_eq = Constraint(time, rule=enthalpy_mixing_rule)

    def add_pipeline_to_inlet(self, pipeline, idx=None):
        """
        Creates an Arc between the outlet port of a pipeline and one
        of this node's inlet ports.

        """
        pipeline_port = pipeline.outlet_port
        if idx is None:
            idx = self.n_inlet_pipelines
        if idx not in self.inlet_set:
            # Either bad index given or no more inlets
            raise ValueError()
        if self.inlets[idx].has_pipeline:
            # Don't want to add two pipelines to the same inlet
            raise RuntimeError()
        self.inlets[idx].arc = Arc(ports=(pipeline_port, self.inlets[idx].port))
        self.inlets[idx].has_pipeline = True
        self._inlet_pipelines[idx] = pipeline
        self.n_inlet_pipelines += 1

    def add_pipeline_to_outlet(self, pipeline, idx=None):
        """
        Adds an Arc between one of this node's outlet ports
        and the inlet port of a pipeline.

        """
        pipeline_port = pipeline.inlet_port
        if idx is None:
            idx = self.n_outlet_pipelines
        if idx not in self.outlet_set:
            # Either bad index given or no more outlets
            raise ValueError()
        if self.outlets[idx].has_pipeline:
            # Don't want to add two pipelines to the same outlet
            raise RuntimeError()
        self.outlets[idx].arc = Arc(ports=(pipeline_port, self.outlets[idx].port))
        self.outlets[idx].has_pipeline = True
        self._outlet_pipelines[idx] = pipeline
        self.n_outlet_pipelines += 1
