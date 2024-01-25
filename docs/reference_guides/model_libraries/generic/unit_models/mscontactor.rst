Multi-Stream Contactor
======================

.. module:: idaes.models.unit_models.mscontactor

The Multi-Stream Contactor is a general purpose model for unit operations involving mass and energy transfer between multiple immiscible streams, such as membrane filtration systems and solvent extraction processes. The Multi-Stream Contactor provides a general framework for writing the necessary material, energy and momentum balances for each stream and includes terms for inter-stream transfer phenomena. The model also supports modeling these systems as a series of finite elements (either representing a series of well-mixed stages or a one-dimensional variation using a 1st order finite difference approximation). Finally, the model also supports the presence of side feeds/draws for each stream.

When adding a multi-stream contactor to a flowsheet, users can define the options they wish to use for the unit and each stream as shown below. A list of all available configuration options is shown later in the class documentation.

.. testcode::

  import pyomo.environ as pyo  # Pyomo environment
  from idaes.core import FlowsheetBlock, FlowDirection
  from idaes.models.unit_models import MSContactor
  from idaes.models.properties import iapws95

  # Create an empty flowsheet and steam property parameter block.
  model = pyo.ConcreteModel()
  model.fs = FlowsheetBlock(dynamic=False)
  model.fs.properties = iapws95.Iapws95ParameterBlock()

  # Add a multi-stream contactor model to the flowsheet.
  model.fs.contactor = MSContactor(
      number_of_finite_elements=2,
      streams={
          "stream1": {
              "property_package": model.fs.properties
          },
          "stream2": {
              "property_package": model.fs.properties,
              "flow_direction": FlowDirection.backward,
          },
      },
  )

Degrees of Freedom
------------------

As a general purpose model, the degrees of freedom of the multi-stream contactor models depend on the options chosen by the user. The potential degrees of freedom are:

* states for feed blocks for each stream,
* material transfer terms (time points :math:`\times` finite elements :math:`\times` interacting streams :math:`\times` common components),
* energy transfer terms if included (time points :math:`\times` finite elements :math:`\times` interacting streams)
* pressure change terms if included (time points :math:`\times` finite elements :math:`\times` streams with pressure change)
* reaction extent terms for rate based and heterogeneous reactions if included (time points :math:`\times` finite elements :math:`\times` number of reactions)

Model Structure
---------------

Due to the custom nature of multi-stream contactors, this model does not make use of control volumes. Instead, a set of StateBlocks (named using the name given in the ``streams`` configuration dictionary) are created for each stream indexed by time and the set of finite elements, with an additional StateBlock for the feed state (named using the stream name appended with ``_inlet_state``) indexed only by time (unless ``has_feed`` is set to ``False`` for that stream). For streams with side streams (feed or draw), an additional set of indexed StateBlocks (named using the stream name appended with ``_side_stream_state``) is created for the side states which are indexed by time and the set of side states for that stream.

If reactions are required for a given stream, a set of indexed ReactionBlocks (named using the stream name appended with ``_reactions``) are created indexed by time and the set of finite elements. The ``MSContactor`` model also supports the concept of heterogeneous reactions which involve species across multiple streams; in these cases users can declare a single "heterogeneous reaction block" for the system in which case a Heterogeneous ReactionBlock will be created for each finite element (see later for more details on defining heterogeneous reaction blocks).

**Note:** due to the custom nature of multi-stream contactors, the ``MSContactor`` model **does not** define a constraint for the extent of reaction as the basis for this is not well defined (e.g. volume can be defined either by the volume of each stream/phase in a finite element or the total volume in each element. Users are required to add the necessary constraints appropriate to their system.

All other variables and constraints are written at the unit model level.

Variables
---------

The multi-stream contactor creates the following variables. Here ``t`` indicates the time domain and ``x`` indicates finite element.

============================= ============================================== ============================================================================================================= =========================================
Variable                      Name                                           Description                                                                                                   Notes
============================= ============================================== ============================================================================================================= =========================================
:math:`M_{t,x,s1,s2,j}`       material_transfer_term                         Material transfer term for component ``j`` between stream ``s1`` and ``s2`` at ``x`` and ``t``
:math:`E_{t,x,s1,s2}`         energy_transfer_term                           Energy transfer term between stream ``s1`` and ``s2`` at ``x`` and ``t``                                      Only if energy balances included
:math:`Q_{t,x,s}`             stream + "_heat"                               External heat transfer into stream ``s`` at ``x`` and ``t``                                                   Only if ``has_heat_transfer`` for stream
:math:`\Delta P_{t,x,s}`      stream + "_deltaP"                             Pressure change in stream ``s`` at ``x`` and ``t``                                                            Only if ``has_pressure_change`` for stream
:math:`G_{rate,t,x,s,p,j}`    stream + "_rate_reaction_generation"           Generation of component ``j`` in phase ``p`` due to rate reactions in stream ``s`` at ``x`` ``t``             Only if rate reactions present for stream
:math:`G_{equil,t,x,s,p,j}`   stream + "_equilibrium_reaction_generation"    Generation of component ``j`` in phase ``p`` due to equilibrium reactions in stream ``s`` at ``x`` and ``t``  Only if equilibrium reactions present for stream
:math:`G_{inher,t,x,s,p,j}`   stream + "_inherent_reaction_generation"       Generation of component ``j`` in phase ``p`` due to inherent reactions in stream ``s`` at ``x`` and ``t``     Only if inherent reactions present for stream
:math:`G_{hetero,t,x,s,p,j}`  stream + "_heterogeneous_reaction_generation"  Generation of component ``j`` in phase ``p`` due to heterogeneous reactions in stream ``s`` at ``x`` ``t``    Only if heterogeneous reactions present
:math:`X_{rate,t,x,s,r}`      stream + "_rate_reaction_extent"               Extent of rate reaction ``r`` in stream ``s`` at ``x`` and ``t``                                              Only if rate reactions present for stream
:math:`X_{equil,t,x,s,r}`     stream + "_equilibrium_reaction_extent"        Extent of equilibrium reaction ``r`` in stream ``s`` at ``x`` and ``t``                                       Only if equilibrium reactions present for stream
:math:`X_{inher,t,x,s,r}`     stream + "_inherent_reaction_extent"           Extent of inherent reaction ``r`` in stream ``s`` at ``x`` and ``t``                                          Only if inherent reactions present for stream
:math:`X_{hetero,t,x,r}`      heterogeneous_reaction_extent                  Extent of heterogeneous reaction ``r`` at ``x`` and ``t``                                                     Only if heterogeneous reactions present
:math:`V_x`                   volume                                         Total volume of element ``x``                                                                                 Only if ``has_holdup``
:math:`f_{t,x,s}`             volume_frac_stream                             Volume fraction of stream ``s`` in element ``x`` at time ``t``                                                Only if ``has_holdup``
:math:`\phi_{t,x,s,p}`        stream + "_phase_fraction"                     Volume fraction of phase ``p`` in stream ``s`` in element ``x`` at time ``t``                                 Only if ``has_holdup``
:math:`N_{t,x,s,p,j}`         stream + "_material_holdup"                    Holdup of component ``j`` in phase ``p`` for stream ``s`` at ``x`` and ``t``                                  Only if ``has_holdup``
:math:`dN/dt_{t,x,s,p,j}`     stream + "_material_accumulation"              Accumulation of component ``j`` in phase ``p`` for stream ``s`` at ``x`` and ``t``                            Only if ``dynamic``
:math:`U_{t,x,s,p}`           stream + "_energy_holdup"                      Holdup of energy in phase ``p`` for stream ``s`` at ``x`` and ``t``                                           Only if ``has_holdup``
:math:`dU/dt_{t,x,s,p}`       stream + "_energy_accumulation"                Accumulation of energy in phase ``p`` for stream ``s`` at ``x`` and ``t``                                     Only if ``dynamic``
============================= ============================================== ============================================================================================================= =========================================

Constraints
-----------

In all cases, the multi-stage contactor model writes a set of material balances for each stream in the model. For component ``j`` in stream ``s`` the following constraint, named ``stream + "_material_balance"``, is written for all finite elements ``x``:

.. math:: dN/dt_{t,x,s,p,j} = \sum_p{F_{t,x-,s,p,j}} - \sum_p{F_{t,x,s,p,j}} + \left[ \sum_p{F_{side,t,x,s,p,j}} \right] + \sum_o{M_{t,x,s,o,j}} + \left[ \sum_p{G_{rate,t,x,s,p,j}} + \sum_p{G_{equil,t,x,s,p,j}} + \sum_p{G_{inher,t,x,s,p,j}} + \sum_p{G_{hetero,t,x,s,p,j}} \right]

where ``F`` is the material flow term, ``x-`` represents the previous finite element (``x-1`` in the case of co-current flow and ``x+1`` in the case of counter-current flow), ``F_side`` is the material flow term for a side stream (if present) and ``o`` represents all other streams in the model (for cases where ``s`` is the second index (i.e., M_{t,x,o,s,j}) the term is multiplied by -1). The reaction generation terms are only included if the appropriate reaction type is supported by the reaction or property package for the stream.

For systems including rate reactions, the following constraint, names ``stream + "_rate_reaction_constraint"``, is written to relate the generation of component ``j`` in phase ``p`` to the extent of each rate reaction as shown below where :math:`\alpha_{r,p,j}` is the stoichiometric coefficient for component ``j`` in phase ``p`` for reaction ``r``.

.. math:: G_{rate,t,x,s,p,j} = \sum_r{\alpha_{rate_r,p,j} \times X_{rate,t,x,s,r}}

Equivalent constraints are written for equilibrium, inherent and heterogeneous reactions as necessary.

For streams including energy balances (``has_energy_balance = True``) the following constraint (named ``stream + "_energy_balance"``) is written at each finite element:

.. math:: \sum_p{dU/dt_{t,x,s,p}} = \sum_p{H_{t,x-,s,p}} - \sum_p{H_{t,x,s,p}} + \biggl[ \sum_p{H_{side,t,x,s,p}} \biggr] + \sum_o{E_{t,x,s,o}} + \biggl[ Q_{t,x,s} \biggr] + \biggl[ \sum_{rate}{\Delta H_{rxn,r} \times X_{rate,t,x,s,r}} + \sum_{equil}{\Delta H_{rxn,r} \times X_{equil,t,x,s,r}} + \sum_{inher}{\Delta H_{rxn,r} \times X_{inher,t,x,s,r}} \biggr]

where ``H`` represent enthalpy flow terms and :math:`\Delta H_{rxn}` represents heat of reaction. The heat of reaction terms are only included if a reaction package is provided for the stream AND the configuration option ``has_heat_of_reaction = True`` is set for the stream. **Note** heterogeneous reactions **do not** support heat of reaction terms as it is uncertain which stream/phase the heat should be added too.

For streams including pressure balances (``has_pressure_balance = True``) the following constraint (named ``stream + "_pressure_balance"``) is written at each finite element:

.. math:: 0 = P_{t,x-,s} - P_{t,x,s} + \biggl[ \Delta P_{t,x,s} \biggr]

where ``P`` represents pressure. For streams with side streams, the following pressure equality constraint (named ``stream + "_side_stream_pressure_balance"``) is also written:

.. math:: P_{t,x,s} = P_{side,t,x,s}

If ``has_holdup`` is true, the following additional constraints are included to calculate holdup terms. First, ``sum_volume_frac`` constrains the sum of all volume fractions to be 1.

.. math:: 1 = \sum_s{f_{t,x,s}}

Additionally, constraints are written for the sum of phase fractions in each stream (named ``stream + "_sum_phase_fractions"``):

.. math:: 1 = \sum_p{\phi_{t,x,s,p}}

The material holdup is defined by the following constraint (named ``stream + "_material_holdup_constraint"``):

.. math:: N_{t,x,s,p,j} = V \times f_{t,x,s} \phi_{t,x,s,p} \times \C_{t,x,s,p,j}

where :math:`C_{t,x,s,p,j}` is the concentration of component ``j`` in phase ``p`` for stream ``s`` at ``x`` and ``t``.

The energy holdup is defined by the following constraint (named ``stream + "_energy_holdup_constraint"``):

.. math:: U_{t,x,s,p} = V*f_{t,x,s} \times \phi_{t,x,s,p} \times \u_{t,x,s,p}

where :math:`u_{t,x,s,p}` is the internal energy density of phase ``p`` for stream ``s`` at ``x`` and ``t``.

Initialization
--------------

.. autoclass:: MSContactorInitializer
   :members: initialization_routine

MSContactor Class
-----------------

.. autoclass:: MSContactor
  :members:
  
Stream Configuration Options
----------------------------

============================= ============================= ======== ==============================================================================================
Argument                      Type                          Default  Description
============================= ============================= ======== ==============================================================================================
property_package              PropertyParameter Block       None     Property package associated with stream
property_package_args         dict                          None     Configuration arguments for State Blocks
reaction_package              Reaction Parameter Block      None     Reaction package associated with stream
reaction_package_args         dict                          None     Configuration arguments for Reaction Blocks
flow_direction                FlowDirection Enum            forward  Direction of flow for stream
has_feed                      bool                          True     Whether stream has a feed Port and inlet state, or if all flow is provided via mass transfer.
has_rate_reactions            bool                          False    Whether rate-based reactions occur in stream.
has_equilibrium_reactions     bool                          False    Whether equilibrium-based reactions occur in stream.
has_energy_balance            bool                          True     Whether to include energy balance for stream.
has_heat_transfer             bool                          False    Whether to include external heat transfer terms in energy balance for stream.
has_heat_of_reaction          bool                          False    Whether heat of reaction terms should be included in energy balance for stream.
has_pressure_balance          bool                          True     Whether to include pressure balance for stream.
has_pressure_change           bool                          False    Whether to include :math:`\Delta P` terms in pressure balance for stream.
side_streams                  list                          None     Finite elements at which a side stream should be included.
heterogeneous_reactions       heterogeneous reaction block  None     Heterogeneous reaction package to use for system
heterogeneous_reactions_args  dict                          None     Configuration arguments for heterogeneous reaction blocks
============================= ============================= ======== ==============================================================================================

Heterogeneous Reaction Blocks
-----------------------------

Heterogeneous reaction blocks are a new feature in IDAES, and are currently still in beta development. Due to this, there is no base class for heterogeneous reaction blocks yet as the API is not yet finalized.

Currently, the requirements for are similar to those for ReactionBlocks and are demonstrated in the example outline below:

.. code-block::

    from pyomo.environ import Constraint, Set, units, Var
    from pyomo.common.config import ConfigValue

    from idaes.core import (
        declare_process_block_class,
        ProcessBlockData,
        ProcessBlock,
    )
    from idaes.core.base import property_meta
    from idaes.core.util.misc import add_object_reference

    # -----------------------------------------------------------------------------
    # Heterogeneous Reaction Parameter Block
    @declare_process_block_class("MyHeterogeneousReactionParameters")
    class MyHeterogeneousReactionParametersData(ProcessBlockData, property_meta.HasPropertyClassMetadata):
        def build(self):
            super().build()

            self._reaction_block_class = MyHeterogeneousReactionsBlock

            self.reaction_idx = Set(
                initialize=[
                    "reaction1",
                    "reaction2",
                    ...
                ]
            )

            self.reaction_stoichiometry = {
                ("reaction1", "phase1", "component1"): stoichiometric_coefficient,
                ("reaction1", "phase1", "component2"): stoichiometric_coefficient,
                ...
            }

            # Add any global parameters here, such as activation energies and Arrhenius constants

        # Define base units for reactions
        @classmethod
        def define_metadata(cls, obj):
            obj.add_default_units(
                {
                    "time": units.hour,
                    "length": units.m,
                    "mass": units.kg,
                    "amount": units.mol,
                    "temperature": units.K,
                }
            )

        # The next few lines are boilerplate - you should be able to just copy these
        @property
        def reaction_block_class(self):
            return self._reaction_block_class

        def build_reaction_block(self, *args, **kwargs):
            """
            Methods to construct a ReactionBlock associated with this
            ReactionParameterBlock. This will automatically set the parameters
            construction argument for the ReactionBlock.

            Returns:
                ReactionBlock

            """
            default = kwargs.pop("default", {})
            initialize = kwargs.pop("initialize", {})

            if initialize == {}:
                default["parameters"] = self
            else:
                for i in initialize.keys():
                    initialize[i]["parameters"] = self

            return self.reaction_block_class(  # pylint: disable=not-callable
                *args, **kwargs, **default, initialize=initialize
            )


    # Define the heterogenous ReactionBlock
    @declare_process_block_class(
        "MyHeterogeneousReactionsBlock", block_class=ProcessBlock
    )
    class MyHeterogeneousReactionsData(ProcessBlockData):
        # Create Class ConfigBlock - this needs to be here
        CONFIG = ProcessBlockData.CONFIG()
        CONFIG.declare(
            "parameters",
            ConfigValue(
                description="""A reference to an instance of the Heterogeneous Reaction Parameter
        Block associated with this property package.""",
            ),
        )
        # Add any additional configuration arguments you want

        def build(self):
            super().build()

            # This creates an easy link back to the parameters
            add_object_reference(self, "_params", self.config.parameters)

            # Need to define a reaction rate variable - make sure units are correct
            self.reaction_rate = Var(
                self.params.reaction_idx,
                initialize=0,
                units=#units
                )

            # Define rule for calculating reaction rate
            def rule_reaction_rate_eq(b, r):
                return b.reaction_rate[r] == #rate expression

            self.reaction_rate_eq = Constraint(self.params.reaction_idx, rule=rule_reaction_rate_eq)

        # This is boilerplate
        @property
        def params(self):
            return self._params
