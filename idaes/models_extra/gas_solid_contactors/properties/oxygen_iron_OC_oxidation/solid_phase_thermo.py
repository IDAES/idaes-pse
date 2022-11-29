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
This package provides the necessary constraints for solid phase properties of
an iron-based oxygen carrier
Components - Fe2O3, Fe3O4, Al2O3

Equations written in this model were primarily derived from:
National Institute of Standards and Technology, NIST Chemistry WebBook,
https://webbook.nist.gov/chemistry/ (accessed March 10, 2018).

"""

# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    Expression,
    Param,
    Reals,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialFlowBasis,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    Component,
    SolidPhase,
)
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
from idaes.core.util.misc import add_object_reference
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_unfixed_variables_in_activated_equalities,
)
import idaes.logger as idaeslog
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver

# Some more information about this module
__author__ = "Chinedu Okoli"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("SolidPhaseParameterBlock")
class PhysicalParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class

    Contains parameters and indexing sets associated with properties for
    methane CLC.
    """

    def build(self):
        """
        Callable method for Block construction.
        """
        super(PhysicalParameterData, self).build()

        self._state_block_class = SolidPhaseStateBlock

        # Create Phase object
        self.Sol = SolidPhase()

        # Create Component objects
        self.Fe2O3 = Component()
        self.Fe3O4 = Component()
        self.Al2O3 = Component()

        # -------------------------------------------------------------------------
        """ Pure solid component properties"""

        # Mol. weights of solid components - units = kg/mol. ref: NIST webbook
        mw_comp_dict = {"Fe2O3": 0.15969, "Fe3O4": 0.231533, "Al2O3": 0.10196}
        self.mw_comp = Param(
            self.component_list,
            mutable=False,
            initialize=mw_comp_dict,
            doc="Molecular weights of solid components [kg/mol]",
            units=pyunits.kg / pyunits.mol,
        )

        # Skeletal density of solid components - units = kg/m3. ref: NIST
        dens_mass_comp_skeletal_dict = {"Fe2O3": 5250, "Fe3O4": 5000, "Al2O3": 3987}
        self.dens_mass_comp_skeletal = Param(
            self.component_list,
            mutable=False,
            initialize=dens_mass_comp_skeletal_dict,
            doc="Skeletal density of solid components [kg/m3]",
            units=pyunits.kg / pyunits.m**3,
        )

        # Ideal gas spec. heat capacity parameters(Shomate) of
        # components - ref: NIST webbook. Shomate equations from NIST.
        # Parameters A-E are used for cp calcs while A-H are used for enthalpy
        # calc.
        # 1e3*cp_comp = A + B*T + C*T^2 + D*T^3 + E/(T^2)
        # where T = Temperature (K)/1000, and cp_comp = (kJ/mol.K)
        # H_comp = H - H(298.15) = A*T + B*T^2/2 + C*T^3/3 +
        # D*T^4/4 - E/T + F - H where T = Temp (K)/1000 and H_comp = (kJ/mol)
        cp_param_dict = {
            ("Al2O3", 1): 102.4290,
            ("Al2O3", 2): 38.74980,
            ("Al2O3", 3): -15.91090,
            ("Al2O3", 4): 2.628181,
            ("Al2O3", 5): -3.007551,
            ("Al2O3", 6): -1717.930,
            ("Al2O3", 7): 146.9970,
            ("Al2O3", 8): -1675.690,
            ("Fe3O4", 1): 200.8320000,
            ("Fe3O4", 2): 1.586435e-7,
            ("Fe3O4", 3): -6.661682e-8,
            ("Fe3O4", 4): 9.452452e-9,
            ("Fe3O4", 5): 3.18602e-8,
            ("Fe3O4", 6): -1174.1350000,
            ("Fe3O4", 7): 388.0790000,
            ("Fe3O4", 8): -1120.8940000,
            ("Fe2O3", 1): 110.9362000,
            ("Fe2O3", 2): 32.0471400,
            ("Fe2O3", 3): -9.1923330,
            ("Fe2O3", 4): 0.9015060,
            ("Fe2O3", 5): 5.4336770,
            ("Fe2O3", 6): -843.1471000,
            ("Fe2O3", 7): 228.3548000,
            ("Fe2O3", 8): -825.5032000,
        }
        self.cp_param_1 = Param(
            self.component_list,
            mutable=False,
            initialize={k: v for (k, j), v in cp_param_dict.items() if j == 1},
            doc="Shomate equation heat capacity coeff 1",
            units=pyunits.J / pyunits.mol / pyunits.K,
        )
        self.cp_param_2 = Param(
            self.component_list,
            mutable=False,
            initialize={k: v for (k, j), v in cp_param_dict.items() if j == 2},
            doc="Shomate equation heat capacity coeff 2",
            units=pyunits.J / pyunits.mol / pyunits.K / pyunits.kK,
        )
        self.cp_param_3 = Param(
            self.component_list,
            mutable=False,
            initialize={k: v for (k, j), v in cp_param_dict.items() if j == 3},
            doc="Shomate equation heat capacity coeff 3",
            units=pyunits.J / pyunits.mol / pyunits.K / pyunits.kK**2,
        )
        self.cp_param_4 = Param(
            self.component_list,
            mutable=False,
            initialize={k: v for (k, j), v in cp_param_dict.items() if j == 4},
            doc="Shomate equation heat capacity coeff 4",
            units=pyunits.J / pyunits.mol / pyunits.K / pyunits.kK**3,
        )
        self.cp_param_5 = Param(
            self.component_list,
            mutable=False,
            initialize={k: v for (k, j), v in cp_param_dict.items() if j == 5},
            doc="Shomate equation heat capacity coeff 5",
            units=pyunits.J / pyunits.mol / pyunits.K * pyunits.kK**2,
        )
        self.cp_param_6 = Param(
            self.component_list,
            mutable=False,
            initialize={k: v for (k, j), v in cp_param_dict.items() if j == 6},
            doc="Shomate equation heat capacity coeff 6",
            units=pyunits.kJ / pyunits.mol,
        )
        self.cp_param_7 = Param(
            self.component_list,
            mutable=False,
            initialize={k: v for (k, j), v in cp_param_dict.items() if j == 7},
            doc="Shomate equation heat capacity coeff 7",
            units=pyunits.J / pyunits.mol / pyunits.K,
        )
        self.cp_param_8 = Param(
            self.component_list,
            mutable=False,
            initialize={k: v for (k, j), v in cp_param_dict.items() if j == 8},
            doc="Shomate equation heat capacity coeff 8",
            units=pyunits.kJ / pyunits.mol,
        )

        # Std. heat of formation of comp. - units = J/(mol comp) - ref: NIST
        enth_mol_form_comp_dict = {
            "Fe2O3": -825.5032e3,
            "Fe3O4": -1120.894e3,
            "Al2O3": -1675.690e3,
        }
        self.enth_mol_form_comp = Param(
            self.component_list,
            mutable=False,
            initialize=enth_mol_form_comp_dict,
            doc="Component molar heats of formation [J/mol]",
            units=pyunits.J / pyunits.mol,
        )

        # -------------------------------------------------------------------------
        """ Mixed solid properties"""
        # These are setup as fixed vars to allow for parameter estimation

        # Particle size
        self.particle_dia = Var(
            domain=Reals,
            initialize=1.5e-3,
            doc="Diameter of solid particles [m]",
            units=pyunits.m,
        )
        self.particle_dia.fix()

        # TODO -provide reference
        # Minimum fluidization velocity - EPAT value used for Davidson model
        self.velocity_mf = Var(
            domain=Reals,
            initialize=0.039624,
            doc="Velocity at minimum fluidization [m/s]",
            units=pyunits.m / pyunits.s,
        )
        self.velocity_mf.fix()

        # Minimum fluidization voidage - educated guess as rough
        # estimate from ergun equation results (0.4) are suspicious
        self.voidage_mf = Var(
            domain=Reals,
            initialize=0.45,
            doc="Voidage at minimum fluidization [-]",
            units=pyunits.dimensionless,
        )
        self.voidage_mf.fix()

        # Voidage of the bed
        self.voidage = Var(
            domain=Reals,
            initialize=0.35,
            doc="Voidage [-]",
            units=pyunits.dimensionless,
        )
        self.voidage.fix()

        # Particle thermal conductivity
        self.therm_cond_sol = Var(
            domain=Reals,
            initialize=12.3,
            doc="Thermal conductivity of solid particles [J/m.K.s]",
            units=pyunits.J / pyunits.m / pyunits.K / pyunits.s,
        )
        self.therm_cond_sol.fix()

        # Set default scaling for state variables
        self.set_default_scaling("flow_mass", 1e-3)
        self.set_default_scaling("particle_porosity", 1e2)
        self.set_default_scaling("temperature", 1e-2)
        for comp in self.component_list:
            self.set_default_scaling("mass_frac_comp", 1e1, index=comp)

        # Set default scaling for thermophysical and transport properties
        self.set_default_scaling("enth_mass", 1e-6)
        self.set_default_scaling("enth_mol_comp", 1e-6)
        self.set_default_scaling("cp_mol_comp", 1e-6)
        self.set_default_scaling("cp_mass", 1e-6)
        self.set_default_scaling("dens_mass_particle", 1e-2)
        self.set_default_scaling("dens_mass_skeletal", 1e-2)

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_mass": {"method": None},
                "particle_porosity": {"method": None},
                "temperature": {"method": None},
                "mass_frac_comp": {"method": None},
                "dens_mass_skeletal": {"method": "_dens_mass_skeletal"},
                "dens_mass_particle": {"method": "_dens_mass_particle"},
                "cp_mol_comp": {"method": "_cp_mol_comp"},
                "cp_mass": {"method": "_cp_mass"},
                "enth_mass": {"method": "_enth_mass"},
                "enth_mol_comp": {"method": "_enth_mol_comp"},
            }
        )

        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )


class _SolidPhaseStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to State Blocks as a
    whole, rather than individual elements of indexed State Blocks.
    """

    def initialize(
        blk,
        state_args=None,
        hold_state=False,
        state_vars_fixed=False,
        outlvl=idaeslog.NOTSET,
        solver="ipopt",
        optarg={"tol": 1e-8},
    ):
        """
        Initialization routine for property package.
        Keyword Arguments:
            state_args : Dictionary with initial guesses for the state vars
                         chosen. Note that if this method is triggered
                         through the control volume, and if initial guesses
                         were not provided at the unit model level, the
                         control volume passes the inlet values as initial
                         guess.
                         Keys for the state_args dictionary are:
                         flow_mass, temperature, and mass_frac_comp
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None)
            solver : str indicating whcih solver to use during
                     initialization (default = "ipopt")
            hold_state : flag indicating whether the initialization routine
                         should unfix any state variables fixed during
                         initialization (default=False).
                         - True - states varaibles are not unfixed, and
                                 a dict of returned containing flags for
                                 which states were fixed during
                                 initialization.
                        - False - state variables are unfixed after
                                 initialization by calling the
                                 relase_state method
        Returns:
            If hold_states is True, returns a dict containing flags for
            which states were fixed during initialization.
        """
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="properties")

        init_log.info_high("Starting initialization")

        # Deactivate the constraints specific for outlet block i.e.
        # when defined state is False
        for k in blk.keys():
            if blk[k].config.defined_state is False:
                blk[k].sum_component_eqn.deactivate()

        # Fix state variables if not already fixed
        if state_vars_fixed is False:
            flags = fix_state_vars(blk, state_args)
        else:
            # Check when the state vars are fixed already result in dof 0
            for k in blk.keys():
                if degrees_of_freedom(blk[k]) != 0:
                    raise Exception(
                        "State vars fixed but degrees of freedom "
                        "for state block is not zero during "
                        "initialization."
                    )

        # ---------------------------------------------------------------------
        # Initialise values
        for k in blk.keys():
            if hasattr(blk[k], "density_skeletal_constraint"):
                calculate_variable_from_constraint(
                    blk[k].dens_mass_skeletal, blk[k].density_skeletal_constraint
                )

            if hasattr(blk[k], "mixture_heat_capacity_eqn"):
                calculate_variable_from_constraint(
                    blk[k].cp_mass, blk[k].mixture_heat_capacity_eqn
                )

            if hasattr(blk[k], "mixture_enthalpy_eqn"):
                calculate_variable_from_constraint(
                    blk[k].enth_mass, blk[k].mixture_enthalpy_eqn
                )

            for j in blk[k]._params.component_list:

                if hasattr(blk[k], "cp_shomate_eqn"):
                    calculate_variable_from_constraint(
                        blk[k].cp_mol_comp[j], blk[k].cp_shomate_eqn[j]
                    )

                if hasattr(blk[k], "enthalpy_shomate_eqn"):
                    calculate_variable_from_constraint(
                        blk[k].enth_mol_comp[j], blk[k].enthalpy_shomate_eqn[j]
                    )

        # Solve property block if non-empty
        free_vars = 0
        for k in blk.keys():
            free_vars += number_unfixed_variables_in_activated_equalities(blk[k])

        if free_vars > 0:
            # Create solver
            opt = get_solver(solver, optarg)
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
        else:
            res = ""
        init_log.info_high(
            "Initialization complete {}.".format(idaeslog.condition(res))
        )

        # ---------------------------------------------------------------------
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags)

    def release_state(blk, flags, outlvl=0):
        """
        Method to relase state variables fixed during initialization.
        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        """
        if flags is None:
            return

        # Unfix state variables
        revert_state_vars(blk, flags)

        # Activate state variable related constraints
        for k in blk.keys():
            if blk[k].config.defined_state is False:
                blk[k].sum_component_eqn.activate()

        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        init_log.info_high("States released.")


@declare_process_block_class("SolidPhaseStateBlock", block_class=_SolidPhaseStateBlock)
class SolidPhaseStateBlockData(StateBlockData):
    """
    Property package for gas phase properties of methane combustion in CLC FR
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(SolidPhaseStateBlockData, self).build()

        # Object reference for molecular weight if needed by CV1D
        # Molecular weights
        add_object_reference(self, "mw_comp", self.config.parameters.mw_comp)

        self._make_state_vars()

    def _make_state_vars(self):
        """List the necessary state variable objects."""

        # create units object to get default units from the param block
        units_meta = self._params.get_metadata().derived_units

        self.flow_mass = Var(
            initialize=1.0,
            domain=Reals,
            doc="Component mass flowrate",
            units=units_meta.FLOW_MASS,
        )
        self.particle_porosity = Var(
            domain=Reals,
            initialize=0.27,
            doc="Porosity of oxygen carrier [-]",
            units=pyunits.dimensionless,
        )
        self.mass_frac_comp = Var(
            self._params.component_list,
            initialize=1 / len(self._params.component_list),
            doc="State component mass fractions [-]",
            units=pyunits.dimensionless,
        )
        self.temperature = Var(
            initialize=298.15,
            domain=Reals,
            doc="State temperature",
            units=units_meta.TEMPERATURE,
        )

        # Create standard constraints
        # Sum mass fractions if not inlet block
        if self.config.defined_state is False:

            def sum_component_eqn(b):
                return 1 == sum(b.mass_frac_comp[j] for j in b._params.component_list)

            self.sum_component_eqn = Constraint(rule=sum_component_eqn)

    def _dens_mass_skeletal(self):
        # Skeletal density of OC solid particles
        units_meta = self._params.get_metadata().derived_units
        self.dens_mass_skeletal = Var(
            domain=Reals,
            initialize=3251.75,
            doc="Skeletal density of OC",
            units=units_meta.DENSITY_MASS,
        )

        def density_skeletal_constraint(b):
            return (
                b.dens_mass_skeletal
                * sum(
                    b.mass_frac_comp[j] / b._params.dens_mass_comp_skeletal[j]
                    for j in b._params.component_list
                )
                == 1
            )

        try:
            # Try to build constraint
            self.density_skeletal_constraint = Constraint(
                rule=density_skeletal_constraint
            )
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.dens_mass_skeletal)
            self.del_component(self.density_skeletal_constraint)
            raise

    def _dens_mass_particle(self):
        # Particle density of OC (includes the OC pores)
        units_meta = self._params.get_metadata().derived_units
        self.dens_mass_particle = Var(
            domain=Reals,
            initialize=3251.75,
            doc="Particle density of oxygen carrier",
            units=units_meta.DENSITY_MASS,
        )

        def density_particle_constraint(b):
            return (
                b.dens_mass_particle == (1 - b.particle_porosity) * b.dens_mass_skeletal
            )

        try:
            # Try to build constraint
            self.density_particle_constraint = Constraint(
                rule=density_particle_constraint
            )
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.dens_mass_particle)
            self.del_component(self.density_particle_constraint)
            raise

    def _cp_mol_comp(self):
        # Pure component solid heat capacities
        units_meta = self._params.get_metadata().derived_units
        units_cp_mol = units_meta.HEAT_CAPACITY_MOLE
        self.cp_mol_comp = Var(
            self._params.component_list,
            domain=Reals,
            initialize=1.0,
            doc="Pure component solid heat capacities",
            units=units_cp_mol,
        )

        def pure_component_cp_mol(b, j):
            t = pyunits.convert(b.temperature, to_units=pyunits.kK)
            return b.cp_mol_comp[j] == pyunits.convert(
                (
                    b._params.cp_param_1[j]
                    + b._params.cp_param_2[j] * t
                    + b._params.cp_param_3[j] * t**2
                    + b._params.cp_param_4[j] * t**3
                    + b._params.cp_param_5[j] / (t**2)
                ),
                to_units=units_cp_mol,
            )

        try:
            # Try to build constraint
            self.cp_shomate_eqn = Constraint(
                self._params.component_list, rule=pure_component_cp_mol
            )
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol_comp)
            self.del_component(self.cp_shomate_eqn)
            raise

    def _cp_mass(self):
        # Mixture heat capacities
        units_meta = self._params.get_metadata().derived_units
        units_cp_mass = units_meta.HEAT_CAPACITY_MASS
        self.cp_mass = Var(
            domain=Reals,
            initialize=1.0,
            doc="Mixture heat capacity, mass-basis",
            units=units_cp_mass,
        )

        def cp_mass(b):
            return b.cp_mass == sum(
                b.cp_mol_comp[j] * b.mass_frac_comp[j] * (1 / b._params.mw_comp[j])
                for j in b._params.component_list
            )

        try:
            # Try to build constraint
            self.mixture_heat_capacity_eqn = Constraint(rule=cp_mass)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mass)
            self.del_component(self.mixture_heat_capacity_eqn)
            raise

    def _enth_mol_comp(self):
        # Pure component vapour enthalpies
        units_meta = self._params.get_metadata().derived_units
        units_enth_mol = units_meta.ENERGY_MOLE
        self.enth_mol_comp = Var(
            self._params.component_list,
            domain=Reals,
            initialize=1.0,
            doc="Pure component enthalpies",
            units=units_enth_mol,
        )

        def pure_comp_enthalpy(b, j):
            t = pyunits.convert(b.temperature, to_units=pyunits.kK)
            return b.enth_mol_comp[j] == pyunits.convert(
                # parameters 1-5 are defined in J
                b._params.cp_param_1[j] * t
                + b._params.cp_param_2[j] * (t**2) / 2
                + b._params.cp_param_3[j] * (t**3) / 3
                + b._params.cp_param_4[j] * (t**4) / 4
                - b._params.cp_param_5[j] / (t),
                to_units=units_enth_mol,
            ) + pyunits.convert(
                # parameters 6 and 8 are defined in kJ, and must be added
                # after converting to the enthalpy units set
                b._params.cp_param_6[j] - b._params.cp_param_8[j],
                to_units=units_enth_mol,
            )

        try:
            # Try to build constraint
            self.enthalpy_shomate_eqn = Constraint(
                self._params.component_list, rule=pure_comp_enthalpy
            )
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.enth_mol_comp)
            self.del_component(self.enthalpy_shomate_eqn)
            raise

    def _enth_mass(self):
        # Mixture mass enthalpy
        units_meta = self._params.get_metadata().derived_units
        units_enth_mass = units_meta.ENERGY_MASS
        self.enth_mass = Var(
            domain=Reals,
            initialize=0.0,
            doc="Mixture specific enthalpy",
            units=units_enth_mass,
        )
        try:
            # Try to build constraint
            self.mixture_enthalpy_eqn = Constraint(
                expr=(
                    self.enth_mass
                    == sum(
                        self.mass_frac_comp[j]
                        * self.enth_mol_comp[j]
                        * (1 / self._params.mw_comp[j])
                        for j in self._params.component_list
                    )
                )
            )
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.enth_mass)
            self.del_component(self.mixture_enthalpy_eqn)
            raise

    def get_material_flow_terms(self, p, j):
        # return b.flow_mass*b.mass_frac_comp[j]
        if not self.is_property_constructed("material_flow_terms"):
            try:

                def rule_material_flow_terms(b, j):
                    return b.flow_mass * b.mass_frac_comp[j]

                self.material_flow_terms = Expression(
                    self._params.component_list, rule=rule_material_flow_terms
                )
            except AttributeError:
                self.del_component(self.material_flow_terms)
        return self.material_flow_terms[j]

    def get_enthalpy_flow_terms(self, p):
        if not self.is_property_constructed("enthalpy_flow_terms"):
            try:

                def rule_enthalpy_flow_terms(b):
                    return self.enth_mass * self.flow_mass

                self.enthalpy_flow_terms = Expression(rule=rule_enthalpy_flow_terms)
            except AttributeError:
                self.del_component(self.enthalpy_flow_terms)
        return self.enthalpy_flow_terms

    def get_material_density_terms(self, p, j):
        # return b.dens_mass_particle * b.mass_frac_comp[j]
        if not self.is_property_constructed("material_density_terms"):
            try:

                def rule_material_density_terms(b, j):
                    return b.dens_mass_particle * b.mass_frac_comp[j]

                self.material_density_terms = Expression(
                    self._params.component_list, rule=rule_material_density_terms
                )
            except AttributeError:
                self.del_component(self.material_density_terms)
        return self.material_density_terms[j]

    def get_energy_density_terms(self, p):
        if not self.is_property_constructed("energy_density_terms"):
            try:

                def rule_energy_density_terms(b):
                    return b.dens_mass_particle * b.enth_mass

                self.energy_density_terms = Expression(rule=rule_energy_density_terms)
            except AttributeError:
                self.del_component(self.energy_density_terms)
        return self.energy_density_terms

    def define_state_vars(b):
        return {
            "flow_mass": b.flow_mass,
            "particle_porosity": b.particle_porosity,
            "temperature": b.temperature,
            "mass_frac_comp": b.mass_frac_comp,
        }

    def get_material_flow_basis(b):
        return MaterialFlowBasis.mass

    def model_check(blk):
        """
        Model checks for property block
        """
        # Check temperature bounds
        if value(blk.temperature) < blk.temperature.lb:
            _log.error("{} Temperature set below lower bound.".format(blk.name))
        if value(blk.temperature) > blk.temperature.ub:
            _log.error("{} Temperature set above upper bound.".format(blk.name))

    def default_material_balance_type(blk):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(blk):
        return EnergyBalanceType.enthalpyTotal

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # scale some variables
        # nothing here

        # scale some constraints
        if self.is_property_constructed("material_flow_terms"):
            for i, c in self.material_flow_terms.items():
                sf1 = iscale.get_scaling_factor(self.mass_frac_comp[i])
                sf2 = iscale.get_scaling_factor(self.flow_mass)
                iscale.set_scaling_factor(c, sf1 * sf2)

        if self.is_property_constructed("material_density_terms"):
            for i, c in self.material_density_terms.items():
                sf1 = iscale.get_scaling_factor(self.mass_frac_comp[i])
                sf2 = iscale.get_scaling_factor(self.dens_mass_particle)
                iscale.set_scaling_factor(c, sf1 * sf2)

        if self.is_property_constructed("energy_density_terms"):
            for i, c in self.energy_density_terms.items():
                sf1 = iscale.get_scaling_factor(self.enth_mass)
                sf2 = iscale.get_scaling_factor(self.dens_mass_particle)
                iscale.set_scaling_factor(c, sf1 * sf2)

        if self.is_property_constructed("enthalpy_flow_terms"):
            for i, c in self.enthalpy_flow_terms.items():
                sf1 = iscale.get_scaling_factor(self.enth_mass)
                sf2 = iscale.get_scaling_factor(self.flow_mass)
                iscale.set_scaling_factor(c, sf1 * sf2)

        # Scale some constraints
        if self.is_property_constructed("sum_component_eqn"):
            iscale.constraint_scaling_transform(
                self.sum_component_eqn,
                iscale.get_scaling_factor(self.mass_frac_comp["Fe2O3"]),
                overwrite=False,
            )
        if self.is_property_constructed("density_particle_constraint"):
            iscale.constraint_scaling_transform(
                self.density_particle_constraint,
                iscale.get_scaling_factor(self.dens_mass_particle),
                overwrite=False,
            )
        if self.is_property_constructed("density_skeletal_constraint"):
            iscale.constraint_scaling_transform(
                self.density_skeletal_constraint,
                iscale.get_scaling_factor(self.dens_mass_skeletal),
                overwrite=False,
            )
        if self.is_property_constructed("cp_shomate_eqn"):
            for i, c in self.cp_shomate_eqn.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.cp_mol_comp[i]), overwrite=False
                )
        if self.is_property_constructed("mixture_heat_capacity_eqn"):
            iscale.constraint_scaling_transform(
                self.mixture_heat_capacity_eqn,
                iscale.get_scaling_factor(self.cp_mass),
                overwrite=False,
            )
        if self.is_property_constructed("enthalpy_shomate_eqn"):
            for i, c in self.enthalpy_shomate_eqn.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.enth_mol_comp[i]), overwrite=False
                )
        if self.is_property_constructed("mixture_enthalpy_eqn"):
            iscale.constraint_scaling_transform(
                self.mixture_enthalpy_eqn,
                iscale.get_scaling_factor(self.enth_mass),
                overwrite=False,
            )
