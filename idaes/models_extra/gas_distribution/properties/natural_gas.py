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
Natural gas property package with a single (pseudo) component for use
in isothermal unit models.

Data sources:
    [1] Stochastic Optimal Control Model for Natural Gas Network
        Operations. V. Zavala, 2014, Comp. Chem. Eng.

"""

from pyomo.core.base.units_container import units as pyunits
from pyomo.core.base.var import Var
from pyomo.core.base.constraint import Constraint
from pyomo.core.base.expression import Expression
from pyomo.core.base.param import Param
from pyomo.core.expr.current import sqrt

# Import IDAES cores
from idaes.core import (
    PhysicalParameterBlock,
)
import idaes.logger as idaeslog

from idaes.core import (
    declare_process_block_class,
    VaporPhase,
    Component,
    MaterialFlowBasis,
    StateBlock,
    StateBlockData,
)
from idaes.core.util.constants import Constants

# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("NaturalGasParameterBlock")
class NaturalGasParameterBlockData(PhysicalParameterBlock):
    """ """

    def build(self):
        super(NaturalGasParameterBlockData, self).build()
        self._state_block_class = NaturalGasStateBlock
        self.Vap = VaporPhase()

        self.dens_nominal = Param(
            initialize=0.72,
            units=pyunits.kg / pyunits.m**3,
            doc="Density of the gas a standard temperature and pressure",
            # Used to convert mass flow rate between kg/hr and SCM/hr
        )

        nat_gas_data = {"mw": (18.0, pyunits.kg / pyunits.kmol)}
        nat_gas_config = {"parameter_data": nat_gas_data}
        self.natural_gas = Component(**nat_gas_config)
        ng = self.natural_gas

        # Set up dictionaries we will rely on to compute heat capacities.
        kJkmolK = pyunits.kJ / pyunits.kmol / pyunits.K
        kJkmolK2 = pyunits.kJ / pyunits.kmol / pyunits.K**2
        kJkmolK3 = pyunits.kJ / pyunits.kmol / pyunits.K**3
        kJkmolK4 = pyunits.kJ / pyunits.kmol / pyunits.K**4
        ng.cp_mol_coef_A = Param(initialize=2.34 * 18.0, units=kJkmolK)
        ng.cp_mol_coef_B = Param(initialize=0.0, units=kJkmolK2)
        ng.cp_mol_coef_C = Param(initialize=0.0, units=kJkmolK3)
        ng.cp_mol_coef_D = Param(initialize=0.0, units=kJkmolK4)

        ng.cv_mol_coef_A = Param(initialize=1.85 * 18.0, units=kJkmolK)
        ng.cv_mol_coef_B = Param(initialize=0.0, units=kJkmolK2)
        ng.cv_mol_coef_C = Param(initialize=0.0, units=kJkmolK3)
        ng.cv_mol_coef_D = Param(initialize=0.0, units=kJkmolK4)

        self.temperature_ref = Param(initialize=298.15, units=pyunits.K)

    @classmethod
    def define_metadata(cls, obj):
        kghr = pyunits.kg / pyunits.hr
        nondim = pyunits.dimensionless
        kgkmolK = pyunits.kg / pyunits.kmol / pyunits.K
        kmolm3 = pyunits.kmol / pyunits.m**3
        kgm3 = pyunits.kg / pyunits.m**3
        mps = pyunits.m / pyunits.s
        kmolhr = pyunits.kmol / pyunits.hr
        obj.add_properties(
            {
                # What do the units in this dict get used for? What if they're
                # different than the units we define for each variable in the
                # StateBlockData methods?
                "flow_mol": {"method": None, "units": kmolhr},
                "pressure": {"method": None, "units": pyunits.bar},
                "temperature": {"method": None, "units": pyunits.K},
                "mole_frac_comp": {"method": None, "units": pyunits.dimensionless},
                "flow_mol_comp": {"method": "_flow_mol_comp", "units": kmolhr},
                "mw": {"method": "_mw", "units": pyunits.kg / pyunits.kmol},
                "flow_mass": {"method": "_flow_mass", "units": kghr},
                "cp_mol_comp": {"method": "_cp_mol_comp", "units": kgkmolK},
                "cp_mol": {"method": "_cp_mol", "units": kgkmolK},
                "cp_mass": {"method": "_cp_mass", "units": kgkmolK},
                "cv_mol_comp": {"method": "_cv_mol_comp", "units": kgkmolK},
                "cv_mol": {"method": "_cv_mol", "units": kgkmolK},
                "cv_mass": {"method": "_cv_mass", "units": kgkmolK},
                "heat_capacity_ratio": {
                    "method": "_heat_capacity_ratio",
                    "units": nondim,
                },
                "heat_capacity_ratio_phase": {
                    "method": "_heat_capacity_ratio_phase",
                    "units": nondim,
                },
                "compressibility": {"method": "_compressibility", "units": nondim},
                "dens_mol": {"method": "_dens_mol", "units": kmolm3},
                "dens_mol_comp": {"method": "_dens_mol_comp", "units": kmolm3},
                "dens_mass": {"method": "_dens_mass", "units": kgm3},
                "speed_of_sound": {"method": "_speed_of_sound", "units": mps},
            }
        )
        # NOTE: We do not implement enthalpy as we are not yet using this
        # property package with a non-isothermal pipeline.

        obj.add_default_units(
            {
                "time": pyunits.hr,
                "length": pyunits.km,
                "mass": pyunits.kg,
                "amount": pyunits.kmol,
                "temperature": pyunits.K,
                # I would like to specify some reasonable units for area,
                # energy, and pressure, but this is not supported in IDAES.
                # "area": pyunits.m**2,
                # "energy": pyunits.kJ,
                # "pressure": pyunits.bar,
            }
        )


class NaturalGasStateBlock(StateBlock):
    # This is confusing.
    # (a) Not sure why this class is necessary when I don't want
    #     to attach any methods.
    # (b) Now is NaturalGasStateBlock defined twice? Once here and once
    #     by NaturalGasStateBlock? This appears to be the case. This
    #     class gets overridden...
    # declare_process_block_class on its own is confusing...
    pass


@declare_process_block_class(
    "NaturalGasStateBlock",
    block_class=NaturalGasStateBlock,
)
class NaturalGasStateBlockData(StateBlockData):
    def build(self):
        super().build()

        # TODO: Initialize to 10 (1e6 SCM)/day
        self.flow_mol = Var(
            initialize=300.0,
            doc="Molar flow rate",
            units=pyunits.kmol / pyunits.hr,
        )
        self.pressure = Var(
            initialize=50.0,
            doc="Gas pressure",
            units=pyunits.bar,
        )
        self.temperature = Var(
            initialize=298.15,
            doc="Gas temperature",
            units=pyunits.K,
        )
        component_list = self.config.parameters.component_list
        self.mole_frac_comp = Var(
            component_list,
            initialize=1.0 / len(component_list),
            doc="Component mole fractions within the gas",
            units=pyunits.dimensionless,
        )

        if not self.config.defined_state:
            # For a single-phase property package, should we just always
            # "fix" mole fraction? Probably not, as it will violate assumptions
            # when we combine multiple units (ports) as long as mole_frac_comp
            # is a state variable, which I would rather not change.
            def sum_component_eq_rule(b):
                return 1.0 == sum(b.mole_frac_comp[j] for j in component_list)

            self.sum_component_eq = Constraint(
                rule=sum_component_eq_rule,
                doc=(
                    "Enforces that the sum of mole fractions equals one when\n"
                    "state variables are not already elsewhere defined."
                ),
            )

    def _flow_mol_comp(self):
        params = self.config.parameters
        component_list = params.component_list

        def flow_mol_comp_rule(b, j):
            return self.flow_mol * self.mole_frac_comp[j]

        self.flow_mol_comp = Expression(
            component_list,
            rule=flow_mol_comp_rule,
            doc="Molar flow rate of a particular component",
        )

    def _mw(self):
        params = self.config.parameters
        component_list = params.component_list
        mole_frac = self.mole_frac_comp
        mw_comp = {j: params.get_component(j).mw for j in component_list}
        self.mw = Expression(
            expr=sum(mole_frac[j] * mw_comp[j] for j in component_list),
            doc="Average molecular weight of the gas",
        )

    def _flow_mass(self):
        # I would like flow_mass to be a variable, so I can fix it.
        # However, I leave it as an expression here for consistency with
        # generic property packages. A flow_mass variable will be added
        # by pipeline unit models so that its derivatives may be computed.
        self.flow_mass = Expression(
            expr=self.mw * self.flow_mol,
            doc="Mass flow rate of the gas",
        )

    def _cp_mol_comp(self):
        params = self.config.parameters
        component_list = params.component_list
        comp = {j: params.get_component(j) for j in component_list}

        def cp_mol_comp_rule(b, j):
            return (
                comp[j].cp_mol_coef_A
                + comp[j].cp_mol_coef_B * b.temperature
                + comp[j].cp_mol_coef_C * b.temperature**2
                + comp[j].cp_mol_coef_D * b.temperature**3
            )

        self.cp_mol_comp = Expression(
            component_list,
            rule=cp_mol_comp_rule,
            doc=(
                "Pure component constant-pressure molar heat capacity "
                "of each component"
            ),
        )

    def _cp_mol(self):
        component_list = self.config.parameters.component_list
        self.cp_mol = Expression(
            expr=sum(
                self.mole_frac_comp[j] * self.cp_mol_comp[j] for j in component_list
            ),
            doc="Constant-pressure molar heat capacity of the gas mixture",
        )

    def _cp_mass(self):
        self.cp_mass = Expression(
            expr=self.cp_mol / self.mw,
            doc="Constant-pressure specific heat capacity of the gas mixture",
        )

    def _cv_mol_comp(self):
        params = self.config.parameters
        component_list = params.component_list
        comp = {j: params.get_component(j) for j in component_list}

        def cv_mol_comp_rule(b, j):
            return (
                comp[j].cv_mol_coef_A
                + comp[j].cv_mol_coef_B * b.temperature
                + comp[j].cv_mol_coef_C * b.temperature**2
                + comp[j].cv_mol_coef_D * b.temperature**3
            )

        self.cv_mol_comp = Expression(
            component_list,
            rule=cv_mol_comp_rule,
            doc=(
                "Pure component constant-volume molar heat capacity "
                "of each component"
            ),
        )

    def _cv_mol(self):
        component_list = self.config.parameters.component_list
        self.cv_mol = Expression(
            expr=sum(
                self.mole_frac_comp[j] * self.cv_mol_comp[j] for j in component_list
            ),
            doc="Constant-volume molar heat capacity of the gas mixture",
        )

    def _cv_mass(self):
        self.cv_mass = Expression(
            expr=self.cv_mol / self.mw,
            doc="Constant-volume specific heat capacity of the gas mixture",
        )

    def _heat_capacity_ratio(self):
        self.heat_capacity_ratio = Expression(
            expr=self.cp_mass / self.cv_mass,
            doc=(
                "Ratio of constant-volume to constant-pressure heat "
                "capacities of the gas mixture"
            ),
        )

    def _heat_capacity_ratio_phase(self):
        # Pipeline unit models require a heat_capacity_ratio_phase attribute
        # for compatibility with generic property packages.
        self.heat_capacity_ratio_phase = Expression(
            self.config.parameters.phase_list,
            expr=self.heat_capacity_ratio.expr,
            doc=(
                "Ratio of constant-volume to constant-pressure heat "
                "capacities of the gas mixture"
            ),
        )

    def _compressibility(self):
        # Compressibility is a param because here it is constant.
        # It could be a variable/expression, however, in a more complicated
        # model.
        self.compressibility = Param(
            initialize=0.80,
            doc="Compressibility factor of the gas",
        )

    def _dens_mol(self):
        # Make molar density a variable as it is important enough that
        # somebody might want to fix it.
        self.dens_mol = Var(
            initialize=1.0,
            units=pyunits.kmol / pyunits.m**3,
            doc="Molar density of the gas phase",
        )
        pressure = self.pressure
        gas_const = Constants.gas_constant
        compressibility = self.compressibility
        temperature = self.temperature
        dens_mol_expr = pyunits.convert(
            pressure / gas_const / compressibility / temperature,
            pyunits.kmol / pyunits.m**3,
        )
        self.dens_mol_eq = Constraint(
            expr=self.dens_mol == dens_mol_expr,
            doc=(
                "Equation used to calculate molar density -- "
                "ideal gas equation with a\ncompressibility factor"
            ),
        )

    def _dens_mol_comp(self):
        component_list = self.config.parameters.component_list

        def dens_mol_comp_rule(b, j):
            return self.mole_frac_comp[j] * self.dens_mol

        self.dens_mol_comp = Expression(
            component_list,
            rule=dens_mol_comp_rule,
            doc="Molar density of a particular component in the gas",
        )

    def _speed_of_sound(self):
        # Use a constraint here to make balance equation expressions more
        # legible.
        self.speed_of_sound = Var(
            initialize=300.0,
            units=pyunits.m / pyunits.s,
            doc="Speed of sound in the gas",
        )
        gas_const = pyunits.convert(
            Constants.gas_constant,
            pyunits.kJ / pyunits.kmol / pyunits.K,
        )
        speed_of_sound_expr = pyunits.convert(
            sqrt(
                self.heat_capacity_ratio
                * self.compressibility
                * gas_const
                * self.temperature
                / self.mw
            ),
            pyunits.m / pyunits.s,
        )
        self.speed_of_sound_eq = Constraint(
            expr=self.speed_of_sound == speed_of_sound_expr,
            doc="Equation to calculate speed of sound",
        )

    def define_state_vars(self):
        return {
            "flow_mol": self.flow_mol,
            "pressure": self.pressure,
            "temperature": self.temperature,
            "mole_frac_comp": self.mole_frac_comp,
        }

    def get_material_flow_terms(self, p, j):
        return self.flow_mol * self.mole_frac_comp[j]

    def get_material_density_terms(self, p, j):
        # Converting to kmol/km^3 is a workaround -- really I would like to
        # change the units of area and material_holdup in the control volume.
        return pyunits.convert(
            self.dens_mol_comp[j],
            pyunits.kmol / pyunits.km**3,
        )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def get_enthalpy_flow_terms(self, p):
        return (
            (self.temperature - self.config.parameters.temperature_ref)
            * self.cp_mol
            * self.flow_mol
        )
