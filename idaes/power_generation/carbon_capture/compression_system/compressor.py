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
Compression Model.
Reference:

Modekurti et al., (2017). "Design, Dynamic Modeling, and
Control of a Multistage CO2 Compresor System." International Journal
of Greenhouse Gas Control. v62., page 31-45
Created: April 2020
__Author__ = "Quang Minh Le, John Eslick, Andrew Lee"
"""
from pyomo.environ import (SolverFactory,
                           NonNegativeReals, Var,
                           Param, value, log,
                           sqrt, exp,
                           units as pyunits)
from pyomo.common.config import ConfigValue, In
from idaes.core import declare_process_block_class

from idaes.generic_models.unit_models.pressure_changer import (
    PressureChangerData, ThermodynamicAssumption)

import idaes.logger as idaeslog
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.constants import Constants as const
from enum import Enum

_log = idaeslog.getLogger(__name__)


class ImpellerType(Enum):
    cover_impeller = 1
    open_impeller = 2
    custom = 3


def _build_cover_impeller_equations(blk):
    blk.mass_flow_coeff = Var(blk.flowsheet().config.time,
                              initialize=0.0735,
                              doc="Compressor Flow Coefficient",
                              bounds=(0.01, 0.15))

    @blk.Constraint(blk.flowsheet().config.time)
    def impeller_work_coeff_eqn(b, t):
        phi = blk.mass_flow_coeff[t]
        return b.impeller_work_coeff[t] * 100 * phi == (
                          0.62 * phi - phi * (phi / 0.4)**3 + 0.0014) * 100

    @blk.Constraint(blk.flowsheet().config.time)
    def polytropic_head_coeff_vaned_diffuser_eqn(b, t):
        phi = blk.mass_flow_coeff[t]
        return b.mu_p_v[t] * 100 * phi == (
                      0.51 * phi + phi**2 - 7.6 * phi**3 - 0.00025) * 100


def _build_open_impeller_equations(blk):
    blk.mass_flow_coeff = Var(blk.flowsheet().config.time,
                              initialize=0.0735,
                              doc="Compressor Flow Coefficient",
                              bounds=(0.01, 0.15))

    @blk.Constraint(blk.flowsheet().config.time)
    def impeller_work_coeff_eqn(b, t):
        phi = blk.mass_flow_coeff[t]
        return b.impeller_work_coeff[t] * phi == (
                      0.68 * phi - phi * (phi / 0.37)**3 + 0.002)

    @blk.Constraint(blk.flowsheet().config.time)
    def polytropic_head_coeff_vaned_diffuser_eqn(b, t):
        phi = blk.mass_flow_coeff[t]
        return b.mu_p_v[t] * phi == 0.5 * (1 + (phi - 0.065) / abs(
            phi - 0.065)) * (
                0.59 * phi + 0.7 * phi**2 - 7.5 * phi**3 - 0.00025
                ) + 0.5 * (1 - (phi - 0.065) / abs(phi - 0.065)) * 0.6 * phi


class VaneDiffuserType(Enum):
    vane_diffuser = 1
    vaneless_diffuser = 2
    custom = 3


def _build_vane_diffuser_equations(blk):
    @blk.Constraint(blk.flowsheet().config.time)
    def polytropic_head_coeff_eqn(b, t):
        return b.mu_p[t] == b.mu_p_v[t]

    @blk.Constraint(blk.flowsheet().config.time)
    def polytropic_efficiency_eqn(b, t):
        return b.eff_p[t] == b.eff_p_v[t]


def _build_vaneless_diffuser_equations(blk):
    blk.mass_flow_coeff = Var(blk.flowsheet().config.time,
                              initialize=0.0735,
                              doc="Compressor Flow Coefficient",
                              bounds=(0.01, 0.15))

    @blk.Constraint(blk.flowsheet().config.time)
    def polytropic_head_coeff_eqn(b, t):
        return b.mu_p[t] == b.impeller_work_coeff[t] * b.eff_p[t]

    @blk.Constraint(blk.flowsheet().config.time)
    def polytropic_efficiency_eqn(b, t):
        phi = blk.mass_flow_coeff[t]
        return (b.eff_p[t] - b.eff_p_v[t]) * (
                      0.04 + 5 * phi + b.eff_p_v[t]**3) == -0.017


@declare_process_block_class("CompressionStage", doc="Compression Stage Model")
class CompressionStageData(PressureChangerData):
    # Pressure changer with isentropic compressor options
    CONFIG = PressureChangerData.CONFIG()
    CONFIG.compressor = True
    CONFIG.get("compressor")._default = True
    CONFIG.get("compressor")._domain = In([True])
    CONFIG.thermodynamic_assumption = ThermodynamicAssumption.isentropic
    CONFIG.get("thermodynamic_assumption")._default = \
        ThermodynamicAssumption.isentropic
    CONFIG.get('thermodynamic_assumption')._domain = \
        In([ThermodynamicAssumption.isentropic])
    CONFIG.declare(
        "impeller_type",
        ConfigValue(
            default=ImpellerType.open_impeller,
            domain=In(ImpellerType),
            description="Impeller type, if custom provide an expression rule",
            doc="""The type of impeller, if custom provide an expression rule
with the impeller_rule argument.
**default** - ImpellerType.open_impeller
**Valid values** - {
ImpellerType.cover_impeller,
ImpellerType.open_impeller,
ImpellerType.custom}""",
        ),
    )
    CONFIG.declare(
        "vane_diffuser_type",
        ConfigValue(
            default=VaneDiffuserType.vane_diffuser,
            domain=In(VaneDiffuserType),
            description="Vane diffuser type,"
            "if custom provide an expression rule",
            doc="""The type of vane diffuser,"
            "if custom provide an expression rule
with the vane_diffuser_rule argument.
**default** - VaneDiffuserType.vane_diffuser
**Valid values** - {
VaneDiffuserType.vane_diffuser,
VaneDiffuserType.vaneless_diffuser,
VaneDiffuserType.custom}""",
        ),
    )
    CONFIG.declare(
      "first_stage",
      ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Flag to define first stage",
        doc="If this is the first stage, user must provide IGV value."
        "Otherwise, user may provide the value of parameter A, B, C"))

    ################################################

    def build(self):
        super().build()

        #####################################################
        # First stage input
        if self.config.first_stage is None:
            raise ConfigurationError('User must provide a value for IGV'
                                     'if first stage, or values of A, B, C'
                                     'if other stages')
        # Declare variables for the model

        # Some shorter refernces to property blocks
        properties_in = self.control_volume.properties_in
        properties_out = self.control_volume.properties_out

        self.mass_flow_coeff = Var(self.flowsheet().config.time,
                                   initialize=0.0735,
                                   doc="Compressor Flow Coefficient",
                                   bounds=(0.01, 0.15),
                                   units=pyunits.dimensionless)
        self.r2 = Var(initialize=0.075,
                      doc="Impeller Tip Radius",
                      units=pyunits.m)
        self.rspeed = Var(self.flowsheet().config.time,
                          initialize=1500,
                          doc='Rotation Speed of The Impeller',
                          bounds=(0.5, 8000),
                          units=1 / pyunits.s)
        self.U2 = Var(self.flowsheet().config.time,
                      initialize=315,
                      doc="Impeller Tip Speed",
                      bounds=(0, 4000),
                      units=pyunits.m / pyunits.s)
        self.delta_enth_polytropic = Var(self.flowsheet().config.time,
                                         initialize=2548.5,
                                         doc="Polytropic Enthalpy Change",
                                         units=pyunits.J / pyunits.mol,
                                         within=NonNegativeReals)
        self.impeller_work_coeff = Var(self.flowsheet().config.time,
                                       initialize=0.6993,
                                       doc="Impeller Work Coefficient",
                                       bounds=(0.001, 1))
        self.mu_p_v = Var(self.flowsheet().config.time,
                          initialize=0.60,
                          doc="Polytropic Head Coefficient Vaned Diffuser",
                          bounds=(0.001, 1))
        self.mu_p = Var(self.flowsheet().config.time,
                        initialize=0.60,
                        doc="Polytropic Head Coefficient",
                        bounds=(0.001, 1))
        self.eff_p_v = Var(self.flowsheet().config.time,
                           initialize=0.85,
                           doc="Polytropic Efficiency Vaned Diffuser",
                           bounds=(0.001, 1))
        self.eff_p = Var(self.flowsheet().config.time,
                         initialize=0.85,
                         doc="Polytropic Efficiency",
                         bounds=(0.001, 1.0))
        self.Ma = Var(self.flowsheet().config.time,
                      initialize=0.1,
                      doc='Rotational Mach Number',
                      bounds=(0, 5))
        self.Ma_max = Var(self.flowsheet().config.time,
                          initialize=0.1,
                          doc='Rotational Mach Number Upper limit',
                          bounds=(0, 5))
        # dynamic related variables
        self.psi_3 = Var(self.flowsheet().config.time,
                         initialize=0.01,
                         doc="Dimensionless Exit Flow Coefficient",
                         within=NonNegativeReals,
                         units=pyunits.dimensionless)
        self.psi_s = Var(self.flowsheet().config.time,
                         initialize=0.01,
                         doc="Dimensionless Isentropic Head Coefficient",
                         within=NonNegativeReals,
                         units=pyunits.dimensionless)
        ########################################################

        # Declare variables for the model
        self.z_s = Var(initialize=0.97373,
                       doc="Compressibility factors at suction")
        self.z_d1 = Var(initialize=0.88949,
                        doc="variable used for calculating"
                        "compressibility factor at discharge pressure")
        self.efficiency_mech = Var(initialize=0.97,
                                   doc="Mechanical Efficiency")
        self.eff_drive = Var(initialize=1.0,
                             doc="Driver efficiency")

        ############################################################

        self.ratioP[:] = 2.0   # make sure these have a number value
        self.deltaP[:] = 0     # to avoid an error later in initialize

        ###########################################################
        # Tip speed, mass flow coefficient, Mach number, and pressure ratio
        # suggested by John and Andrew, I assume that the property package
        # calculates the speed of sound
        @self.Constraint(self.flowsheet().config.time, doc="Mach Number")
        def Ma_con(b, t):
            # not sure if there are two phases at the inlet
            # but I assume that this contains single vapor phase.
            speed_of_sound = self.control_volume.properties_in[
                t].speed_sound_phase['Vap']
            return b.Ma[t] == b.U2[t] / speed_of_sound

        @self.Constraint(self.flowsheet().config.time, doc='Pressure Ratio')
        def Pratio_con(b, t):
            return properties_out[t].pressure == b.ratioP[t] * \
                  properties_in[t].pressure

        @self.Constraint(self.flowsheet().config.time, doc="Rotation Speed "
                         "of the Impeller")
        def rspeed_con(b, t):
            return b.U2[t] == 2 * const.pi * b.r2 * b.rspeed[t]

        # set up the vane diffuser rule.
        vdselect = self.config.vane_diffuser_type
        if vdselect is not VaneDiffuserType.custom:
            _log.warning("A vane diffuser callback was provided"
                         "but the valve diffuser type is not custom.")
        if vdselect == VaneDiffuserType.vane_diffuser:
            _build_vane_diffuser_equations(self)
        elif vdselect == VaneDiffuserType.vaneless_diffuser:
            _build_vaneless_diffuser_equations(self)

        # set up the impeller rule.
        iselect = self.config.impeller_type
        if iselect is not ImpellerType.custom:
            _log.warning("An impeller callback was provided but the impeller"
                         "type is not custom.")
        if iselect == ImpellerType.cover_impeller:
            _build_cover_impeller_equations(self)
        elif iselect == ImpellerType.open_impeller:
            _build_open_impeller_equations(self)

        @self.Constraint(self.flowsheet().config.time, doc="Polytropic "
                         "Efficiency")
        def eff_p_v_cons(b, t):
            return b.eff_p_v[t] * b.impeller_work_coeff[t] == b.mu_p_v[t]

        # --------------------------------------------------------------------

        # Calculate total enthalpy and entropy change through the stage
        @self.Expression(self.flowsheet().config.time, doc="Specific "
                         "Enthalpy Change of Isentropic Process")
        def delta_enth_isentropic(b, t):
            return b.properties_isentropic[t].enth_mol - \
                properties_in[t].enth_mol

        @self.Expression(self.flowsheet().config.time, doc="Actual Enthalpy "
                         "Change")
        def delta_enth_actual(b, t):
            return properties_out[t].enth_mol - properties_in[t].enth_mol

        @self.Expression(self.flowsheet().config.time, doc="Entropy change in"
                         "Compressor Stage")
        def deltaS(b, t):
            return properties_out[t].entr_mol - properties_in[t].entr_mol

        @self.Constraint(self.flowsheet().config.time, doc="Isentropic"
                         "efficiency")
        def eff_isen_eqn(b, t):
            eff_isen = b.efficiency_isentropic[t]
            return eff_isen * b.delta_enth_actual[t] == \
                b.delta_enth_isentropic[t]

        @self.Constraint(self.flowsheet().config.time, doc="Polytropic "
                         "Enthalpy Correlation")
        def polytropic_correlation(b, t):
            return b.delta_enth_polytropic[t] == b.mu_p[t] * \
                b.delta_enth_actual[t] / b.impeller_work_coeff[t]

        # This equation of polytropic enthalpy is derived from Eq. (2.13)
        # From Aungier textbook (2000)
        @self.Constraint(self.flowsheet().config.time, doc="Equation for"
                         "Polytropic Enthalpy")
        def delta_enth_polytropic_con(b, t):
            Tout = properties_out[t].temperature
            Tin = properties_in[t].temperature
            return b.delta_enth_polytropic[t] == b.delta_enth_actual[t] \
                - b.deltaS[t] * (Tout - Tin)/log(Tout / Tin)

        # --------------------------------------------------------------------
        # Design variable constraints
        @self.Expression(self.flowsheet().config.time, doc="Maximum Allowable"
                         "Impeller Tip Speed")
        def U2max(b, t):
            phi = b.mass_flow_coeff[t]
            return sqrt((1984.1 * phi**2 - 616.88 * phi + 215.97
                         ) * 0.7 * 830)

        @self.Expression(self.flowsheet().config.time,
                         doc="Thermodynamic Power")
        def fluid_pow(b, t):
            flow = properties_in[t].flow_mol
            return b.delta_enth_actual[t] * flow

        @self.Expression(self.flowsheet().config.time, doc="Shaft Power")
        def elec_pow(b, t):
            return b.fluid_pow[t] / (b.efficiency_mech * b.eff_drive)

        # -------------------------------------------------------------------
        @self.Expression(self.flowsheet().config.time,
                         doc="Static Volumetric Flow at Suction Inlet")
        def Vsdot(b, t):
            flow_in = properties_in[t].flow_mol
            rho_v0 = properties_in[t].dens_mol
            return flow_in / rho_v0

        @self.Constraint(self.flowsheet().config.time, doc="Mass Flow "
                         "Coefficient")
        def mass_flow_coeff_eqn(b, t):
            phi = b.mass_flow_coeff[t]
            return phi == b.Vsdot[t] / (const.pi * b.r2**2 * b.U2[t])

        @self.Expression(self.flowsheet().config.time,
                         doc="Static volumetric flow at impeller exit")
        def V3dot(b, t):
            Tin = properties_in[t].temperature
            Tout = properties_out[t].temperature
            Pratio = b.ratioP[t]
            b.z_d = 0.5 * (b.z_s + b.z_d1)
            return b.Vsdot[t] * (b.z_d / b.z_s) * (Tout / Tin) * (1 / Pratio)

        @self.Constraint(self.flowsheet().config.time,
                         doc="Dimensionless Exit Flow Coefficient")
        def psi_3_eqn(b, t):
            b2 = 0.076 * 2 * b.r2
            return b.psi_3[t] == b.V3dot[t] / (
                const.pi * 2 * b.r2 * b2 * b.U2[t])

        @self.Expression(self.flowsheet().config.time,
                         doc="Isentropic Head")
        def ys_model(b, t):
            k_v = self.control_volume.properties_in[t].heat_capacity_ratio
            a = (k_v - 1) / k_v
            Pratio = b.ratioP[t]
            Tin = properties_in[t].temperature
            mw = properties_in[t].mw
            gas_const = const.gas_constant
            return b.z_s * (gas_const / mw) * Tin * (1 / a) * (Pratio**a - 1)

        @self.Constraint(self.flowsheet().config.time,
                         doc="Dimensionless Isentropic Head Coefficient")
        def psi_s_eqn(b, t):
            return b.psi_s[t] == 2 * b.ys_model[t] / (b.U2[t]**2)

        # Ang is used to estimate coeffcients a, b, c for first stage
        if self.config.first_stage is True:
            self.Ang = Var(self.flowsheet().config.time,
                           initialize=0,
                           doc="Inlet Guide Vanes Angle",
                           bounds=(-15, 90), units=pyunits.dimensionless)

            @self.Expression(self.flowsheet().config.time)
            def coeff_a(b, t):
                return -38.46 * exp(-0.009058 * b.Ang[t]) + (
                    - 0.1043) * exp(0.09821 * b.Ang[t])

            @self.Expression(self.flowsheet().config.time)
            def coeff_b(b, t):
                return 6.863 * exp(-0.01329 * b.Ang[t]) \
                    + 0.02018 * exp(0.09443 * b.Ang[t])

            @self.Expression(self.flowsheet().config.time)
            def coeff_c(b, t):
                return -0.00891 * exp(0.06757 * b.Ang[t]) \
                    + 0.9661 * exp(0.006228 * b.Ang[t])
        else:
            self.coeff_a = Param(default=-428.64,
                                 mutable=True,
                                 doc="Performance Curve Parameter A")
            self.coeff_b = Param(default=29.178,
                                 mutable=True,
                                 doc="Perforamance Curve Parameter B")
            self.coeff_c = Param(default=0.9737,
                                 mutable=True,
                                 doc="Performance Curve parameter C")

        @self.Constraint(self.flowsheet().config.time,
                         doc="Correlation between psi_s and psi_3")
        def psi_s_stage_eqn(b, t):
            if self.config.first_stage is True:
                return b.psi_s[t] == b.coeff_a[t] * b.psi_3[t]**2 + \
                    b.coeff_b[t] * b.psi_3[t] + b.coeff_c[t]
            else:
                return b.psi_s[t] == b.coeff_a * b.psi_3[t]**2 + \
                    b.coeff_b * b.psi_3[t] + b.coeff_c

    def initialize(self, state_args={}, outlvl=idaeslog.NOTSET, solver='ipopt',
                   optarg={'tol': 1e-6, 'max_iter': 50}):
        """
        Initialize the inlet compressor stage model.
        This deactivates the specialized constraints,
        then does the isentropic compressor initialization,
        then reactivates the constraints and solves.

        Args:
            state_args (dict): Initial state for property initialization
            outlvl (int): Amount of output (0 to 3) 0 is lowest
            solver (str): Solver to use for initialization
            optarg (dict): Solver arguments dictionary
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        opt = SolverFactory(solver)
        opt.options = optarg

        in_flags = self.control_volume.properties_in.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
            state_args=state_args,
        )
        out_flags = self.control_volume.properties_out.initialize(
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=False,
            state_args=state_args,
        )

        for t in self.flowsheet().config.time:
            # if there is not a good guess for efficiency or outlet pressure
            # provide something reasonable.
            eff = self.efficiency_isentropic[t]
            eff.fix(eff.value if value(eff) > 0.3 and value(eff
                                                            ) < 1.0 else 0.85)
            # check for alternate pressure specs
            #
            if self.outlet.pressure[t].fixed:
                self.ratioP[t] = value(
                    self.outlet.pressure[t] / self.inlet.pressure[t])
            elif self.ratioP[t].fixed:
                self.outlet.pressure[t] = value(
                    self.ratioP[t] * self.inlet.pressure[t])

            # intialize isentropic property block
            isen_state = {}
            t_init = self.flowsheet().config.time.first()
            for i in self.properties_isentropic[
                    t_init].define_state_vars().keys():
                state_var = getattr(self.control_volume.properties_out[
                    t_init], i)
                if state_var.is_indexed:
                    isen_state[i] = {}
                    for j in state_var:
                        isen_state[i][j] = value(state_var[j])
                else:
                    isen_state[i] = value(state_var)
            self.properties_isentropic.initialize(
                state_args=isen_state, solver=solver, optarg=optarg)

        # Deactivate special constraints
        self.eff_isen_eqn.deactivate()
        self.Pratio_con.deactivate()
        self.Ma_con.deactivate()
        self.rspeed_con.deactivate()
        self.eff_p_v_cons.deactivate()
        self.polytropic_correlation.deactivate()
        self.delta_enth_polytropic_con.deactivate()
        self.mass_flow_coeff_eqn.deactivate()
        self.psi_3_eqn.deactivate()
        self.psi_s_eqn.deactivate()
        self.psi_s_stage_eqn.deactivate()

        # call pressure changer initialzation
        super().initialize(
            state_args=state_args, outlvl=outlvl, solver=solver, optarg=optarg
        )
        # Fix variables
        self.Ma.fix()
        self.rspeed.fix()
        self.impeller_work_coeff.fix()
        self.delta_enth_polytropic.fix()
        self.mass_flow_coeff.fix()
        self.psi_3.fix()
        self.psi_s.fix()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)

        init_log.info_high(
                "Initialization Step 1 {}.".format(idaeslog.condition(res))
            )

        # Activate special constraints
        self.eff_isen_eqn.activate()
        self.Pratio_con.activate()
        self.Ma_con.activate()
        self.rspeed_con.activate()
        self.eff_p_v_cons.activate()
        self.polytropic_correlation.activate()
        self.delta_enth_polytropic_con.activate()
        self.mass_flow_coeff_eqn.activate()
        self.psi_3_eqn.activate()
        self.psi_s_eqn.activate()
        self.psi_s_stage_eqn.activate()

        # Unfix variables
        self.efficiency_isentropic.unfix()
        self.outlet.pressure.unfix()
        self.Ma.unfix()
        self.rspeed.unfix()
        self.impeller_work_coeff.unfix()
        self.delta_enth_polytropic.unfix()
        self.mass_flow_coeff.unfix()
        self.psi_3.unfix()
        self.psi_s.unfix()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high(
                "Initialization Step 2 {}.".format(idaeslog.condition(res))
            )

        self.control_volume.properties_in.release_state(
            flags=in_flags, outlvl=outlvl)
        self.control_volume.properties_out.release_state(
            flags=out_flags, outlvl=outlvl)
        # self.control_volume.release_state(flags=flags, outlvl=outlvl)
        init_log.info(
            "Initialization Complete: {}".format(idaeslog.condition(res))
        )
