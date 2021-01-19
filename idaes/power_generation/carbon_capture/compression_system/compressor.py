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
Centrifugal compressor stage model.  This model is based on:

Modekurti et al., (2017). "Design, Dynamic Modeling, and Control of a Multistage CO2
				Compresor System." International Journal of Greenhouse Gas Control. v62.
Created: April 2020
__Author__ = "Quang Minh Le "
"""

from __future__ import division
# import Pyomo
from pyomo.environ import (ConcreteModel, SolverFactory, Reals, NonNegativeReals,
    Var, Param, Expression, Constraint, SolverFactory, value, exp, log, sqrt)
from pyomo.common.config import ConfigBlock, ConfigValue, In


# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        FlowsheetBlock)
from idaes.generic_models.unit_models.pressure_changer import (PressureChangerData,
																ThermodynamicAssumption)

from idaes.generic_models.properties.helmholtz.helmholtz import (
    HelmholtzThermoExpressions as ThermoExpr
)
import idaes.logger as idaeslog
from idaes.core.util import from_json, to_json, StoreSpec

from idaes.core.util.model_statistics import degrees_of_freedom

import idaes
import math


_log = idaeslog.getLogger(__name__)

@declare_process_block_class("InletStage", doc="Inlet stage of compressor model")
class InletStageData(PressureChangerData):
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

    ################################################

    def build(self):
        super(InletStageData, self).build()
        config = self.config # shorter config pointer

        # The thermodynamic expression writer object, te, writes expressions
        # including external function calls to calculate thermodynamic quantities
        # from a set of state variables.
        te = ThermoExpr(blk=self, parameters=config.property_package)

        if SolverFactory('ipopt').available():
          solver = SolverFactory('ipopt')
          solver.options = {'tol':1e-6,
                            'max_iter': 100,
                            'linear_solver': 'ma27'}
        else:
          solver = None

        #####################################################
        # Declare variables for the model

        # Some shorter refernces to property blocks
        properties_in = self.control_volume.properties_in
        properties_out = self.control_volume.properties_out

        
        self.kappaT = Var(self.flowsheet().config.time,
                            initialize=5.0e-7)
        self.c0 = Var(self.flowsheet().config.time,
                        initialize = 260,
                        doc="Speed of sound in the inlet gas stream",
                        bounds=(0,4000)) # m/s
        self.mass_flow_coeff = Var(self.flowsheet().config.time,
                            initialize=0.0735,
                            doc="Compressor flow coefficient",
                            bounds=(0.01,0.15))        
        self.I = Var(self.flowsheet().config.time,
                            initialize=0.6993,
                            doc="Impeller work coefficient",
                            bounds=(0.001,1))
        self.mu_p_v = Var(self.flowsheet().config.time,
                            initialize = 0.60,
                            doc="Polytropic head coefficient vaned diffuser",
                            bounds=(0.001,1))        
        self.mu_p = Var(self.flowsheet().config.time,
                            initialize = 0.60,
                            doc="Polytropic head coefficient",
                            bounds=(0.001,1))
        self.eff_p_v = Var(self.flowsheet().config.time,
                            initialize=0.85,
                            doc="Polytropic stage efficiency vaned diffuser",
                            bounds=(0.001,1))
        self.eff_p = Var(self.flowsheet().config.time,
                            initialize=0.85,
                            doc="Polytropic stage efficiency",
                            bounds = (0.001, 1.0))       
        
        self.delta_enth_polytropic = Var(self.flowsheet().config.time,
                            initialize=2548.5,
                            doc="Polytropic enthalpy change [J/mol]",
                            within=NonNegativeReals)
        
        self.U2 = Var(self.flowsheet().config.time,
                            initialize = 315, 
                            doc="Impeller tip speed", 
                            bounds= (0,4000)) # upper: 4000, lower: 0 // m/s
        self.Ma = Var(self.flowsheet().config.time,
                            initialize = 0.1,
                            doc = 'Rotational Ma number',
                            bounds = (0,5))
        self.rspeed = Var(self.flowsheet().config.time,
                            initialize = 1500,
                            doc = 'Rotation speed of the impeller',
                            bounds = (0.5, 8000))      
        self.r2 = Var(initialize= 0.075, 
                            doc="Impeller tip radius")  
        
        ########################################################     
        
        self.dp = Param(default= 1, within=Reals, doc='a pressure difference for finite difference derivative calculation')
        self.diffuser_type = Param(default=1.0, 
                                    mutable=True,
                                    doc="1 for vane diffuser, 2 for vaneless diffuser")
        self.impeller_type = Param(default=2.0,
                                    mutable=True,
                                    doc="1 for covered impeller, 2 for open impeller")
        self.efficiency_mech = Param(default= 0.97,
                                    doc="Compressor mechanical efficiency")
        self.eff_drive = Param(default= 1.0,
                                doc="Compressor driver efficiency")        

        ############################################################ 
        # Fix variables 
        self.ratioP[:] = 2.0   # make sure these have a number value  
        self.deltaP[:] = 0     #   to avoid an error later in initialize
    
        ###########################################################
        # Speed of sound stuff
        # finite difference derivative for kappa_T calculation
        @self.Constraint(self.flowsheet().config.time, doc="Isothermal compressibility [m2/N]")
        def kappaT_con(b,t):
            Ps = properties_in[t].pressure            
            P0_f = Ps + b.dp
            rho_v0 = properties_in[t].dens_mol
            # try thermo expression to estimate density fraction, assume it be a vapor phase
            # rho_v0f = te.rho_mol_vap(T= value(properties_in[t].pressure), p= value(P0_f), x=1)
            rho_v0f = (P0_f/Ps)*rho_v0                             # TODo: how to call this incremental density
            return b.kappaT[t]*b.dp == 1 - rho_v0/rho_v0f
        
        @self.Constraint(self.flowsheet().config.time, doc="Speed of sound")
        def c0_con(b,t):
            Cp0 = properties_in[t].cp_mol
            Cv0 = properties_in[t].cv_mol
            mw =  properties_in[t].mw
            rho_v0 = properties_in[t].dens_mol
            return b.c0[t]**2*(Cv0 * mw * b.kappaT[t] * rho_v0) == Cp0

        # ---------------------------------------------------------------------------------
        # Tip speed, mass flow coefficient, Mach number, and pressure ratio
        @self.Constraint(self.flowsheet().config.time, doc="Mach number")
        def Ma_con(b,t):
            return b.Ma[t] == b.U2[t]/b.c0[t]

        @self.Constraint(self.flowsheet().config.time, doc='Pressure ratio')
        def Pratio_con(b,t):
            return properties_out[t].pressure == b.ratioP[t] * properties_in[t].pressure

        @self.Constraint(self.flowsheet().config.time, doc="Rotation speed of the impeller")
        def rspeed_con(b,t):
            return  b.U2[t] == b.rspeed[t] * b.r2 

        @self.Constraint(self.flowsheet().config.time, doc="Mass flow coefficient")
        def mass_flow_coeff_eqn(b,t):
            flow = properties_in[t].flow_mol
            rho_v0 = properties_in[t].dens_mol
            phi = b.mass_flow_coeff[t]
            return phi*(3.14159*rho_v0*b.r2**2*b.U2[t]) == flow

        # -----------------------------------------------------------------------------------        
        
        @self.Constraint(self.flowsheet().config.time)
        def specification_impeller_type_rule_1(b,t):
            phi = b.mass_flow_coeff[t]
            if value(b.impeller_type) < 1.5:
                return b.I[t]*100*phi == (0.62*phi - phi*(phi/0.4)**3 + 0.0014)*100
            else:
                return 100*b.I[t]*phi == (0.68*phi - phi*(phi/0.37)**3+0.002)*100

        @self.Constraint(self.flowsheet().config.time)
        def specification_impeller_type_rule_2(b,t):
            phi = b.mass_flow_coeff[t]
            if value(b.impeller_type) < 1.5:
                return b.mu_p_v[t]*100*phi == (0.51*phi + phi**2 -7.6*phi**3-0.00025)*100
            else:
                # return b.mu_p_v[t]*phi == 0.5*(1+(phi-0.065)/abs(phi-0.065))*(0.59*phi+0.7*phi**2-7.5*phi**3-0.00025)+\
                #                     0.5*(1-(phi-0.065)/abs(phi-0.065))*0.6*phi
                return 100*b.mu_p_v[t]*phi == (0.59*phi+0.7*phi**2-7.5*phi**3-0.00025)*100

        @self.Constraint(self.flowsheet().config.time)
        def specification_diffuser_type_rule_1(b,t):
            phi = b.mass_flow_coeff[t]
            if value(b.diffuser_type) < 1.5:
                return b.mu_p[t] == b.mu_p_v[t]
            else:
                return b.mu_p[t] == b.I[t]*b.eff_p[t]

        @self.Constraint(self.flowsheet().config.time)
        def specification_diffuser_type_rule_2(b,t):
            if value(b.diffuser_type) < 1.5:
                return b.eff_p[t] == b.eff_p_v[t]
            else:
                return (b.eff_p[t]-b.eff_p_v[t])*(0.04+5*phi+b.eff_p_v[t]**3) == -0.017

        @self.Constraint(self.flowsheet().config.time, doc="Polytropic efficiency")
        def eff_p_v_cons(b,t):
            return b.eff_p_v[t]*b.I[t] == b.mu_p_v[t]

        #------------------------------------------------------------------------------------
        # Need to extend for mixture
        # calculate total enthalpy and entropy change through the stage
         

        @self.Expression(self.flowsheet().config.time, doc="Specific enthalpy change of isentropic process")
        def delta_enth_isentropic(b,t):
            return b.properties_isentropic[t].enth_mol - properties_in[t].enth_mol  

        @self.Expression(self.flowsheet().config.time, doc="Actual enthalpy change")
        def delta_enth_actual(b,t):
            return properties_out[t].enth_mol - properties_in[t].enth_mol 

        @self.Expression(self.flowsheet().config.time, doc="Entropy change in compressor stage")
        def deltaS(b,t):
            return properties_out[t].entr_mol - properties_in[t].entr_mol   

        @self.Constraint(self.flowsheet().config.time, doc="Isentropic efficiency")
        def eff_isen_eqn(b,t):   
            eff = b.efficiency_isentropic[t]     
            return b.delta_enth_isentropic[t] == eff * b.delta_enth_actual[t]      

        @self.Constraint(self.flowsheet().config.time)
        def deltaH_t_con(b,t):
            mw = properties_in[t].mw
            return b.delta_enth_actual[t] == b.I[t] * b.U2[t]**2 * mw

        @self.Constraint(self.flowsheet().config.time)
        def deltaH_p_con(b,t):
            mw = properties_in[t].mw
            return b.delta_enth_polytropic[t] == b.mu_p[t] * self.U2[t]**2 * mw      

        # This equation of polytropic enthalpy is based on Eq. (2.13) From Aungier textbook (2000)
        @self.Expression(self.flowsheet().config.time, doc="Equation for polytropic enthalpy") 
        def delta_enth_polytropic_con(b,t):
            Tout = properties_out[t].temperature
            Tin = properties_in[t].temperature
            return b.delta_enth_polytropic[t] == b.delta_enth_actual[t]-b.deltaS[t]*(Tout-Tin)/log(Tout/Tin)     
        
        
        # -------------------------------------------------------------------------------------------------------
        # design variable constraints
        @self.Constraint(self.flowsheet().config.time, doc="Maximum allowable impeller tip speed") 
        def U2max_cons(b,t):
            phi = b.mass_flow_coeff[t]
            return b.U2[t] <= sqrt((1984.1*phi**2 - 616.88*phi + 215.97)*0.7*830)

        # @self.Constraint(self.flowsheet().config.time, doc='Maximum Ma number')
        # def Mamax_cons(b,t):
        #     # n_stage = 8
        #     return b.Ma[t] <= -0.202*log(5) + 1.25

        @self.Expression(self.flowsheet().config.time, doc="Thermodynamic power")
        def fluid_pow(b,t):
            flow = properties_in[t].flow_mol 
            return b.delta_enth_actual[t]*flow

        @self.Expression(self.flowsheet().config.time, doc="Shaft power")
        def elec_pow(b,t):
            return b.fluid_pow[t]/(b.efficiency_mech * b.eff_drive)

    def initialize(self, state_args={}, outlvl=6, solver='ipopt',
        optarg={'tol': 1e-6, 'max_iter':1000}):
        """
        Initialize the inlet compressor stage model.  This deactivates the
        specialized constraints, then does the isentropic compressor initialization,
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

        flags = self.control_volume.initialize(outlvl=outlvl-1,
                                       optarg=optarg,
                                       solver=solver,
                                       state_args=state_args)

        print('dof for step 1 =', degrees_of_freedom(self))
        init_log.info_high("Initialization Step 1 Complete.")
        
        # sp is what to save to make sure state after init is same as the start
        #   saves value, fixed, and active state, doesn't load originally free
        #   values, this makes sure original problem spec is same but initializes
        #   the values of free vars
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)           
        
        for t in self.flowsheet().config.time:
            for k, v in self.inlet.vars.items():
                v[t].fix()
            for k, v in self.outlet.vars.items():
                v[t].unfix()

            # if there is not a good guess for efficiency or outlet pressure
            # provide something reasonable.
            eff = self.efficiency_isentropic[t]
            eff.fix(eff.value if value(eff) > 0.3 and value(eff) < 1.0 else 0.85)
            # for outlet pressure try outlet pressure, pressure ratio
            # then if none of those look reasonable use a pressure ratio of 2.0
            # to calculate outlet pressure
            Pout = self.outlet.pressure[t]
            Pin = self.inlet.pressure[t]
            prdp = value((self.deltaP[t] - Pin)/Pin)
            if value(Pout / Pin) < 4.0 or value(Pout / Pin) > 1.0:
                if value(self.ratioP[t]) < 4.0 and value(self.ratioP[t]) > 1.0:
                    Pout.fix(value(Pin * self.ratioP[t]))
                elif prdp < 4.0 and prdp > 1.0:
                    Pout.fix(value(prdp*Pin))
                else:
                    Pout.fix(value(Pin * 2.0))
            else:
                Pout.fix()

        self.deltaP[:] = value(Pout - Pin)
        self.ratioP[:] = value(Pout / Pin)
        
        for t in self.flowsheet().config.time:
            self.properties_isentropic[t].pressure.value = value(
                self.outlet.pressure[t]
            )
            self.properties_isentropic[t].flow_mol.value = value(self.inlet.flow_mol[t])
            self.properties_isentropic[t].enth_mol.value = value(
                self.inlet.enth_mol[t] * 1.2
            )
            self.outlet.flow_mol[t].value = value(self.inlet.flow_mol[t])
            self.outlet.enth_mol[t].value = value(self.inlet.enth_mol[t] * 1.2)

        # fix variable
        self.mass_flow_coeff.fix()

        # Deactivate special constraints  
        self.mass_flow_coeff_eqn.deactivate()  
        self.eff_isen_eqn.deactivate()
        self.Pratio_con.deactivate()        
        
        print('dof for step 2 =', degrees_of_freedom(self))
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high(
                "Initialization Step 2 {}.".format(idaeslog.condition(res))
            ) 

        # free eff_isen and activate special constraints
        self.mass_flow_coeff_eqn.activate()
        self.eff_isen_eqn.activate()
        self.Pratio_con.activate()    

        self.efficiency_isentropic.unfix()
        self.outlet.pressure.unfix()
        # self.mass_flow_coeff.unfix()
       
        
        print('dof for step 3 =', degrees_of_freedom(self))
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(self, tee=slc.tee)
        init_log.info_high(
                "Initialization Step 3 {}.".format(idaeslog.condition(res))
            ) 

        init_log.info(
            "Initialization Complete: {}".format(idaeslog.condition(res))
        )
        # reload original spec
        from_json(self, sd=istate, wts=sp)
