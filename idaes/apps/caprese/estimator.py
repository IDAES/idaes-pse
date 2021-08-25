# -*- coding: utf-8 -*-
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
""" Block-like object meant for Moving Horizon Estimator (MHE) models.
"""

__author__ = "Kuan-Han Lin"

# import idaes.logger as idaeslog
from idaes.apps.caprese.categorize import (
        # categorize_dae_variables,
        _identify_derivative_if_differential,
        categorize_dae_variables_and_constraints,
        CATEGORY_TYPE_MAP,
        )
from idaes.apps.caprese.dynamic_var import (
        DynamicVar,
        _DynamicVector,
        DiffVar,
        AlgVar,
        InputVar,
        DerivVar,
        FixedVar,
        MeasuredVar,
        ActualMeasurementVar,
        MeasurementErrorVar,
        ModelDisturbanceVar,
        )
from idaes.apps.caprese.dynamic_block import (
        _DynamicBlockData,
        IndexedDynamicBlock,
        DynamicBlock,
        )
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.apps.caprese.common.config import (
        VariableCategory as VC,
        ConstraintCategory as CC,
        )
from pyomo.environ import (
        Var,
        Objective,
        Constraint,
        Block,
        Reference,
        Set,
        TerminationCondition,
        )
from pyomo.core.base.block import _BlockData
from pyomo.common.collections import ComponentMap, ComponentSet
from pyomo.dae.flatten import (
        flatten_dae_components,
        )
from pyomo.core.base.indexed_component import UnindexedComponent_set
from pyomo.dae.contset import ContinuousSet
from pyomo.dae.set_utils import deactivate_model_at


MHE_CATEGORY_TYPE_MAP = {
    VC.ACTUALMEASUREMENT: ActualMeasurementVar,
    VC.MEASUREMENTERROR: MeasurementErrorVar,
    VC.MODELDISTURBANCE: ModelDisturbanceVar,
                        }

class _EstimatorBlockData(_DynamicBlockData):
    """ This class adds methods for building MHE and for working with dynamic
    models. 
    """

    def _construct(self):
        """This function builds necessary variables and constraints for MHE.
            e.g. actual measurements, measurement errors, model disturbances
                error constraints, disturbed differential equations
        """
        self._categorize_constraints = True
        
        if self._sample_time is None:
            raise RuntimeError("MHE needs the sample time to be provided!")

        # Call base class construct method
        super(_EstimatorBlockData, self)._construct()
            
        # This processing should not be necessary here; we should have done it
        # in the base class constructor
        #
        ## MHE requires not only variable categories but also constraint categories.
        ## This if statement check which is given, which is not.
        #if self._category_dict is not None and self._con_category_dict is not None:
        #    unref_category_dict = self._use_user_given_var_categ_dict(inputs, measurements)
        #    unref_con_category_dict = self._con_category_dict

        ## This is the case that I expect to be the most often.
        #elif self._category_dict is None and self._con_category_dict is None:
        #    unref_category_dict, unref_con_category_dict = \
        #        self._categorize_var_con_for_MHE(
        #            mod, time, inputs, measurements
        #        )
        #else:
        #    raise RuntimeError(
        #        "Please give both cateogry_dict and con_category_dict for MHE. "
        #        "Note that the order of differential variables and "
        #        "differential equations should match."
        #    )

        #self.reference_var_category_dict(unref_category_dict)
        self.reference_var_category_dict(self.category_dict)
        #At this stage, no need to use reference on constraints
        #self.con_category_dict = unref_con_category_dict

        #Set a flag so that classification won't do again in _DynamicBlockData
        self.already_categorized_for_MHE = True
        category_dict = self.category_dict

        self._add_sample_point_set()
        # Variables and constraints for MHE will be added under this block.
        self.MHE_VARS_CONS_BLOCK = Block()

        n_measurement = len(category_dict[VC.MEASUREMENT])
        n_diffvar_con = len(category_dict[VC.DIFFERENTIAL])
        self.MHE_VARS_CONS_BLOCK.MEASUREMENT_SET = Set(initialize = range(n_measurement))
        self.MHE_VARS_CONS_BLOCK.DIFFERENTIAL_SET = Set(initialize = range(n_diffvar_con))

        self._add_actual_measurement_param()
        self._add_measurement_error()
        self._add_model_disturbance()

        #self._add_MHE_vars_to_category_dict()
        # Here I inline _add_MHE_vars_to_category_dict, for now, so
        # I can create a new category dict with MHE vars
        #
        # These are new categories we need to add to the "vectors" block
        mhe_var_category_dict = {}
        MHEBlock = self.MHE_VARS_CONS_BLOCK
        # target_set = ComponentSet((self.SAMPLEPOINT_SET,))
        blo_list = [MHEBlock.ACTUAL_MEASUREMENT_BLOCK, 
                    MHEBlock.MEASUREMENT_ERROR_BLOCK, 
                    MHEBlock.MODEL_DISTURBANCE_BLOCK]
        cate_list = [VC.ACTUALMEASUREMENT, 
                      VC.MEASUREMENTERROR, 
                      VC.MODELDISTURBANCE]
        MCTM = MHE_CATEGORY_TYPE_MAP
        for bloitem, category in zip(blo_list, cate_list):
            ctype = MCTM[category]
            newvars = []
            for i in bloitem:
                varlist = list(bloitem[i].component_objects(Var))
                assert len(varlist) == 1
                newvars.append(Reference(varlist[0], ctype=ctype))
            mhe_var_category_dict[category] = newvars

        CATEGORY_TYPE_MAP[VC.ACTUALMEASUREMENT] = ActualMeasurementVar
        CATEGORY_TYPE_MAP[VC.MEASUREMENTERROR] = MeasurementErrorVar
        CATEGORY_TYPE_MAP[VC.MODELDISTURBANCE] = ModelDisturbanceVar
        #super(_EstimatorBlockData, self)._construct()
        # All we need from the base class constructor is the following lines:
        self._add_category_blocks(mhe_var_category_dict)
        self._add_category_references(mhe_var_category_dict)

        self.category_dict.update(mhe_var_category_dict)
        self.categories.update(mhe_var_category_dict)

        self._add_mea_moddis_componentmap()

        self._add_measurement_constraint()
        self._add_disturbance_to_differential_cons()
        # at this stage, we don't need to add new constraints to con_categ_dict

        # Pop the mhe categories. This is important if any plant or controller is
        # defined afterward.
        CATEGORY_TYPE_MAP.pop(VC.ACTUALMEASUREMENT)
        CATEGORY_TYPE_MAP.pop(VC.MEASUREMENTERROR)
        CATEGORY_TYPE_MAP.pop(VC.MODELDISTURBANCE)

        #Set initial values for actual measurements
        # for ind in self.MEASUREMENT_SET:
        #     init_mea_block = self.MEASUREMENT_BLOCK[ind]
        #     init_val = {idx: 0.0 if init_mea_block.var[idx].value is None else init_mea_block.var[idx].value 
        #                 for idx in self.SAMPLEPOINT_SET} #Not sure why it doesn't work when the value is None
        #     self.ACTUALMEASUREMENT_BLOCK[ind].var.set_values(init_val)

    def _use_user_given_var_categ_dict(self, inputs, measurements):
        unref_category_dict = self._category_dict
        if (VC.INPUT not in unref_category_dict and inputs is not None):
            unref_category_dict[VC.INPUT] = inputs
        if (VC.MEASUREMENT not in unref_category_dict and measurements is not None):
            unref_category_dict[VC.MEASUREMENT] = measurements
        self.dae_vars = []
        for categ, varlist in unref_category_dict.items():
            if categ not in [VC.MEASUREMENT, VC.UNUSED]:
                # Assume that measurements are duplicates
                self.dae_vars.extend(varlist)
        
        return unref_category_dict

        
    def _categorize_var_con_for_MHE(self, mod, time, inputs, measurements):
        '''This function category variables and constraints for MHE.
        1. Flatten vars and constraints
        2. Convert the given _GeneralVarData in inputs and measurements to IndexedVar
        3. Categorize variables and constraints
        '''
        
        t0 = time.first()
        
        scalar_vars, dae_vars = flatten_dae_components(mod, time, ctype=Var)
        self.scalar_vars = scalar_vars
        self.dae_vars = dae_vars
        
        scalar_cons, dae_cons = flatten_dae_components(mod, time, ctype = Constraint)
        
        #change items in inputs and measurements to IndexedVar with only time index
        dae_map = ComponentMap((var[t0], var) for var in dae_vars)
        t0_vardata = list(dae_map.keys())
        input_set = ComponentSet(inputs)
        measurement_set = ComponentSet(measurements)
        input_vars = []
        measurement_vars = []
        #This loop is important in order to maintain the same order of 
        #measurements and inputs between MHE and plant
        for var0 in t0_vardata: 
            if var0 in input_set:
                time_slice = dae_map[var0]
                input_vars.append(time_slice)
            if var0 in measurement_set:
                time_slice = dae_map[var0]
                measurement_vars.append(time_slice)

        unref_category_dict, unref_con_category_dict = categorize_dae_variables_and_constraints(
                                                                                    mod,
                                                                                    dae_vars,
                                                                                    dae_cons,
                                                                                    time,
                                                                                    input_vars = input_vars,
                                                                                    )
        
        unref_category_dict[VC.MEASUREMENT] = measurement_vars
        
        return unref_category_dict, unref_con_category_dict
    
    
    def reference_var_category_dict(self, unref_category_dict):
        category_dict = {
                        category: [
                            Reference(ref.referent, ctype=ctype)
                            for ref in unref_category_dict[category]
                            ]
                        for category, ctype in CATEGORY_TYPE_MAP.items()
                        if category in unref_category_dict #skip disturbance and unused
                        }
        
        self.category_dict.update(category_dict)
        

    def _add_sample_point_set(self):
        '''
        Measurement errors, actual measurements, and model disturbances are indexed by
        sample points. 
        '''
        set_name = "SAMPLEPOINT_SET"
        self.add_component(set_name, ContinuousSet(initialize = self.sample_points))
        

    def _add_actual_measurement_param(self):
        """This function creates a indexed block and "fixed" variables to allocate 
        actual measurements that come from the plant. 
        """
        MHEBlock = self.MHE_VARS_CONS_BLOCK
        
        block_name = "ACTUAL_MEASUREMENT_BLOCK"
        mea_set = MHEBlock.MEASUREMENT_SET
        actmea_block = Block(mea_set)
        MHEBlock.add_component(block_name, actmea_block)

        sps_set = self.SAMPLEPOINT_SET       
        for i in mea_set:
            var_name = "actual_measurement"
            actmea_block[i].add_component(var_name, Var(sps_set, initialize = 0.0))
            #Variables in this block should always be fixed!
            actmea_block[i].find_component(var_name).fix() 
            
    
    def _add_measurement_error(self):
        """This function creates a indexed block for measurement errors.
        """
        MHEBlock = self.MHE_VARS_CONS_BLOCK
        
        block_name = "MEASUREMENT_ERROR_BLOCK"
        mea_set = MHEBlock.MEASUREMENT_SET
        meaerr_block = Block(mea_set)
        MHEBlock.add_component(block_name, meaerr_block)
        
        sps_set = self.SAMPLEPOINT_SET
        for i in mea_set:
            var_name = "measurement_error"
            meaerr_block[i].add_component(var_name, Var(sps_set, initialize = 0.0))
            
        
    def _add_model_disturbance(self):
        """This function creates a indexed block for model disturbances.
        The order of disturbances is the same as differential vars and derivative vars.
        """
        MHEBlock = self.MHE_VARS_CONS_BLOCK
        
        block_name = "MODEL_DISTURBANCE_BLOCK"
        diffvar_set = MHEBlock.DIFFERENTIAL_SET
        moddis_block = Block(diffvar_set)
        MHEBlock.add_component(block_name, moddis_block)
        
        sps_set = self.SAMPLEPOINT_SET
        t0 = self.time.first()
        for i in diffvar_set:
            var_name = "model_disturbance"
            moddis_block[i].add_component(var_name, Var(sps_set, initialize = 0.0))
            #Fix model disturbance at t = 0 as 0.0
            moddis_block[i].find_component(var_name)[t0].fix(0.0) 
            
            
    def _add_MHE_vars_to_category_dict(self):        
        MHEBlock = self.MHE_VARS_CONS_BLOCK
        # target_set = ComponentSet((self.SAMPLEPOINT_SET,))        
        blo_list = [MHEBlock.ACTUAL_MEASUREMENT_BLOCK, 
                    MHEBlock.MEASUREMENT_ERROR_BLOCK, 
                    MHEBlock.MODEL_DISTURBANCE_BLOCK]
        cate_list = [VC.ACTUALMEASUREMENT, 
                      VC.MEASUREMENTERROR, 
                      VC.MODELDISTURBANCE]
        MCTM = MHE_CATEGORY_TYPE_MAP
        for bloitem, category in zip(blo_list, cate_list):
            ctype = MCTM[category]
            newvars = []
            for i in bloitem:                        
                varlist = list(bloitem[i].component_objects(Var))
                assert len(varlist) == 1
                newvars.append(Reference(varlist[0], ctype=ctype))
            self.category_dict[category] = newvars


    def _add_mea_moddis_componentmap(self):
        '''Add component map for:
            1. diffvar @ t0 --> model disturbance var
            2. derivar @ t0 --> model disturbance var
            3. meavar @ t0 --> measurement error var
            4. meavar @ t0 --> actual measurement var
        '''
        
        t0 = self.time.first()
        
        differential_block = self.find_component("DIFFERENTIAL_BLOCK")
        moddis_block = self.find_component("MODELDISTURBANCE_BLOCK")
        diffvar_set = self.DIFFERENTIAL_SET
        self.diffvar_map_moddis = ComponentMap((differential_block[ind].var[t0], moddis_block[ind].var)
                                                for ind in diffvar_set)
        
        derivative_block = self.find_component("DERIVATIVE_BLOCK")
        derivar_set = self.DERIVATIVE_SET
        self.derivar_map_moddis = ComponentMap((derivative_block[ind].var[t0], moddis_block[ind].var)
                                                for ind in derivar_set)
        #Note that order of self.DIFFERENTIAL_SET and the order of self.DERIVATIVE_SET are corresponding.
                                                
              
        measurement_block = self.find_component("MEASUREMENT_BLOCK")
        meaerr_block = self.find_component("MEASUREMENTERROR_BLOCK")
        actmea_block = self.find_component("ACTUALMEASUREMENT_BLOCK")
        mea_set = self.MEASUREMENT_SET
        self.meavar_map_meaerr = ComponentMap((measurement_block[ind].var[t0], meaerr_block[ind].var)
                                                for ind in mea_set)
        self.meavar_map_actmea = ComponentMap((measurement_block[ind].var[t0], actmea_block[ind].var)
                                                for ind in mea_set)
        
        
    def _add_measurement_constraint(self):
        '''
        This function creates a indexed block, containing measurement error equations.
        '''
        
        MHEBlock = self.MHE_VARS_CONS_BLOCK
        block_name = "MEASUREMENT_CONSTRAINT_BLOCK" 
        mea_set = MHEBlock.MEASUREMENT_SET
        meacon_block = Block(mea_set)
        MHEBlock.add_component(block_name, meacon_block)
        
        con_name = "con_mea_err" #actual_mea = measured_state + mea_err
        time = self.time
        sample_points = self.sample_points
        for i in mea_set:
            def _con_mea_err(m, s):
                # Measurement equation is only defined at sample points, but 
                # it needs to be constructed with time set. 
                # Otherwise, other functions wouldn't work. 
                if s not in sample_points:
                    return Constraint.Skip
                else:
                    ind = m.index()
                    top_parent = m.parent_block().parent_block()
                    t0 = top_parent.time.first()
                    measured_stat = top_parent.MEASUREMENT_BLOCK[ind].var
                    actual_mea = top_parent.meavar_map_actmea[measured_stat[t0]]
                    mea_err = top_parent.meavar_map_meaerr[measured_stat[t0]]
                    return actual_mea[s] == measured_stat[s] + mea_err[s]
            meacon_block[i].add_component(con_name, Constraint(time, 
                                                               rule = _con_mea_err))
    
        
    def _add_disturbance_to_differential_cons(self):
        """
        This function creates a indexed block, containing model differential equations
        with model disturbances.
        """
        MHEBlock = self.MHE_VARS_CONS_BLOCK
        block_name = "DISTURBED_DIFFERENTIAL_CONSTRAINT_BLOCK"
        diff_equs = self.con_category_dict[CC.DIFFERENTIAL]
        n_diff_equs = len(diff_equs)
        assert n_diff_equs == len(self.DIFFERENTIAL_SET) #This should be correct!
        ddc_block = Block(range(n_diff_equs))
        MHEBlock.add_component(block_name, ddc_block)
        
        #Because model disturbances is indexed by sample time, this function finds
        #the corresponding sample time from a given time point.
        def curr_sample_point(time_pt, sample_points):
            for j in sample_points:
                if j >= time_pt:
                    return j
                
        time = self.time
        sample_points = self.sample_points
        moddistlist = self.category_dict[VC.MODELDISTURBANCE]
        for i in range(n_diff_equs):
            curr_difeq = diff_equs[i]
            mod_dist = moddistlist[i]
            con_name = "disturbed_diff_con"
            def _dd_rule(m, tp):
                curr_sampt = curr_sample_point(tp, sample_points)
                newbody = (curr_difeq[tp].body - mod_dist[curr_sampt])
                #Differential constraints are confirmed as equalities in categorize_dae_variables_and_constraints
                #I just set it equal to the original upper bound for convience.
                newexpr = newbody == curr_difeq[tp].upper.value 
                return newexpr
            _dd_con = Constraint(time, rule=_dd_rule)
            ddc_block[i].add_component(con_name, _dd_con)
            
            #Deactivate the original/undisturbed differential equs
            curr_difeq.deactivate() 
            
    def add_steady_state_objective(self, desired_ss, ss_weights):
        '''
        Add an objective function for solving a steady state, which is used to 
        initialize past information (measurements, states, etc.)

        Parameters
        ----------
        desired_ss : List of vardata, value tuples describing the
                     desired steady stsates of these specified variables.
        ss_weights : List of vardata, value tuples describing the
                     weightss of these specified variables.
        '''
        
        vardata_map = self.vardata_map
        for vardata, weight in ss_weights:
            nmpc_var = vardata_map[vardata]
            nmpc_var.weight = weight

        weight_vector = []
        for vardata, sp in desired_ss:
            nmpc_var = vardata_map[vardata]
            if nmpc_var.weight is None:
                self.logger.warning('Weight not supplied for %s' % var.name)
                nmpc_var.weight = 1.0
            weight_vector.append(nmpc_var.weight)

        obj_expr = sum(
            weight_vector[i]*(var - sp)**2 for
            i, (var, sp) in enumerate(desired_ss))
        self.steadystate_objective = Objective(expr=obj_expr)
        
    def solve_steady_state(self, solver):
        '''
        This function solve for a steady state, which will be used to initialize
        past information.
        '''
        
        model = self.mod
        time = self.time
        t0 = time.first()
        
        was_originally_active = ComponentMap([(comp, comp.active) for comp in 
                model.component_data_objects((Constraint, Block))])
        non_initial_time = list(time)[1:]
        deactivated = deactivate_model_at(
                model,
                time,
                non_initial_time,
                allow_skip=True,
                suppress_warnings=True,
                )
        
        self.vectors.differential[:, t0].unfix()
        self.vectors.input[:, t0].unfix() #dof
        self.vectors.derivative[:, t0].fix(0.)
        # Solving for steady state doesn't need MHE block.
        self.MHE_VARS_CONS_BLOCK.deactivate()
        # Activate the original/undisturbed differential equations at t0
        for indexcon in self.con_category_dict[CC.DIFFERENTIAL]:
            indexcon[t0].activate()
            
        self.steadystate_objective.activate()
        
        dof = degrees_of_freedom(model)
        # This should be True for solving a steady state.
        assert dof == len(self.INPUT_SET)
        
        results = solver.solve(self, tee=True)
        if results.solver.termination_condition == TerminationCondition.optimal:
            pass
        else:
            msg = 'Failed to solve for full steady state values'
            raise RuntimeError(msg)

        self.steadystate_objective.deactivate()
        
        for indexcon in self.con_category_dict[CC.DIFFERENTIAL]:
            indexcon[t0].deactivate()
        self.MHE_VARS_CONS_BLOCK.activate()
        self.vectors.derivative[:, t0].unfix()
        self.vectors.input[:, t0].fix()
        self.vectors.differential[:, t0].fix()
        
        for t, complist in deactivated.items():
            for comp in complist:
                if was_originally_active[comp]:
                    comp.activate()
                    
    def initialize_actualmeasurements_at_t0(self):
        t0 = self.sample_points[0]
        for ind in self.MEASUREMENT_SET:
            actmea_block = self.ACTUALMEASUREMENT_BLOCK[ind]
            mea_block = self.MEASUREMENT_BLOCK[ind]
            actmea_block.var[0].set_value(mea_block.var[0].value)  
        
    def initialize_past_info_with_steady_state(self, 
                                               desired_ss, 
                                               ss_weights, 
                                               solver):
        '''
        Before running MHE, we need to set up past information (measurements, states, etc.). 
        Here we solve for a steady state and set it as past information.
        
        Parameters
        ----------
        desired_ss : List of vardata, value tuples describing the
                     desired steady stsates of these specified variables.
        ss_weights : List of vardata, value tuples describing the
                     weightss of these specified variables.

        '''
        self.add_steady_state_objective(desired_ss, ss_weights)
        self.solve_steady_state(solver)
        
        self.initialize_actualmeasurements_at_t0()
        self.initialize_to_initial_conditions(ctype=(DiffVar, AlgVar, DerivVar, InputVar,))

          
    def add_noise_minimize_objective(self,
                                     model_disturbance_weights,
                                     measurement_noise_weights,
                                     givenform = "weight"):
        '''
        Set up the objective function for MHE. The given form can be either 'weight'
        or 'variance'.
        
        Parameters
        ----------        
        model_disturbance_weights: A list of vardata-value tuples describing 
                                    the weight for model disturbance to be 
                                    given to the objective term containing each variable.
        measurement noise_weights: A list of vardata-value tuples describing 
                                    the weight for measurement error to be 
                                    given to the objective term containing each variable.
        givenform: The form of given lists. It should be either "weight" (default)
                    or "variance".
        '''
        
        if givenform not in ["weight", "variance"]:
            raise RuntimeError("Wrong argument 'givenform' is given. "
                   "Please assign either 'weight' or 'variance'.")
        
        diff_id_list = [id(var[0]) for var in self.differential_vars]
        for var, val in model_disturbance_weights:  
            #Check whether the given variable is classified under diffvar
            if id(var) not in diff_id_list:
                raise RuntimeError(var.name, 
                                   " is not a differential variable.")
                
            if givenform == "weight":
                self.diffvar_map_moddis[var].weight = val
            elif givenform == "variance":
                self.diffvar_map_moddis[var].weight = 1./val

        mea_id_list = [id(var[0]) for var in self.measurement_vars]
        for var, val in measurement_noise_weights:
            #Check whether the given variable is declared as a measurement before
            if id(var) not in mea_id_list:
                raise RuntimeError(var.name, 
                                   " is not declared as a measurement.")
                
            if givenform == "weight":
                self.meavar_map_meaerr[var].weight = val
            elif givenform == "variance":
                self.meavar_map_meaerr[var].weight = 1./val
        
        moddis_block = self.MODELDISTURBANCE_BLOCK
        moddis_list = [moddis_block[ind].var for ind in self.DIFFERENTIAL_SET]
        
        meaerr_block = self.MEASUREMENTERROR_BLOCK
        meaerr_list = [meaerr_block[ind].var for ind in self.MEASUREMENT_SET]
        
        wQw = sum(
            moddis.weight * moddis[sampt]**2
            for moddis in moddis_list
            for sampt in self.sample_points
            if sampt != self.time.first() #Skip model disturbnace at time.first()
            )
        
        vRv = sum(
            meaerr.weight * meaerr[sampt]**2
            for meaerr in meaerr_list
            for sampt in self.sample_points #Here should include the time.first()
            )
        
        self.noise_minimize_objective = Objective(expr = (wQw) + (vRv))  
        
    def check_var_con_dof(self, skip_dof_check = False):
        '''This function does the final check for fixed and unfixed variables
        as well as the deactivated differential equations. Then check the degree of freedom.
        
        Parameters
        ----------
        skip_dof_check: SKip to check the degree of freedom if it's True. 
                        Default is False.
        '''
        
        self.vectors.input[...].fix()
        self.vectors.differential[...].unfix()
        self.vectors.modeldisturbance[:,0].fix(0.0)
        
        for diff_eq in self.con_category_dict[CC.DIFFERENTIAL]:
            diff_eq.deactivate()
        
        n_diffvars = len(self.DIFFERENTIAL_SET.ordered_data())
        n_steps_in_horizon = len(self.sample_points) - 1
        correct_dof = n_diffvars * (n_steps_in_horizon + 1)
        
        if not skip_dof_check:
            dof = degrees_of_freedom(self)
            assert dof == correct_dof  
        
    def load_inputs_for_MHE(self, inputs):
        last_sampt = self.sample_points[-1]
        secondlast_sampt = self.sample_points[-2]
        time_list = [tp for tp in self.time if tp > secondlast_sampt 
                                                 and tp <= last_sampt]
        
        for var, val in zip(self.INPUT_BLOCK[:].var, inputs):
            for tind in time_list:
                var[tind].set_value(val)
                
    def generate_estimates_at_time(self, t):
        return [val for val in self.vectors.differential[:, t].value]
                    

class EstimatorBlock(DynamicBlock):
    """ This is a user-facing class to be instantiated when one
    wants to work with a block capable of the methods of
    `_EstimatorBlockData`.
    """
    _ComponentDataClass = _EstimatorBlockData

    def __new__(cls, *args, **kwds):
        # Decide what class to allocate
        if cls != EstimatorBlock:
            target_cls = cls
        elif not args or (args[0] is UnindexedComponent_set and len(args) == 1):
            target_cls = SimpleEstimatorBlock
        else:
            target_cls = IndexedEstimatorBlock
        return super(EstimatorBlock, cls).__new__(target_cls)


class SimpleEstimatorBlock(_EstimatorBlockData, EstimatorBlock):
    def __init__(self, *args, **kwds):
        _EstimatorBlockData.__init__(self, component=self)
        EstimatorBlock.__init__(self, *args, **kwds)

    # Pick up the display() from Block and not BlockData
    display = EstimatorBlock.display


class IndexedEstimatorBlock(EstimatorBlock):

    def __init__(self, *args, **kwargs):
        EstimatorBlock.__init__(self, *args, **kwargs)
