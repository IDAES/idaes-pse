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
from idaes.apps.caprese.common.config import VariableCategory as VC
from idaes.apps.caprese.categorize import (
        categorize_dae_variables,
        _identify_derivative_if_differential,
        categorize_dae_variables_and_constraints,
        )
from idaes.apps.caprese.nmpc_var import (
        _NmpcVector,
        NmpcVar,
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
        )
from pyomo.core.base.block import _BlockData
from pyomo.common.collections import ComponentMap, ComponentSet
from pyomo.dae.flatten import (
        flatten_components_along_sets,
        flatten_dae_components,
        )
from pyomo.core.base.indexed_component import UnindexedComponent_set
from pyomo.dae.contset import ContinuousSet


MHE_CATEGORY_TYPE_MAP = {
    VC.ACTUALMEASUREMENT: ActualMeasurementVar,
    VC.MEASUREMENTERROR: MeasurementErrorVar,
    VC.MODELDISTURBANCE: ModelDisturbanceVar,
                        }

class _EstimatorBlockData(_DynamicBlockData):
    """ This class adds methods useful for working with dynamic
    models to be used by a estimator. These include methods for
    """

    def _construct(self):
        if self._sample_time is None:
            raise RuntimeError("MHE needs the sample time to be provided!")
        else:
            sample_time = self._sample_time
            self.set_sample_time(sample_time)
            
        mod = self.mod
        time = self.time
        
        try:
            inputs = self._inputs
        except AttributeError:
            inputs = self._inputs = None
        try:
            measurements = self._measurements
        except AttributeError:
            measurements = self._measurements = None
        
        scalar_vars, dae_vars = flatten_dae_components(
                mod,
                time,
                ctype=Var,
                )
        self.scalar_vars = scalar_vars
        self.dae_vars = dae_vars
        
        category_dict = categorize_dae_variables(
                dae_vars,
                time,
                inputs,
                measurements=measurements,
                )
        self.category_dict = category_dict
            
        scalar_cons, dae_cons = flatten_dae_components(mod, time, ctype = Constraint)
        #need to add reference if use the following var_cactegory_dict
        not_use_category_dict, con_category_dict = categorize_dae_variables_and_constraints(
                                                                                mod,
                                                                                dae_vars,
                                                                                dae_cons,
                                                                                time,
                                                                                )
        self.con_category_dict = con_category_dict  
        
        #set a flag so that classification won't do again in _DynamicBlockData
        self.already_categorized_for_MHE = True
        

        self._add_sample_point_set()
        self.MHE_VARS_CONS_BLOCK = Block()
        
        n_measurement = len(category_dict[VC.MEASUREMENT])
        n_diffvar_con = len(category_dict[VC.DIFFERENTIAL])
        self.MHE_VARS_CONS_BLOCK.MEASUREMENT_SET = Set(initialize = range(n_measurement))
        self.MHE_VARS_CONS_BLOCK.DIFFERENTIAL_SET = Set(initialize = range(n_diffvar_con))
        
        self._add_actual_measurement_param()
        self._add_measurement_error()
        self._add_model_disturbance()
        
        self._add_MHE_vars_to_category_dict()
        
        super(_EstimatorBlockData, self)._construct()
        if VC.ACTUALMEASUREMENT in self.category_dict:
            self.actualmeasurement_vars = category_dict[VC.ACTUALMEASUREMENT]
        if VC.MEASUREMENTERROR in self.category_dict:
            self.measurementerror_vars = category_dict[VC.MEASUREMENTERROR]
        if VC.MODELDISTURBANCE in self.category_dict:
            self.modeldisturbance_vars = category_dict[VC.MODELDISTURBANCE]
        
        self._add_mea_moddis_componentmap()
        
        # Add the new constraint to the con_category_dict?
        self._add_measurement_constraint()
        self._add_disturbance_to_differential_cons()
        
        #Set initial values for actual measurements
        for ind in self.MEASUREMENT_SET:
            init_mea_block = self.MEASUREMENT_BLOCK[ind]
            init_val = {idx: 0.0 if init_mea_block.var[idx].value is None else init_mea_block.var[idx].value 
                        for idx in self.SAMPLEPOINT_SET} #Not sure why it doesn't work when the value is None
            self.ACTUALMEASUREMENT_BLOCK[ind].var.set_values(init_val)

    def _add_sample_point_set(self):
        set_name = "SAMPLEPOINT_SET"
        self.add_component(set_name, ContinuousSet(initialize = self.sample_points))
        

    def _add_actual_measurement_param(self):
        """This function creates a indexed block and "fixed" variables to allocate 
        actual measurement measured from the plant. 
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
            actmea_block[i].find_component(var_name).fix() #Variables in this block should always be fixed!
    
    def _add_measurement_error(self):
        """This function creates a indexed block, including measurement errors and 
        measurement constraints: actual mea = measured state + mea_err
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
        """This function creates a indexed block, including model disturbances.
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
            moddis_block[i].find_component(var_name)[t0].fix(0.0) #fix model disturbance at t = 0 as 0.0
            
    def _add_MHE_vars_to_category_dict(self):        
        MHEBlock = self.MHE_VARS_CONS_BLOCK
        target_set = ComponentSet((self.SAMPLEPOINT_SET,))        
        blo_list = [MHEBlock.ACTUAL_MEASUREMENT_BLOCK, 
                    MHEBlock.MEASUREMENT_ERROR_BLOCK, 
                    MHEBlock.MODEL_DISTURBANCE_BLOCK]
        cate_list = [VC.ACTUALMEASUREMENT, 
                     VC.MEASUREMENTERROR, 
                     VC.MODELDISTURBANCE]
        MHE_category_list_map = {}
        for bloitem, cateitem in zip(blo_list, cate_list):
            save_list = []
            for ind in bloitem.index_set():
                sets_list, comps_list = flatten_components_along_sets(bloitem[ind], 
                                                                      target_set, 
                                                                      Var)
                save_list.append(comps_list[0][0])
            MHE_category_list_map[cateitem] = save_list
            
        for category, ctype in MHE_CATEGORY_TYPE_MAP.items():
            self.category_dict[category] = [Reference(var, ctype=ctype) 
                                            for var in MHE_category_list_map[category]]
        
            
    def _add_mea_moddis_componentmap(self):
        
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
                if s not in sample_points:
                    return Constraint.Skip
                else:
                    ind = m.index()
                    top_parent = m.parent_block().parent_block()
                    measured_stat = top_parent.MEASUREMENT_BLOCK[ind].var
                    t0 = top_parent.time.first()
                    actual_mea = top_parent.meavar_map_actmea[measured_stat[t0]]
                    mea_err = top_parent.meavar_map_meaerr[measured_stat[t0]]
                    return actual_mea[s] == measured_stat[s] + mea_err[s]
            meacon_block[i].add_component(con_name, Constraint(time, 
                                                               rule = _con_mea_err))
    
        
    def _add_disturbance_to_differential_cons(self):
        """
        This function creates a indexed block, including model differential equations
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

        #Now, var_category and con_category don't come from the same function.
        #(var_category is from categorize_dae_variables; con_category is from categorize_dae_variables_and_constraints)
        #Thus, the order of differential vars may be different from the order of differential equations.
        #(Note that the order of differential vars, derivative vars, and model disturbances are corresponding.)
        #Use the following function to get the correct model disturbance.
        def find_correspondsing_model_disturbance(curr_difeq, time, derivar_map_moddis):
            t0 = time.first()
            is_diff, target_deriv = _identify_derivative_if_differential(curr_difeq[t0], time)
            if not is_diff:
                raise RuntimeError(
                "Fail to find the corresponding model disturbance from a differential equation."
                )
            return derivar_map_moddis[target_deriv]
        #if var_category and con_category are all from "categorize_dae_variables_and_constraints",
        #"find_correspondsing_model_disturbance" can be removed!
                
        time = self.time
        sample_points = self.sample_points
        derivativelist = self.category_dict[VC.DERIVATIVE]
        for i in range(n_diff_equs):
            curr_difeq = diff_equs[i]
            # Get corresponding disturbance from index i
            # And gettting the time point is easy.
            mod_dist = find_correspondsing_model_disturbance(curr_difeq, 
                                                             time, 
                                                             self.derivar_map_moddis)
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
                
            curr_difeq.deactivate() #deactivate the original/undisturbed differential equs
            
    def add_noise_minimize_objective(self,
                                     model_disturbance_weights,
                                     measurement_noise_weights):
        
        #give weights or variances?
        for var, val in model_disturbance_weights:
            self.diffvar_map_moddis[var].weight = val
            
        for var, val in measurement_noise_weights:
            self.meavar_map_meaerr[var].weight = val
            
        #check if given measurement is not declared as measurement before
        
        moddis_block = self.MODELDISTURBANCE_BLOCK
        moddis_list = [moddis_block[ind].var for ind in self.DIFFERENTIAL_SET]
        
        meaerr_block = self.MEASUREMENTERROR_BLOCK
        meaerr_list = [meaerr_block[ind].var for ind in self.MEASUREMENT_SET]
        
        wQw = sum(
            moddis.weight * moddis[sampt]**2
            for moddis in moddis_list
            for sampt in self.sample_points
            if sampt != self.time.first()
            )
        
        vRv = sum(
            meaerr.weight * meaerr[sampt]**2
            for meaerr in meaerr_list
            for sampt in self.sample_points #here should include the time.first()
            )
        
        self.noise_minimize_objective = Objective(expr = (wQw) + (vRv))  
        
    def check_var_con_dof(self, skip_dof_check = False):
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
     
    def initialize_actualmeasurements_at_t0(self):
        t0 = self.sample_points[0]
        for ind in self.MEASUREMENT_SET:
            actmea_block = self.ACTUALMEASUREMENT_BLOCK[ind]
            mea_block = self.MEASUREMENT_BLOCK[ind]
            actmea_block.var[0].set_value(mea_block.var[0].value)    
        
    def load_inputs_for_MHE(self, inputs):
        
        last_sampt = self.sample_points[-1]
        secondlast_sampt = self.sample_points[-2]
        time_list = [tp for tp in self.time if tp > secondlast_sampt 
                                                 and tp <= last_sampt]
        
        for var, val in zip(self.INPUT_BLOCK[:].var, inputs):
            for tind in time_list:
                var[tind].set_value(val)
                    

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
