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
# from idaes.apps.caprese.util import initialize_by_element_in_range
# from idaes.apps.caprese.common.config import (
#         ControlPenaltyType,
#         )
from idaes.apps.caprese.common.config import VariableCategory as VC
from idaes.apps.caprese.categorize import (
        # categorize_dae_variables,
        _identify_derivative_if_differential,
        # CATEGORY_TYPE_MAP,
        )
# from idaes.apps.caprese.nmpc_var import (
#         NmpcVar,
#         DiffVar,
#         AlgVar,
#         InputVar,
#         DerivVar,
#         FixedVar,
#         MeasuredVar,
#         )
from idaes.apps.caprese.dynamic_block import (
        _DynamicBlockData,
        IndexedDynamicBlock,
        DynamicBlock,
        )
# from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.apps.caprese.common.config import (
        VariableCategory as VC,
        ConstraintCategory as CC,
        )

from pyomo.environ import (
        Var,
#         Objective,
#         TerminationCondition,
        Constraint,
        Block,
        Reference,
        Set,
        )
from pyomo.core.base.block import _BlockData
# from pyomo.common.collections import ComponentMap
# from pyomo.core.base.range import remainder
# from pyomo.dae.set_utils import deactivate_model_at
# from pyomo.dae.flatten import flatten_dae_components
from pyomo.core.base.indexed_component import UnindexedComponent_set

import idaes.apps.caprese.nmpc_var as nmpc_var
# from idaes.apps.caprese.nmpc_var import (
#         NmpcVar,
#         )


class _EstimatorBlockData(_DynamicBlockData):
    """ This class adds methods useful for working with dynamic
    models to be used by a estimator. These include methods for
    """
    
    def empty_now(self):
        print("dummy fun = ", 123)
        
    def _add_sample_point_set(self):
        set_name = "SAMPLEPOINT_SET"
        self.add_component(set_name, Set(initialize = self.sample_points))
        
    def _add_actual_measurement_param(self):
        """This function creates a indexed block and "fixed" variables to allocate 
        actual measurement measured from the plant. 
        """
        block_name = "ACTUALMEASUREMENT_BLOCK"
        mea_set = self.MEASUREMENT_SET
        actmea_block = Block(mea_set)
        self.add_component(block_name, actmea_block)
        # meaerr_block.deactivate() #keep it activate? Yes!
        sample_points = self.sample_points
        sps_set = self.SAMPLEPOINT_SET
        for i in mea_set:
            var_name = "var"
            init_mea_block = self.MEASUREMENT_BLOCK[i]
            init_val = {idx: init_mea_block.var[idx].value for idx in sample_points}
            actmea_block[i].add_component(var_name, Var(sps_set, initialize = init_val))
            actmea_block[i].find_component(var_name).fix()
    
    def _add_measurement_error(self):
        """This function creates a indexed block, including measurement errors and 
        measurement constraints: actual mea = measured state + mea_err
        """
        block_name = "MEASUREMENTERROR_BLOCK"
        mea_set = self.MEASUREMENT_SET
        meaerr_block = Block(mea_set)
        self.add_component(block_name, meaerr_block)
        # meaerr_block.deactivate() #keep it activate? Yes!
        sample_points = self.sample_points
        sps_set = self.SAMPLEPOINT_SET
        for i in mea_set:
            var_name = "var"
            meaerr_block[i].add_component(var_name, Var(sps_set, initialize = 0.0))
            
            # con_name = "con_mea_err" #actual_mea = measured_state + mea_err
            # def _con_mea_err(m, s):
            #     ind = m.index()
            #     parent = m.parent_block()
            #     actual_mea = parent.ACTUALMEASUREMENT_BLOCK[ind].actual_mea
            #     measured_stat = parent.MEASUREMENT_BLOCK[ind].var
            #     mea_err = parent.MEASUREMENTERROR_BLOCK[ind].mea_err
            #     return actual_mea[s] == measured_stat[s] + mea_err[s]
            # meaerr_block[i].add_component(con_name, Constraint(sample_points, rule = _con_mea_err))
            
    def _add_measurement_constraint(self):
        block_name = "MEASUREMENT_CONSTRAINT_BLOCK" 
        mea_set = self.MEASUREMENT_SET
        meacon_block = Block(mea_set)
        self.add_component(block_name, meacon_block)
        
        con_name = "con_mea_err" #actual_mea = measured_state + mea_err
        time = self.time
        sample_points = self.sample_points
        for i in mea_set:
            def _con_mea_err(m, s):
                if s not in sample_points:
                    return Constraint.Skip
                else:
                    ind = m.index()
                    parent = m.parent_block()
                    actual_mea = parent.ACTUALMEASUREMENT_BLOCK[ind].var
                    measured_stat = parent.MEASUREMENT_BLOCK[ind].var
                    mea_err = parent.MEASUREMENTERROR_BLOCK[ind].var
                    return actual_mea[s] == measured_stat[s] + mea_err[s]
            meacon_block[i].add_component(con_name, Constraint(time, 
                                                               rule = _con_mea_err))
    
    def _add_model_disturbance(self):
        """This function creates a indexed block, including model disturbances.
        The order of disturbances is the same as differential vars and derivative vars.
        """
        block_name = "MODELDISTURBANCE_BLOCK"
        diffvar_set = self.DIFFERENTIAL_SET
        moddis_block = Block(diffvar_set)
        self.add_component(block_name, moddis_block)
        sample_points = self.sample_points
        sps_set = self.SAMPLEPOINT_SET
        for i in diffvar_set:
            var_name = "var"
            moddis_block[i].add_component(var_name, Var(sps_set, initialize = 0.0))
            moddis_block[i].find_component(var_name)[0].fix(0.0) #fix model disturbance at t = 0 as 0.0
        
    def _add_disturbance_to_differential_cons(self):
        """
        This function creates a indexed block, including model differential equations
        with model disturbances.
        """
        block_name = "DISTURBED_DIFFERENTIAL_CONSTRAINT_BLOCK"
        diff_equs = self.con_category_dict[CC.DIFFERENTIAL]
        n_diff_equs = len(diff_equs)
        assert n_diff_equs == len(self.DIFFERENTIAL_SET) #This should be correct!
        ddc_block = Block(range(n_diff_equs))
        self.add_component(block_name, ddc_block)
        
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
        def find_correspondsing_model_disturbance(curr_difeq, time, derivativelist):
            t0 = time.first()
            check_list = [id(var[t0]) for var in derivativelist]
            is_diff, target_deriv = _identify_derivative_if_differential(curr_difeq[t0], time)
            if not is_diff:
                raise RuntimeError(
                "Fail to find the corresponding model disturbance from a differential equation."
                )
            target_id = id(target_deriv)
            model_dist_block_ind = check_list.index(target_id)
            return model_dist_block_ind
        #if var_category and con_category are all from "categorize_dae_variables_and_constraints",
        #"find_correspondsing_model_disturbance" can be removed!
                
        time = self.time
        sample_points = self.sample_points
        derivativelist = self.category_dict[VC.DERIVATIVE]
        for i in range(n_diff_equs):
            curr_difeq = diff_equs[i]
            mod_dist_block_ind = find_correspondsing_model_disturbance(curr_difeq, 
                                                                       time, 
                                                                       derivativelist)
            con_name = "disturbed_diff_con"
            ddc_block[i].add_component(con_name, Constraint(time))
            #Here run the iterations to set the constraint expression one by one. Other better way to do that?
            for tp in time:
                curr_sampt = curr_sample_point(tp, sample_points)
                newbody = (curr_difeq[tp].body - 
                           self.MODELDISTURBANCE_BLOCK[mod_dist_block_ind].var[curr_sampt])
                #Differential constraints are confirmed as equalities in categorize_dae_variables_and_constraints
                #I just set it equal to the original upper bound for convience.
                newexpr = newbody == curr_difeq[tp].upper.value 
                con_attr = ddc_block[i].find_component(con_name)
                con_attr.add(tp, expr = newexpr)
            curr_difeq.deactivate()
            
    def _add_MHE_variables_into_vectors(self):
        MHE_CATEGORY_TYPE_MAP = {
            VC.ACTUALMEASUREMENT: nmpc_var.ActualMeasurementVar,
            VC.MEASUREMENTERROR: nmpc_var.MeasurementErrorVar,
            VC.MODELDISTURBANCE: nmpc_var.ModelDisturbanceVar,
                                }
        
        MHE_category_list_map = {
            VC.ACTUALMEASUREMENT: [self.ACTUALMEASUREMENT_BLOCK[i].var
                                   for i in self.MEASUREMENT_SET],
            VC.MEASUREMENTERROR: [self.MEASUREMENTERROR_BLOCK[i].var
                                  for i in self.MEASUREMENT_SET],
            VC.MODELDISTURBANCE: [self.MODELDISTURBANCE_BLOCK[i].var
                                  for i in self.DERIVATIVE_SET],
                                }
        
        MHE_category_dict = {
            category: [
                Reference(var, ctype=ctype)
                for var in MHE_category_list_map[category]
                ]
            for category, ctype in MHE_CATEGORY_TYPE_MAP.items()
                            }
        
        self.MHE_category_dict = MHE_category_dict
            
        for categ in MHE_category_dict:
            ctype = MHE_CATEGORY_TYPE_MAP.get(categ, nmpc_var.NmpcVar)
            block_name = self.get_category_block_name(categ)
            var_name = self._var_name
            block = getattr(self, block_name)
            _slice = block[:]
            _slice = getattr(_slice, var_name)[:]
            ref = Reference(_slice, ctype=nmpc_var._NmpcVector)
            self.vectors.add_component(
                    categ.name.lower(), # Lowercase of the enum name
                    ref,
                    )
        
        

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