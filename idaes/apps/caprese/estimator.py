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
# from idaes.apps.caprese.categorize import (
#         categorize_dae_variables,
#         CATEGORY_TYPE_MAP,
#         )
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

from pyomo.environ import (
        Var,
#         Objective,
#         TerminationCondition,
        Constraint,
        Block,
        )
from pyomo.core.base.block import _BlockData
# from pyomo.common.collections import ComponentMap
# from pyomo.core.base.range import remainder
# from pyomo.dae.set_utils import deactivate_model_at
# from pyomo.dae.flatten import flatten_dae_components
from pyomo.core.base.indexed_component import UnindexedComponent_set


class _EstimatorBlockData(_DynamicBlockData):
    """ This class adds methods useful for working with dynamic
    models to be used by a estimator. These include methods for
    """
    
    def empty_now(self):
        print("dummy fun = ", 123)
        
    def _add_actual_measurement_param(self):
        """This function creates a indexed block and "fixed" variables to allocate 
        actual measurement measured from the plant. 
        """
        block_name = "ACTUAL_MEASUREMENT_BLOCK"
        mea_set = self.MEASUREMENT_SET
        actmea_block = Block(mea_set)
        self.add_component(block_name, actmea_block)
        # meaerr_block.deactivate() #keep it activate? Yes!
        sample_points = self.sample_points
        for i in mea_set:
            var_name = "actual_mea"
            init_mea_block = self.MEASUREMENT_BLOCK[i]
            init_val = {idx: init_mea_block.var[idx].value for idx in sample_points}
            actmea_block[i].add_component(var_name, Var(sample_points, initialize = init_val))
            actmea_block[i].find_component(var_name).fix()
    
    def _add_measurement_error(self):
        """This function creates a indexed block, including measurement errors and 
        measurement constraints: actual mea = measured state + mea_err
        """
        block_name = "MEASUREMENT_ERROR_BLOCK"
        mea_set = self.MEASUREMENT_SET
        meaerr_block = Block(mea_set)
        self.add_component(block_name, meaerr_block)
        # meaerr_block.deactivate() #keep it activate? Yes!
        sample_points = self.sample_points
        for i in mea_set:
            var_name = "mea_err"
            meaerr_block[i].add_component(var_name, Var(sample_points, initialize = 0.0))
            
            con_name = "con_mea_err" #actual_mea = measured_state + mea_err
            def _con_mea_err(m, s):
                ind = m.index()
                parent = m.parent_block()
                actual_mea = parent.ACTUAL_MEASUREMENT_BLOCK[ind].actual_mea
                measured_stat = parent.MEASUREMENT_BLOCK[ind].var
                mea_err = parent.MEASUREMENT_ERROR_BLOCK[ind].mea_err
                return actual_mea[s] == measured_stat[s] + mea_err[s]
            meaerr_block[i].add_component(con_name, Constraint(sample_points, rule = _con_mea_err))
            
    # def _add_measurement_constraint(self):
    #     block_name = "MEASUREMENT_CONSTRAINT_BLOCK" #should be change after define MEA_ERR category
    #     mea_set = self.MEASUREMENT_SET #better to use get_category_set_name
    #     meacon_block = Block(mea_set)
    #     self.add_component(block_name, meacon_block)
        
    #     con_name = "con_mea_err" #actual_mea = measured_state + mea_err
    #     sample_points = self.sample_points
    #     for i in mea_set:
    #         def _con_mea_err(m, s):
    #             ind = m.index()
    #             parent = m.parent_block()
    #             actual_mea = parent.ACTUAL_MEASUREMENT_BLOCK[ind].actual_mea
    #             measured_stat = parent.MEASUREMENT_BLOCK[ind].var
    #             mea_err = parent.MEASUREMENT_ERROR_BLOCK[ind].mea_err
    #             return actual_mea[s] == measured_stat[s] + mea_err[s]
    #         meacon_block[i].add_component(con_name, Constraint(sample_points, rule = _con_mea_err))
    
    def _add_model_disturbance(self):
        """This function creates a indexed block, including model disturbances.
        The order of disturbances is the same as differential vars and derivative vars.
        """
        block_name = "MODEL_DISTURBANCE_BLOCK"
        diffvar_set = self.DIFFERENTIAL_SET
        moddis_block = Block(diffvar_set)
        self.add_component(block_name, moddis_block)
        sample_points = self.sample_points
        for i in diffvar_set:
            var_name = "mod_disturb"
            moddis_block[i].add_component(var_name, Var(sample_points, initialize = 0.0))
            moddis_block[i].find_component(var_name)[0].fix(0.0) #fix model disturbance at t = 0 as 0.0
        
        
        

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