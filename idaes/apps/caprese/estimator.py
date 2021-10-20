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

        # Update category dict with references to proper ctypes.
        # Not necessary for constraints currently.
        self.category_dict.update(
            self.reference_var_category_dict(self.category_dict)
        )

        # Assign local category_dict variable for convenience
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

        # These are new categories we need to add to the "vectors" block
        mhe_var_category_dict = self._get_mhe_var_category_dict()

        CATEGORY_TYPE_MAP[VC.ACTUALMEASUREMENT] = ActualMeasurementVar
        CATEGORY_TYPE_MAP[VC.MEASUREMENTERROR] = MeasurementErrorVar
        CATEGORY_TYPE_MAP[VC.MODELDISTURBANCE] = ModelDisturbanceVar

        # Add category blocks and references for MHE variables
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

    def reference_var_category_dict(self, unref_category_dict):
        category_dict = {
                        category: [
                            Reference(ref.referent, ctype=ctype)
                            for ref in unref_category_dict[category]
                            ]
                        for category, ctype in CATEGORY_TYPE_MAP.items()
                        if category in unref_category_dict #skip disturbance and unused
                        }
        return category_dict

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

    def _get_mhe_var_category_dict(self):
        MHEBlock = self.MHE_VARS_CONS_BLOCK
        # target_set = ComponentSet((self.SAMPLEPOINT_SET,))
        blo_list = [MHEBlock.ACTUAL_MEASUREMENT_BLOCK,
                    MHEBlock.MEASUREMENT_ERROR_BLOCK,
                    MHEBlock.MODEL_DISTURBANCE_BLOCK]
        cate_list = [VC.ACTUALMEASUREMENT,
                      VC.MEASUREMENTERROR,
                      VC.MODELDISTURBANCE]
        MCTM = MHE_CATEGORY_TYPE_MAP
        mhe_category_dict = {}
        for bloitem, category in zip(blo_list, cate_list):
            ctype = MCTM[category]
            newvars = []
            for i in bloitem:
                varlist = list(bloitem[i].component_objects(Var))
                assert len(varlist) == 1
                newvars.append(Reference(varlist[0], ctype=ctype))
            mhe_category_dict[category] = newvars
        return mhe_category_dict

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
        self.diffvar_map_moddis = ComponentMap(
            (differential_block[ind].var[t0], moddis_block[ind].var)
                                                for ind in diffvar_set)

        derivative_block = self.find_component("DERIVATIVE_BLOCK")
        derivar_set = self.DERIVATIVE_SET
        self.derivar_map_moddis = ComponentMap(
            (derivative_block[ind].var[t0], moddis_block[ind].var)
                                                for ind in derivar_set)
        # Note that order of self.DIFFERENTIAL_SET and the order of
        # self.DERIVATIVE_SET are corresponding.

        measurement_block = self.find_component("MEASUREMENT_BLOCK")
        meaerr_block = self.find_component("MEASUREMENTERROR_BLOCK")
        actmea_block = self.find_component("ACTUALMEASUREMENT_BLOCK")
        mea_set = self.MEASUREMENT_SET
        self.meavar_map_meaerr = ComponentMap(
            (measurement_block[ind].var[t0], meaerr_block[ind].var)
                                                for ind in mea_set)
        self.meavar_map_actmea = ComponentMap(
            (measurement_block[ind].var[t0], actmea_block[ind].var)
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
                # Differential constraints are confirmed as equalities 
                # in categorize_dae_variables_and_constraints.
                # I just set it equal to the original upper bound for convience.
                newexpr = newbody == curr_difeq[tp].upper.value
                return newexpr
            _dd_con = Constraint(time, rule=_dd_rule)
            ddc_block[i].add_component(con_name, _dd_con)

            #Deactivate the original/undisturbed differential equs
            curr_difeq.deactivate() 

    def initialize_actualmeasurements_at_t0(self):
        t0 = self.time.first()
        for ind in self.MEASUREMENT_SET:
            actmea_block = self.ACTUALMEASUREMENT_BLOCK[ind]
            mea_block = self.MEASUREMENT_BLOCK[ind]
            actmea_block.var[t0].set_value(mea_block.var[t0].value)

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
        # self.add_steady_state_objective(desired_ss, ss_weights)
        # self.solve_steady_state(solver)

        self.add_single_time_optimization_objective(desired_ss, ss_weights)
        self.solve_single_time_optimization(solver, 
                                            ic_type = "differential_var",
                                            require_steady = True,
                                            load_setpoints = False,
                                            restore_ic_input_after_solve = False,
                                            isMHE_block = True,)
        # This function solves an optimization problem, then uses the results
        # to initialize (a) measurement parameters and (b) variables in the
        # model.
        # We would rather (a) solve the optimization problem, then (b) get
        # scalar data ({cuid: value}) from the result.
        # We would then EITHER (a) use this scalar data to initialize
        # variables at initial points and copy the values forward OR
        # (b) convert this scalar data to time indexed data (effectively
        # doing the copying) and initialize the model from time indexed data.

        self.initialize_actualmeasurements_at_t0()
        self.initialize_to_initial_conditions(
            ctype=(DiffVar, AlgVar, DerivVar, InputVar, ActualMeasurementVar)
        )

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
        t0 = self.time.first()

        if givenform not in ["weight", "variance"]:
            raise RuntimeError("Wrong argument 'givenform' is given. "
                   "Please assign either 'weight' or 'variance'.")

        diff_id_set = {id(var[t0]) for var in self.differential_vars}
        for var, val in model_disturbance_weights:  
            #Check whether the given variable is classified under diffvar
            if id(var) not in diff_id_set:
                raise RuntimeError(var.name,
                                   " is not a differential variable.")

            if givenform == "weight":
                self.diffvar_map_moddis[var].weight = val
            elif givenform == "variance":
                self.diffvar_map_moddis[var].weight = 1./val

        mea_id_set = {id(var[t0]) for var in self.measurement_vars}
        for var, val in measurement_noise_weights:
            #Check whether the given variable is declared as a measurement before
            if id(var) not in mea_id_set:
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

    def load_inputs_into_last_sample(self, inputs):
        sample_points = self.sample_points
        # Accessing these entries is valid because sample_points
        # should always have length of at least two.
        t1 = sample_points[-2]
        t2 = sample_points[-1]
        time_subset = [t for t in self.time if t > t1 and t <= t2]
        self.inject_inputs(inputs, time_subset=time_subset)

    def generate_estimates_at_time(self, t):
        return [val for val in self.vectors.differential[:, t].value]

    def load_measurements(self, measured, timepoint=None):
        time = self.time
        t = timepoint if timepoint is not None else time.first()
        for var, val in zip(self.vectors.actualmeasurement[:, t], measured):
            var.fix(val)


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
