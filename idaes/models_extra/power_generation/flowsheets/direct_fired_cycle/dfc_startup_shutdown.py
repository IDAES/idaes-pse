"""
Contains functions that implement startup and shutdown constraints for the power cycle.

In this version, we implement the following and startup and shutdown procedure.

Cycle startup:
t = 0 - 1 => Ramp up fuel from 0 to 5.8518 kg/s. Power production = 0 MW
t = 1 - 2 => Ramp up fuel from 5.8515 to 10.36 kg/s. Power production = 0 MW
t = 2 - 3 => Cycle operates at P_min, (20% of the maximum load)
t = 3 - 4 => Normal operation of the plant. Cycle could ramp up to P_max

Cycle shutdown:
t = 0 - 1 => Cycle needs to ramp down to P_min
t = 1 - 2 => Cycle operates at P_min, (20% of the maximum load)
t = 2 - 3 => Ramp down fuel from 10.36 to 5.85 kg/s. Power production = 0 MW
t = 3 - 4 => Ramp down fuel from 5.8518 kg/s to 0. Power production = 0 MW
"""

from pyomo.environ import Constraint


DFC_OFFLOAD = 0.2


def dfc_startup_shutdown_constraints(m):
    # The startup and shutdown constraints will be added in a separate block 
    # for eCach of the processes.

    # Create an alias for the parent block. 
    pb = m.parent_block()
    set_period = pb.set_period
    num_time_periods = len(set_period)

    # Cycle can operate at time t, only if one of the following two is true
    #     - Cycle was in operation at time t - 1
    #     - Cycle startup was initiated at time t - 2
    @m.Constraint(set_period)
    def cycle_operation(blk, t):
        if t == 1:
            # Constraint is not applicable at t = 1. Initial condition is left free
            # But it can choose at most one of the three
            return (
                pb.period[t].fs.dfc.op_mode + pb.period[t].fs.dfc.startup +
                pb.period[t].fs.dfc.shutdown <= 1
            )

        if t == 2:
            # Slightly modify the constraint at t = 2, since startup variable is not
            # defined at t = 0
            return pb.period[t].fs.dfc.op_mode == pb.period[t - 1].fs.dfc.op_mode

        return (
            pb.period[t].fs.dfc.op_mode <= pb.period[t - 1].fs.dfc.op_mode +
            pb.period[t - 2].fs.dfc.startup
        )

    # If cycle startup is initiated at time t i.e., startup[t] = 1, then
    # at t + 1, op_mode, startup and shutdown must be 0
    # at t + 2, op_mode must be 1 and power output must be P_min
    # at t + 3, op_mode must be 1
    @m.Constraint(set_period)
    def startup_con_1(blk, t):
        if t + 1 > num_time_periods:
            return Constraint.Skip

        # Also ensures that at most one of the three operations (operation, startup, shutdown)
        # takes place at each time period
        return (
            pb.period[t + 1].fs.dfc.startup + pb.period[t + 1].fs.dfc.op_mode + 
            pb.period[t + 1].fs.dfc.shutdown <= 1 - pb.period[t].fs.dfc.startup
        )

    @m.Constraint(set_period)
    def startup_con_2(blk, t):
        if t + 2 > num_time_periods:
            return Constraint.Skip

        return pb.period[t].fs.dfc.startup <= pb.period[t + 2].fs.dfc.op_mode

    @m.Constraint(set_period)
    def startup_con_3(blk, t):
        if t + 3 > num_time_periods:
            return Constraint.Skip

        return pb.period[t].fs.dfc.startup <= pb.period[t + 3].fs.dfc.op_mode

    # Ensure that the power output is P_min at t + 2
    # old NL constraint linearized
    #  pb.period[t + 2].fs.dfc.power - DFC_OFFLOAD * des_mdl.dfc_design.capacity <=
    #  (1 - DFC_OFFLOAD) * des_mdl.dfc_design.capacity * (1 - pb.period[t].fs.dfc.startup)
    @m.Constraint(set_period)
    def startup_con_4(blk, t):
        if t + 2 > num_time_periods:
            return Constraint.Skip

        # Design models are contained in the parent block of pb
        des_mdl = pb.parent_block()

        return (
            pb.period[t + 2].fs.dfc.power <= des_mdl.dfc_design.capacity + (DFC_OFFLOAD - 1) 
            * des_mdl.dfc_design.capacity_ub * pb.period[t].fs.dfc.startup
        )
    
    @m.Constraint(set_period)
    def shutdown_con_4_1(blk, t):
        if t + 2 > num_time_periods:
          return Constraint.Skip

        des_mdl = pb.parent_block()

        return (
            pb.period[t + 2].fs.dfc.power <= (des_mdl.dfc_design.capacity * DFC_OFFLOAD 
            - des_mdl.dfc_design.capacity_lb + des_mdl.dfc_design.capacity_lb * 
            pb.period[t].fs.dfc.startup)/(DFC_OFFLOAD -1)
        )
    # Fuel requirement during the first hour of the startup
    # FIXME: The fuel consumption is for a fixed capacity of 838 MW. This needs to be modified 
    # if the design the DFC is a variable
    @m.Constraint(set_period)
    def startup_con_5(blk, t):
        # During the first hour, the fuel flowrate ramps up from 0 to 5.8518 kg/s.
        # In the model, we work with the average fuel flowrate: (0 + 5.8518) /2 = 2.9529
        return pb.period[t].fs.dfc.su_sd_ng_flow >= 2.9529 * pb.period[t].fs.dfc.startup

    @m.Constraint(set_period)
    def startup_con_6(blk, t):
        # During the second hour, the fuel flowrate ramps up from 5.8518 to 10.889 kg/s.
        # In the model, we work with the average fuel flowrate: (5.8518 + 10.889) / 2 = 8.3704
        # 10.889 kg/s is chosen instead of 10.36 to be consistent with the surrogate model
        if t + 1 > num_time_periods:
            return Constraint.Skip

        return pb.period[t + 1].fs.dfc.su_sd_ng_flow >= 8.3704 * pb.period[t].fs.dfc.startup

    # TODO: I suspect this constraint is implied by previous constraints, but I have not verified
    # it. So, I'm adding this constraint for now. If we can show that it is implied, we can remove
    # it later. Even if it is not implied, it is less likely that we will find a solution violating this
    # constraint. So, we can add it as a lazy constraint, if pyomo supports one.
    # Cycle must not be operating at time t - 1, if startup is initiated at time t
    @m.Constraint(set_period)
    def startup_con_7(blk, t):
        if t <= 1:
            return Constraint.Skip

        return pb.period[t].fs.dfc.startup <= 1 - pb.period[t - 1].fs.dfc.op_mode

    # If cycle shutdown is initiated at time t, i.e., shutdown[t] = 1, then
    # at time t - 2, op_mode must be 1 and power produced must be P_min
    # at time t - 1, op_mode must be 1 and power produced must be P_min
    # at time t, shutdown is initiated so power must be zero
    # at time t + 1, op_mode, startup and shutdown must be zero
    @m.Constraint(set_period)
    def shutdown_con_1(blk, t):
        if t <= 2:
            return Constraint.Skip

        return pb.period[t - 2].fs.dfc.op_mode >= pb.period[t].fs.dfc.shutdown

    @m.Constraint(set_period)
    def shutdown_con_2(blk, t):
        if t <= 1:
            return Constraint.Skip

        return pb.period[t - 1].fs.dfc.op_mode >= pb.period[t].fs.dfc.shutdown

    @m.Constraint(set_period)
    def shutdown_con_3(blk, t):
        if t + 1 > num_time_periods:
            return Constraint.Skip

        return (
            pb.period[t + 1].fs.dfc.startup + pb.period[t + 1].fs.dfc.op_mode + 
            pb.period[t + 1].fs.dfc.shutdown <= 1 - pb.period[t].fs.dfc.shutdown
        )

    # Ensure that the power output is P_min at t - 2
    # old cosntraint linearized
    # pb.period[t - 2].fs.dfc.power - DFC_OFFLOAD * des_mdl.dfc_design.capacity <=
    # (1 - DFC_OFFLOAD) * des_mdl.dfc_design.capacity * (1 - pb.period[t].fs.dfc.shutdown)
    @m.Constraint(set_period)
    def shutdown_con_4(blk, t):
        if t <= 2:
            return Constraint.Skip

        des_mdl = pb.parent_block()

        return (
            pb.period[t - 2].fs.dfc.power <= des_mdl.dfc_design.capacity + (DFC_OFFLOAD - 1) 
            * des_mdl.dfc_design.capacity_ub * pb.period[t].fs.dfc.shutdown
        )
    @m.Constraint(set_period)
    def shutdown_con_4_1(blk, t):
        if t <= 2:
          return Constraint.Skip

        des_mdl = pb.parent_block()

        return (
            pb.period[t - 2].fs.dfc.power <= (des_mdl.dfc_design.capacity * DFC_OFFLOAD 
            - des_mdl.dfc_design.capacity_lb + des_mdl.dfc_design.capacity_lb * 
            pb.period[t].fs.dfc.shutdown)/(DFC_OFFLOAD -1)
        )
    # Ensure that the power output is P_min at t - 1

    #            pb.period[t - 1].fs.dfc.power - DFC_OFFLOAD * des_mdl.dfc_design.capacity <=
    #            (1 - DFC_OFFLOAD) * des_mdl.dfc_design.capacity * (1 - pb.period[t].fs.dfc.shutdown)

    @m.Constraint(set_period)
    def shutdown_con_5(blk, t):
        if t <= 1:
            return Constraint.Skip

        des_mdl = pb.parent_block()

        return (
            pb.period[t - 1].fs.dfc.power <= des_mdl.dfc_design.capacity + (DFC_OFFLOAD - 1) 
            * des_mdl.dfc_design.capacity_ub * pb.period[t].fs.dfc.shutdown
        )
    

    @m.Constraint(set_period)
    def shutdown_con_5_1(blk, t):
        if t <= 2:
          return Constraint.Skip

        des_mdl = pb.parent_block()

        return (
            pb.period[t - 1].fs.dfc.power <= (des_mdl.dfc_design.capacity * DFC_OFFLOAD 
            - des_mdl.dfc_design.capacity_lb + des_mdl.dfc_design.capacity_lb * 
            pb.period[t].fs.dfc.shutdown)/(DFC_OFFLOAD -1)
        )
    # When shutdown is initiated, the fuel flowrate is ramped down from 10.889 to 5.8515
    # In the model, we use the average value: (10.889 + 5.8515) / 2 = 8.3704
    @m.Constraint(set_period)
    def shutdown_con_6(blk, t):
        return pb.period[t].fs.dfc.su_sd_ng_flow >= 8.3704 * pb.period[t].fs.dfc.shutdown

    # During the second hour, the fuel flowrate is ramped down from 5.8515 to 0
    # In the model, we use the average value: (5.8515 + 0) / 2 = 2.9529
    @m.Constraint(set_period)
    def shutdown_con_7(blk, t):
        if t + 1 > num_time_periods:
            return Constraint.Skip
        
        return pb.period[t + 1].fs.dfc.su_sd_ng_flow >= 2.9529 * pb.period[t].fs.dfc.shutdown 

    # If the plant was in operation at time t and *not* in operation at t + 1, then plant shutdown
    # must have been initiated at time t + 1
    @m.Constraint(set_period)
    def shutdown_con_8(blk, t):
        if t + 1 > num_time_periods:
            return Constraint.Skip

        return (
            pb.period[t].fs.dfc.op_mode - pb.period[t + 1].fs.dfc.op_mode <= 
            pb.period[t + 1].fs.dfc.shutdown
        ) 
