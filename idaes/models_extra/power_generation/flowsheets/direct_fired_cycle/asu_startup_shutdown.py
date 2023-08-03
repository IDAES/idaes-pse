"""
Contains functions that implement startup and shutdown constraints for the ASU

In this version, we implement the following and startup and shutdown procedure.

Cycle startup:
t = 0 - 8 => Consumes 80% of the max_power, but does not produce any oxygen
t > 8 => Can operate anywhere between P_min and P_max

Cycle shutdown:
t = 0 - 1 => Consumes 50% of the max_power and does not produce any oxygen
"""

from pyomo.environ import Constraint

# Startup time for the ASU
SU_TIME = 8


def asu_startup_shutdown_constraints(m):
    # The startup and shutdown constraints will be added in a separate block 
    # for eCach of the processes.

    # Create an alias for the parent block. 
    pb = m.parent_block()
    set_period = pb.set_period
    num_time_periods = len(set_period)

    # Cycle can operate at time t, only if one of the following two is true
    #     - Cycle was in operation at time t - 1
    #     - Cycle startup was initiated at time t - 8
    @m.Constraint(set_period)
    def cycle_operation(blk, t):
        if t == 1:
            # Constraint is not applicable at t = 1. Initial condition is left free
            # But it can choose at most one of the three
            return (
                pb.period[t].fs.asu.op_mode + pb.period[t].fs.asu.startup +
                pb.period[t].fs.asu.shutdown <= 1
            )

        # Slightly modify the constraint at t <= 8, since startup variable is not
        # defined for t <= 0
        return (
            pb.period[t].fs.asu.op_mode <= pb.period[t - 1].fs.asu.op_mode +
            (pb.period[t - SU_TIME].fs.asu.startup if t-SU_TIME > 0 else 0)
        )

    # After plant startup, the unit cannot produce oxygen in the next 8 hours, but
    # consumes 80% of P_max power during that period
    @m.Constraint(set_period, [j for j in range(1, SU_TIME)])
    def startup_con_1(blk, t, i):
        if t + i > num_time_periods:
            return Constraint.Skip

        return (
            pb.period[t + i].fs.asu.startup + pb.period[t + i].fs.asu.op_mode +
            pb.period[t + i].fs.asu.shutdown <= 1 - pb.period[t].fs.asu.startup
        )



    @m.Constraint(set_period, [j for j in range(SU_TIME)])
    def startup_con_2(blk, t, i):
        if t + i > num_time_periods:
            return Constraint.Skip

        des_mdl = pb.parent_block()
        return (
            pb.period[t + i].fs.asu.su_sd_power >= 0.8 * des_mdl.asu_design.max_power
            - 0.8 * des_mdl.asu_design.max_power_ub * (1-pb.period[t].fs.asu.startup)
        )

    # ASU must operate at the ninth hour after startup is initiated
    @m.Constraint(set_period)
    def startup_con_3(blk, t):
        if t + SU_TIME > num_time_periods:
            return Constraint.Skip

        return pb.period[t + SU_TIME].fs.asu.op_mode >= pb.period[t].fs.asu.startup

    # Startup can be initiated at time t only if op_mode[t - 1] = 0
    @m.Constraint(set_period)
    def startup_con_4(blk, t):
        if t == 1:
            return Constraint.Skip

        return pb.period[t].fs.asu.startup <= 1 - pb.period[t - 1].fs.asu.op_mode

    # Shutdown can be initiated at time t only if the ASU was operating at time t - 1
    @m.Constraint(set_period)
    def shutdown_con_1(blk, t):
        if t == 1:
            return Constraint.Skip

        # I suspect this will always be satisfied by the optimal solution. Need to verify it though
        return pb.period[t].fs.asu.shutdown <= pb.period[t - 1].fs.asu.op_mode

    # During the shutdown, power consumption is 50% of the P_max and O2 is not produced
    @m.Constraint(set_period)
    def shutdown_con_2(blk, t):
        des_mdl = pb.parent_block()
        return (
            pb.period[t].fs.asu.su_sd_power >= 0.5 * des_mdl.asu_design.max_power 
            - 0.5 *  des_mdl.asu_design.max_power_ub * (1-pb.period[t].fs.asu.shutdown)
        )

    # If op_mode[t] = 1, and op_mode[t + 1] = 0, then shutdown must be initated at t + 1
    @m.Constraint(set_period)
    def shutdown_con_3(blk, t):
        if t + 1 > num_time_periods:
            return Constraint.Skip

        return (
            pb.period[t].fs.asu.op_mode - pb.period[t + 1].fs.asu.op_mode <= pb.period[t + 1].fs.asu.shutdown
        )
