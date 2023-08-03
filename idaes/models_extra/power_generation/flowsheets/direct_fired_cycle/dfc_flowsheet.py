from pyomo.environ import (
    Var,
    NonNegativeReals,
    Constraint,
    Expression,
)
from idaes.core import FlowsheetBlock
from idaes.models_extra.power_generation.flowsheets.direct_fired_cycle.unit_models import (
    DFCOperation,
    MonoASUOperation,
    NLUOperation,
    OxygenTankOperation,
)
from idaes.models_extra.power_generation.flowsheets.direct_fired_cycle.default_model_parameters import NG_HHV, HR_TO_SEC


def build_dfc_flowsheet(
    m,
    dfc_design,
    asu_design,
):
    """
    Build a flowsheet with DFC and ASU units only. ASU provides all the oxygen required
    for the power cycle.
    """

    m.fs = FlowsheetBlock(dynamic=False)

    # Append unit models
    m.fs.dfc = DFCOperation(design_blk=dfc_design)
    m.fs.asu = MonoASUOperation(design_blk=asu_design)

    # Declare new variables
    m.fs.power_dfc_to_grid = Var(within=NonNegativeReals)
    m.fs.power_dfc_to_asu = Var(within=NonNegativeReals)
    m.fs.power_grid_to_asu = Var(within=NonNegativeReals)

    m.fs.oxygen_asu_to_dfc = Var(within=NonNegativeReals)
    m.fs.oxygen_asu_to_vent = Var(within=NonNegativeReals)

    # Power balance across the power cycle
    m.fs.dfc_power_balance = Constraint(
        expr=m.fs.dfc.power == m.fs.power_dfc_to_grid + m.fs.power_dfc_to_asu,
    )

    # Power balance across the ASU
    m.fs.asu_power_balance = Constraint(
        expr=m.fs.asu.total_power == m.fs.power_dfc_to_asu + m.fs.power_grid_to_asu,
    )

    # Oxygen balance across the power cycle
    m.fs.dfc_oxygen_balance = Constraint(
        expr=m.fs.dfc.o2_flow == m.fs.oxygen_asu_to_dfc,
    )

    # Oxygen balance across the ASU
    m.fs.asu_oxygen_balance = Constraint(
        expr=m.fs.asu.o2_flow == m.fs.oxygen_asu_to_dfc + m.fs.oxygen_asu_to_vent,
    )


def build_dfc_flowsheet_with_lox(
    m,
    dfc_design,
    asu_design,
    tank_design,
):
    """Build a flowsheet with DFC, ASU unit with a provision for LOx withdrawal,
    and an oxygen storage tank.
    """
    m.fs = FlowsheetBlock(dynamic=False)
    
    # Append unit models
    m.fs.dfc = DFCOperation(design_blk=dfc_design)
    m.fs.asu = MonoASUOperation(design_blk=asu_design)
    m.fs.tank = OxygenTankOperation(design_blk=tank_design)

    # Declare new variables
    m.fs.power_dfc_to_grid = Var(within=NonNegativeReals)
    m.fs.power_dfc_to_asu = Var(within=NonNegativeReals)
    m.fs.power_dfc_to_tank = Var(within=NonNegativeReals)
    m.fs.power_grid_to_asu = Var(within=NonNegativeReals)
    m.fs.power_grid_to_tank = Var(within=NonNegativeReals)

    m.fs.gox_fraction = Var(
        within=NonNegativeReals, 
        bounds=(0.9, 1),
        doc="Fraction of O2 withdrawn in gasoues state from the ASU",
    )
    m.fs.oxygen_asu_to_dfc = Var(within=NonNegativeReals)
    m.fs.oxygen_asu_to_vent = Var(within=NonNegativeReals)

    # Power balance across the power cycle
    m.fs.dfc_power_balance = Constraint(
        expr=m.fs.dfc.power == m.fs.power_dfc_to_grid + m.fs.power_dfc_to_asu + 
        m.fs.power_dfc_to_tank,
    )

    # Power balance across the ASU
    m.fs.asu_power_balance = Constraint(
        expr=2 * m.fs.asu.total_power - m.fs.gox_fraction * m.fs.asu.total_power == 
        m.fs.power_dfc_to_asu + m.fs.power_grid_to_asu,
    )

    # Power balance across the tank
    m.fs.tank_power_balance = Constraint(
        expr=m.fs.tank.power == m.fs.power_dfc_to_tank + m.fs.power_grid_to_tank,
    )

    # Oxygen balance across the power cycle
    m.fs.dfc_oxygen_balance = Constraint(
        expr=m.fs.dfc.o2_flow == m.fs.oxygen_asu_to_dfc + m.fs.tank.lox_out,
    )

    # Gasoues oxygen balance across the ASU
    m.fs.gox_balance = Constraint(
        expr=m.fs.asu.o2_flow * m.fs.gox_fraction == m.fs.oxygen_asu_to_dfc +
        m.fs.oxygen_asu_to_vent,
    )

    # Liquid oxygen balance across the ASU
    m.fs.lox_balance = Constraint(
        expr=m.fs.asu.o2_flow - m.fs.asu.o2_flow * m.fs.gox_fraction ==
        m.fs.tank.lox_in,
    )


def build_dfc_flowsheet_with_nlu(
    m,
    dfc_design,
    asu_design,
    nlu_design,
    tank_design
):
    """
    Build a flowsheet with DFC, ASU, NLU and an oxygen storage tank.
    """
    m.fs = FlowsheetBlock(dynamic=False)

    # Append unit models
    m.fs.dfc = DFCOperation(design_blk=dfc_design)
    m.fs.asu = MonoASUOperation(design_blk=asu_design)
    m.fs.nlu = NLUOperation(design_blk=nlu_design)
    m.fs.tank = OxygenTankOperation(design_blk=tank_design)

    # Declare new variables
    m.fs.power_dfc_to_grid = Var(within=NonNegativeReals)
    m.fs.power_dfc_to_asu = Var(within=NonNegativeReals)
    m.fs.power_dfc_to_nlu = Var(within=NonNegativeReals)
    m.fs.power_dfc_to_tank = Var(within=NonNegativeReals)
    m.fs.power_grid_to_asu = Var(within=NonNegativeReals)
    m.fs.power_grid_to_nlu = Var(within=NonNegativeReals)
    m.fs.power_grid_to_tank = Var(within=NonNegativeReals)

    m.fs.oxygen_asu_to_dfc = Var(within=NonNegativeReals)
    m.fs.oxygen_asu_to_vent = Var(within=NonNegativeReals)

    # Power balance across the power cycle
    m.fs.dfc_power_balance = Constraint(
        expr=m.fs.dfc.power == m.fs.power_dfc_to_grid + m.fs.power_dfc_to_asu + 
        m.fs.power_dfc_to_nlu + m.fs.power_dfc_to_tank,
    )

    # Power balance across the ASU
    m.fs.asu_power_balance = Constraint(
        expr=m.fs.asu.total_power == m.fs.power_dfc_to_asu + m.fs.power_grid_to_asu
    )

    # Power balance across the NLU
    m.fs.nlu_power_balance = Constraint(
        expr=m.fs.nlu.power == m.fs.power_dfc_to_nlu + m.fs.power_grid_to_nlu
    )

    # Power balance across the tank
    m.fs.tank_power_balance = Constraint(
        expr=m.fs.tank.power == m.fs.power_dfc_to_tank + m.fs.power_grid_to_tank
    )

    # Oxygen balance across the power cycle
    m.fs.dfc_oxygen_balance = Constraint(
        expr=m.fs.dfc.o2_flow == m.fs.oxygen_asu_to_dfc + m.fs.tank.lox_out,
    )

    # Oxygen balance across the ASU
    m.fs.asu_oxygen_balance = Constraint(
        expr=m.fs.asu.o2_flow == m.fs.oxygen_asu_to_dfc + m.fs.nlu.o2_flow +
        m.fs.oxygen_asu_to_vent,
    )

    # Oxygen balance across the NLU
    m.fs.nlu_oxygen_balance = Constraint(
        expr=m.fs.nlu.o2_flow == m.fs.tank.lox_in,
    ) 


def append_op_costs_dfc(
    m,
    lmp,
    penalty,
    cost_ng,
    carbon_price,
):
    """Append cost and revenue expressions for each time step

    Args:
        m: Object containing the flowsheet
        lmp: LMP 
        penalty: To avoid degeneracy, we impose penalty for borrowing electricity
                 from the grid. (default: electricity borrowed from the grid costs $5/MWh
                 more than the LMP at that time instant)
        cost_ng: Cost of natural gas [in $/MMBtu]
    """

    m.fs.electricity_revenue = Expression(expr=lmp * m.fs.power_dfc_to_grid) 

    m.fs.electricity_cost = Expression(    
        expr=(lmp + penalty) * (m.fs.power_grid_to_asu
        + (m.fs.power_grid_to_nlu if hasattr(m.fs, "nlu") else 0)
        + (m.fs.power_grid_to_tank if hasattr(m.fs, "tank") else 0)
        ),
        doc="Net revenue generated by selling power to the grid [in $1000]",
    )

    m.fs.fuel_cost = Expression(
        expr=cost_ng * NG_HHV * (HR_TO_SEC / 1000) * m.fs.dfc.total_ng_flow,
        doc="Operational cost associated with fuel consumption [in $1000]",
    )

    # TODO: Add CO2 emissions for the startup and shutdown (likely negligible)
    m.fs.co2_price = Expression(
        expr=m.fs.dfc.co2_emissions * carbon_price / 1000,
        doc="Carbon price per hour [in $1000]",
    )

    m.fs.net_cash_flow = Expression(
        expr=m.fs.electricity_revenue - m.fs.electricity_cost
        - m.fs.fuel_cost
        - m.fs.dfc.non_fuel_vom - m.fs.asu.non_fuel_vom 
        - (m.fs.nlu.non_fuel_vom if hasattr(m.fs, "nlu") else 0)
        - m.fs.co2_price,
    )


def append_cashflows(m, params):
    plant_life = params["plant_life"]
    tax_rate = params["tax_rate"]
    discount_rate = params["discount_rate"]

    m.CAPEX = Var(
        within=NonNegativeReals,
        doc="Total CAPEX of the plant [in $1000]",
    )
    m.FOM = Var(
        within=NonNegativeReals,
        doc="Total fixed O&M costs per year [in $1000]",
    )
    m.DEPRECIATION = Var(
        within=NonNegativeReals,
        doc="Depreciation value per year [in $1000]"
    )
    m.CORP_TAX = Var(
        within=NonNegativeReals,
        doc="Net corporate tax per year [in $1000]",
    )
    m.NET_PROFIT = Var(
        doc="Net profit per year [in $1000]",
    )

    m.capex_calculation = Constraint(
        expr=m.CAPEX == m.dfc_design.capex + m.asu_design.capex + 
        (m.nlu_design.capex if hasattr(m, "nlu_design") else 0) +
        (m.tank_design.capex if hasattr(m, "tank_design") else 0),
        doc="Calculates the total CAPEX [in $1000]",
    )
    m.fom_calculation = Constraint(
        expr=m.FOM == m.dfc_design.fom + m.asu_design.fom +
        (m.nlu_design.fom if hasattr(m, "nlu_design") else 0) +
        (m.tank_design.fom if hasattr(m, "tank_design") else 0),
        doc="Calculates the total FOM [in $1000]",
    )
    m.depreciation_calculation = Constraint(
        expr=m.DEPRECIATION == m.CAPEX / plant_life,
        doc="Straight line depreciation with zero salvage value [in $1000]",
    )
    
    # Tax = max{0, Formula below}. We will relax it to Tax >= max{0, Formula below}
    # The inequality will be binding for the optimal solution. 
    m.corp_tax_calculation = Constraint(
        expr=m.CORP_TAX >= tax_rate * (
            sum(m.mp_model.period[t].fs.net_cash_flow for t in m.set_time)
            - m.FOM - m.DEPRECIATION
        ),
        doc="Calculates the total corporate tax [in $1000]",
    )

    m.net_profit_calculation = Constraint(
        expr=m.NET_PROFIT == 
        + sum(m.mp_model.period[t].fs.net_cash_flow for t in m.set_time)
        - m.FOM - m.CORP_TAX,
        doc="Calculate the net profit generated in a year [in $1000]",
    )

    constant_cf_factor = (1 - (1 + discount_rate) ** (- plant_life)) / discount_rate

    # NPV Calculation
    m.npv = Expression(
        expr=constant_cf_factor * m.NET_PROFIT - m.CAPEX,
        doc="Calculates the net present value [in $1000]",
    )
