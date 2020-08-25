##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
Waterwall section model test

main equations:

* Heat is given by fire-side boiler model
* Calculate pressure change due to friction and gravity
* Calculate slag layer wall temperature
* Consider a layer of metal and a layer of slag


Created on Thu Aug 24 2020 by Boiler Team (J. Ma, M. Zamarripa)
"""
import pytest
# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.network import Arc

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom


# Import Unit Model Modules
from idaes.generic_models.properties import iapws95

# Import IDAES standard unit model
import idaes.logger as idaeslog

from idaes.power_generation.unit_models.waterwall_section import \
    WaterwallSection
from idaes.core.util.testing import (get_default_solver)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()

# -----------------------------------------------------------------------------


@pytest.fixture
def build_waterwall():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop_water = iapws95.Iapws95ParameterBlock()
    m.fs.ww_zones = pyo.RangeSet(10)
    m.fs.Waterwalls = WaterwallSection(m.fs.ww_zones,
                                       default={
                                           "dynamic": False,
                                           "has_holdup": False,
                                           "property_package": m.fs.prop_water,
                                           "has_heat_transfer": True,
                                           "has_pressure_change": True,
                                           },
                                       )


    # build connections
    m.fs.stream_ww1 = Arc(
        source=m.fs.Waterwalls[1].outlet, destination=m.fs.Waterwalls[2].inlet
    )
    m.fs.stream_ww2 = Arc(
        source=m.fs.Waterwalls[2].outlet, destination=m.fs.Waterwalls[3].inlet
    )
    m.fs.stream_ww3 = Arc(
        source=m.fs.Waterwalls[3].outlet, destination=m.fs.Waterwalls[4].inlet
    )
    m.fs.stream_ww4 = Arc(
        source=m.fs.Waterwalls[4].outlet, destination=m.fs.Waterwalls[5].inlet
    )
    m.fs.stream_ww5 = Arc(
        source=m.fs.Waterwalls[5].outlet, destination=m.fs.Waterwalls[6].inlet
    )
    m.fs.stream_ww6 = Arc(
        source=m.fs.Waterwalls[6].outlet, destination=m.fs.Waterwalls[7].inlet
    )
    m.fs.stream_ww7 = Arc(
        source=m.fs.Waterwalls[7].outlet, destination=m.fs.Waterwalls[8].inlet
    )
    m.fs.stream_ww8 = Arc(
        source=m.fs.Waterwalls[8].outlet, destination=m.fs.Waterwalls[9].inlet
    )
    m.fs.stream_ww9 = Arc(
        source=m.fs.Waterwalls[9].outlet,
        destination=m.fs.Waterwalls[10].inlet
    )

    # Pyomo expands arcs writing constraints outlet unit 1 = inlet unit 2
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)

    return m


@pytest.mark.unit
def test_basic_build(build_waterwall):
    """Make a turbine model and make sure it doesn't throw exception"""
    m = build_waterwall


@pytest.mark.skipif(not iapws95.iapws95_available(),
                    reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_waterwall(build_waterwall):
    m = build_waterwall
    # fix inputs
    # 10 waterwall sections
    for i in m.fs.ww_zones:
        m.fs.Waterwalls[i].tube_di.fix(0.047)
        m.fs.Waterwalls[i].tube_thickness.fix(0.00350)
        m.fs.Waterwalls[i].fin_thickness.fix(0.00455)
        m.fs.Waterwalls[i].slag_thickness[:].fix(0.001)
        m.fs.Waterwalls[i].fin_length.fix(0.0115)
        m.fs.Waterwalls[i].number_tubes.fix(610)
        m.fs.Waterwalls[i].fcorrection_dp.fix(1.2)

    # water wall section height (must be equal to the fire side model zones)
    m.fs.Waterwalls[1].height.fix(6.150)
    m.fs.Waterwalls[2].height.fix(3.150)
    m.fs.Waterwalls[3].height.fix(1.5)
    m.fs.Waterwalls[4].height.fix(1.450)
    m.fs.Waterwalls[5].height.fix(1.350)
    m.fs.Waterwalls[6].height.fix(1.250)
    m.fs.Waterwalls[7].height.fix(1.150)
    m.fs.Waterwalls[8].height.fix(1.350)
    m.fs.Waterwalls[9].height.fix(3.250)
    m.fs.Waterwalls[10].height.fix(3.450)

    # water wall section projected area
    m.fs.Waterwalls[1].area_proj_total.fix(320.0)
    m.fs.Waterwalls[2].area_proj_total.fix(150.3)
    m.fs.Waterwalls[3].area_proj_total.fix(70.8)
    m.fs.Waterwalls[4].area_proj_total.fix(70.0)
    m.fs.Waterwalls[5].area_proj_total.fix(58.6)
    m.fs.Waterwalls[6].area_proj_total.fix(58.6)
    m.fs.Waterwalls[7].area_proj_total.fix(50.1)
    m.fs.Waterwalls[8].area_proj_total.fix(65.6)
    m.fs.Waterwalls[9].area_proj_total.fix(145.6)
    m.fs.Waterwalls[10].area_proj_total.fix(165.5)

    # Heat loss to waterwall Q in W
    m.fs.Waterwalls[1].heat_fireside[:].fix(2.3e7)
    m.fs.Waterwalls[2].heat_fireside[:].fix(1.5e7)
    m.fs.Waterwalls[3].heat_fireside[:].fix(6.9e6)
    m.fs.Waterwalls[4].heat_fireside[:].fix(1.2e7)
    m.fs.Waterwalls[5].heat_fireside[:].fix(1.2e7)
    m.fs.Waterwalls[6].heat_fireside[:].fix(1.2e7)
    m.fs.Waterwalls[7].heat_fireside[:].fix(1.1e7)
    m.fs.Waterwalls[8].heat_fireside[:].fix(9.9e6)
    m.fs.Waterwalls[9].heat_fireside[:].fix(2.2e7)
    m.fs.Waterwalls[10].heat_fireside[:].fix(1.9e7)

    outlvl = idaeslog.NOTSET
    init_log = idaeslog.getInitLogger(m.name, outlvl, tag="flowsheet")

    optarg = {"tol": 1e-7,
              "linear_solver": "ma27",
              "max_iter": 40}
    solver = pyo.SolverFactory("ipopt")
    solver.options = optarg
    init_log.info_low("Starting initialization...")

    # Set inlet and operating conditions, and some initial conditions.
    m.fs.Waterwalls[1].inlet.flow_mol[0].fix(150055.0)  # mol/s
    m.fs.Waterwalls[1].inlet.enth_mol[0].fix(31000.0)   # J/mol
    m.fs.Waterwalls[1].inlet.pressure[0].fix(1.750e7)   # Pa

    m.fs.Waterwalls[1].initialize(
        state_args={
            "flow_mol": m.fs.Waterwalls[1].inlet.flow_mol[0].value,
            "pressure": m.fs.Waterwalls[1].inlet.pressure[0].value,
            "enth_mol": m.fs.Waterwalls[1].inlet.enth_mol[0].value,
        },
        outlvl=outlvl,
        optarg=solver.options,
    )

    for i in range(2, 11):
        m.fs.Waterwalls[i].inlet.flow_mol[:].fix(m.fs.Waterwalls[i - 1].
                                                 outlet.flow_mol[0].value)
        m.fs.Waterwalls[i].inlet.pressure[:].fix(m.fs.Waterwalls[i - 1].outlet.
                                                 pressure[0].value)
        m.fs.Waterwalls[i].inlet.enth_mol[:].fix(m.fs.Waterwalls[i - 1].outlet.
                                                 enth_mol[0].value)

        m.fs.Waterwalls[i].initialize(
            state_args={
                "flow_mol": m.fs.Waterwalls[i - 1].outlet.flow_mol[0].value,
                "pressure": m.fs.Waterwalls[i - 1].outlet.pressure[0].value,
                "enth_mol": m.fs.Waterwalls[i - 1].outlet.enth_mol[0].value,
                },
            outlvl=0,
            optarg=solver.options,
            )

    assert degrees_of_freedom(m) == -27  # all inputs fixed at this point



@pytest.mark.skipif(not iapws95.iapws95_available(),
                    reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_waterwall(build_waterwall):
    m = build_waterwall
    # fix inputs
    # 10 waterwall sections
    for i in m.fs.ww_zones:
        m.fs.Waterwalls[i].tube_di.fix(0.047)
        m.fs.Waterwalls[i].tube_thickness.fix(0.00350)
        m.fs.Waterwalls[i].fin_thickness.fix(0.00455)
        m.fs.Waterwalls[i].slag_thickness[:].fix(0.001)
        m.fs.Waterwalls[i].fin_length.fix(0.0115)
        m.fs.Waterwalls[i].number_tubes.fix(610)
        m.fs.Waterwalls[i].fcorrection_dp.fix(1.2)

    # water wall section height (must be equal to the fire side model zones)
    m.fs.Waterwalls[1].height.fix(6.150)
    m.fs.Waterwalls[2].height.fix(3.150)
    m.fs.Waterwalls[3].height.fix(1.5)
    m.fs.Waterwalls[4].height.fix(1.450)
    m.fs.Waterwalls[5].height.fix(1.350)
    m.fs.Waterwalls[6].height.fix(1.250)
    m.fs.Waterwalls[7].height.fix(1.150)
    m.fs.Waterwalls[8].height.fix(1.350)
    m.fs.Waterwalls[9].height.fix(3.250)
    m.fs.Waterwalls[10].height.fix(3.450)

    # water wall section projected area
    m.fs.Waterwalls[1].area_proj_total.fix(320.0)
    m.fs.Waterwalls[2].area_proj_total.fix(150.3)
    m.fs.Waterwalls[3].area_proj_total.fix(70.8)
    m.fs.Waterwalls[4].area_proj_total.fix(70.0)
    m.fs.Waterwalls[5].area_proj_total.fix(58.6)
    m.fs.Waterwalls[6].area_proj_total.fix(58.6)
    m.fs.Waterwalls[7].area_proj_total.fix(50.1)
    m.fs.Waterwalls[8].area_proj_total.fix(65.6)
    m.fs.Waterwalls[9].area_proj_total.fix(145.6)
    m.fs.Waterwalls[10].area_proj_total.fix(165.5)

    # Heat loss to waterwall Q in W
    m.fs.Waterwalls[1].heat_fireside[:].fix(2.3e7)
    m.fs.Waterwalls[2].heat_fireside[:].fix(1.5e7)
    m.fs.Waterwalls[3].heat_fireside[:].fix(6.9e6)
    m.fs.Waterwalls[4].heat_fireside[:].fix(1.2e7)
    m.fs.Waterwalls[5].heat_fireside[:].fix(1.2e7)
    m.fs.Waterwalls[6].heat_fireside[:].fix(1.2e7)
    m.fs.Waterwalls[7].heat_fireside[:].fix(1.1e7)
    m.fs.Waterwalls[8].heat_fireside[:].fix(9.9e6)
    m.fs.Waterwalls[9].heat_fireside[:].fix(2.2e7)
    m.fs.Waterwalls[10].heat_fireside[:].fix(1.9e7)

    outlvl = idaeslog.NOTSET
    init_log = idaeslog.getInitLogger(m.name, outlvl, tag="flowsheet")

    optarg = {"tol": 1e-7,
              "linear_solver": "ma27",
              "max_iter": 40}
    solver = pyo.SolverFactory("ipopt")
    solver.options = optarg
    init_log.info_low("Starting initialization...")

    # Set inlet and operating conditions, and some initial conditions.
    m.fs.Waterwalls[1].inlet.flow_mol[0].fix(150055.0)  # mol/s
    m.fs.Waterwalls[1].inlet.enth_mol[0].fix(31000.0)   # J/mol
    m.fs.Waterwalls[1].inlet.pressure[0].fix(1.750e7)   # Pa

    m.fs.Waterwalls[1].initialize(
        state_args={
            "flow_mol": m.fs.Waterwalls[1].inlet.flow_mol[0].value,
            "pressure": m.fs.Waterwalls[1].inlet.pressure[0].value,
            "enth_mol": m.fs.Waterwalls[1].inlet.enth_mol[0].value,
        },
        outlvl=outlvl,
        optarg=solver.options,
    )

    for i in range(2, 11):
        m.fs.Waterwalls[i].inlet.flow_mol[:].fix(m.fs.Waterwalls[i - 1].
                                                 outlet.flow_mol[0].value)
        m.fs.Waterwalls[i].inlet.pressure[:].fix(m.fs.Waterwalls[i - 1].outlet.
                                                 pressure[0].value)
        m.fs.Waterwalls[i].inlet.enth_mol[:].fix(m.fs.Waterwalls[i - 1].outlet.
                                                 enth_mol[0].value)
        m.fs.Waterwalls[i].initialize(
            state_args={
                "flow_mol": m.fs.Waterwalls[i - 1].outlet.flow_mol[0].value,
                "pressure": m.fs.Waterwalls[i - 1].outlet.pressure[0].value,
                "enth_mol": m.fs.Waterwalls[i - 1].outlet.enth_mol[0].value,
                },
            outlvl=0,
            optarg=solver.options,
            )

    #unfix inputs before solving the flowsheet
    for i in m.fs.ww_zones:
        if i == 1:
            m.fs.Waterwalls[i].inlet.flow_mol[:].fix()
            m.fs.Waterwalls[i].inlet.pressure[:].fix()
            m.fs.Waterwalls[i].inlet.enth_mol[:].fix()
        else:
            m.fs.Waterwalls[i].inlet.flow_mol[:].unfix()
            m.fs.Waterwalls[i].inlet.pressure[:].unfix()
            m.fs.Waterwalls[i].inlet.enth_mol[:].unfix()

    solver.solve(m, tee=True)

    assert (pytest.approx(0.08353, abs=1e-3) ==
            pyo.value(m.fs.Waterwalls[10].control_volume.
                      properties_out[0].vapor_frac))
    assert (pytest.approx(150055.0, abs=1e-3) ==
            pyo.value(m.fs.Waterwalls[10].control_volume.
                      properties_out[0].flow_mol))
    assert (pytest.approx(31951.65106, abs=1e-3) ==
            pyo.value(m.fs.Waterwalls[10].control_volume.
                      properties_out[0].enth_mol))
    assert degrees_of_freedom(m) == 0
