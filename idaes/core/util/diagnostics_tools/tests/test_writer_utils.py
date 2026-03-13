#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
This module contains tests for diagnostics writer utils.
"""

from io import StringIO

import pytest

from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    Objective,
    Var,
)
from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxBlock,
    ExternalGreyBoxModel,
)


from idaes.core.util.diagnostics_tools.writer_utils import (
    collect_model_statistics,
    write_report_section,
)

__author__ = "Alex Dowling, Douglas Allan, Andrew Lee"


@pytest.mark.unit
class TestStatsWriter:
    def test_blocks(self):
        m = ConcreteModel()
        m.b1 = Block()
        m.b1.b2 = Block()
        m.b3 = Block()
        m.b3.b4 = Block()
        m.b5 = Block()
        m.b5.b6 = Block()

        m.b3.b4.deactivate()
        m.b5.deactivate()

        stats = collect_model_statistics(m)

        assert len(stats) == 9
        assert stats[0] == "    Activated Blocks: 4 (Deactivated: 3)"
        assert (
            stats[1] == "    Free Variables in Activated Constraints: 0 (External: 0)"
        )
        assert stats[2] == "        Free Variables with only lower bounds: 0"
        assert stats[3] == "        Free Variables with only upper bounds: 0"
        assert stats[4] == "        Free Variables with upper and lower bounds: 0"
        assert (
            stats[5] == "    Fixed Variables in Activated Constraints: 0 (External: 0)"
        )
        assert stats[6] == "    Activated Equality Constraints: 0 (Deactivated: 0)"
        assert stats[7] == "    Activated Inequality Constraints: 0 (Deactivated: 0)"
        assert stats[8] == "    Activated Objectives: 0 (Deactivated: 0)"

    def test_constraints(self):
        m = ConcreteModel()
        m.b1 = Block()
        m.b2 = Block()

        m.b1.v1 = Var()
        m.b1.v2 = Var()
        m.b1.v3 = Var()
        m.b1.v4 = Var()

        m.b1.c1 = Constraint(expr=m.b1.v1 == m.b1.v2)
        m.b1.c2 = Constraint(expr=m.b1.v1 >= m.b1.v2)
        m.b1.c3 = Constraint(expr=m.b1.v1 == m.b1.v2)
        m.b1.c4 = Constraint(expr=m.b1.v1 >= m.b1.v2)

        m.b2.c1 = Constraint(expr=m.b1.v1 == m.b1.v2)
        m.b2.c2 = Constraint(expr=m.b1.v1 >= m.b1.v2)

        m.b2.deactivate()
        m.b1.c3.deactivate()
        m.b1.c4.deactivate()

        stats = collect_model_statistics(m)

        assert len(stats) == 9
        assert stats[0] == "    Activated Blocks: 2 (Deactivated: 1)"
        assert (
            stats[1] == "    Free Variables in Activated Constraints: 2 (External: 0)"
        )
        assert stats[2] == "        Free Variables with only lower bounds: 0"
        assert stats[3] == "        Free Variables with only upper bounds: 0"
        assert stats[4] == "        Free Variables with upper and lower bounds: 0"
        assert (
            stats[5] == "    Fixed Variables in Activated Constraints: 0 (External: 0)"
        )
        assert stats[6] == "    Activated Equality Constraints: 1 (Deactivated: 1)"
        assert stats[7] == "    Activated Inequality Constraints: 1 (Deactivated: 1)"
        assert stats[8] == "    Activated Objectives: 0 (Deactivated: 0)"

    def test_objectives(self):
        m = ConcreteModel()
        m.b1 = Block()
        m.b2 = Block()

        m.b1.o1 = Objective(expr=1)
        m.b1.o2 = Objective(expr=1)

        m.b2.o1 = Objective(expr=1)

        m.b2.deactivate()
        m.b1.o2.deactivate()

        stats = collect_model_statistics(m)

        assert len(stats) == 9
        assert stats[0] == "    Activated Blocks: 2 (Deactivated: 1)"
        assert (
            stats[1] == "    Free Variables in Activated Constraints: 0 (External: 0)"
        )
        assert stats[2] == "        Free Variables with only lower bounds: 0"
        assert stats[3] == "        Free Variables with only upper bounds: 0"
        assert stats[4] == "        Free Variables with upper and lower bounds: 0"
        assert (
            stats[5] == "    Fixed Variables in Activated Constraints: 0 (External: 0)"
        )
        assert stats[6] == "    Activated Equality Constraints: 0 (Deactivated: 0)"
        assert stats[7] == "    Activated Inequality Constraints: 0 (Deactivated: 0)"
        assert stats[8] == "    Activated Objectives: 1 (Deactivated: 1)"

    def test_free_variables(self):
        m = ConcreteModel()
        m.b1 = Block()

        m.v1 = Var(bounds=(0, 1))

        m.b1.v2 = Var(bounds=(0, None))
        m.b1.v3 = Var(bounds=(None, 0))
        m.b1.v4 = Var(bounds=(0, 1))
        m.b1.v5 = Var()
        m.b1.v6 = Var(bounds=(0, 1))

        m.b1.c1 = Constraint(expr=0 == m.v1 + m.b1.v2 + m.b1.v3 + m.b1.v4 + m.b1.v5)

        stats = collect_model_statistics(m.b1)

        assert len(stats) == 9
        assert stats[0] == "    Activated Blocks: 1 (Deactivated: 0)"
        assert (
            stats[1] == "    Free Variables in Activated Constraints: 5 (External: 1)"
        )
        assert stats[2] == "        Free Variables with only lower bounds: 1"
        assert stats[3] == "        Free Variables with only upper bounds: 1"
        assert stats[4] == "        Free Variables with upper and lower bounds: 2"
        assert (
            stats[5] == "    Fixed Variables in Activated Constraints: 0 (External: 0)"
        )
        assert stats[6] == "    Activated Equality Constraints: 1 (Deactivated: 0)"
        assert stats[7] == "    Activated Inequality Constraints: 0 (Deactivated: 0)"
        assert stats[8] == "    Activated Objectives: 0 (Deactivated: 0)"

    def test_fixed_variables(self):
        m = ConcreteModel()
        m.b1 = Block()

        m.v1 = Var(bounds=(0, 1))

        m.b1.v2 = Var(bounds=(0, None))
        m.b1.v3 = Var(bounds=(None, 0))
        m.b1.v4 = Var(bounds=(0, 1))
        m.b1.v5 = Var()
        m.b1.v6 = Var(bounds=(0, 1))

        m.b1.c1 = Constraint(expr=0 == m.v1 + m.b1.v2 + m.b1.v3 + m.b1.v4 + m.b1.v5)

        m.v1.fix(0.5)
        m.b1.v2.fix(-0.5)
        m.b1.v5.fix(-0.5)

        stats = collect_model_statistics(m.b1)

        assert len(stats) == 9
        assert stats[0] == "    Activated Blocks: 1 (Deactivated: 0)"
        assert (
            stats[1] == "    Free Variables in Activated Constraints: 2 (External: 0)"
        )
        assert stats[2] == "        Free Variables with only lower bounds: 0"
        assert stats[3] == "        Free Variables with only upper bounds: 1"
        assert stats[4] == "        Free Variables with upper and lower bounds: 1"
        assert (
            stats[5] == "    Fixed Variables in Activated Constraints: 3 (External: 1)"
        )
        assert stats[6] == "    Activated Equality Constraints: 1 (Deactivated: 0)"
        assert stats[7] == "    Activated Inequality Constraints: 0 (Deactivated: 0)"
        assert stats[8] == "    Activated Objectives: 0 (Deactivated: 0)"

    def test_with_greybox_variables(self):
        """non functional graybox model added to m fixture, to test DOFs

        GreyBoxModel has 3 inputs and 2 outputs calculated an unknown function,
        input a1 and a2 are bound by equality constraint through internal graybox model
        """

        class BasicGrayBox(ExternalGreyBoxModel):
            def input_names(self):
                return ["a1", "a2", "a3"]

            def output_names(self):
                return ["o1", "o2"]

            def equality_constraint_names(self):
                return ["a_sum"]

            def evaluate_equality_constraints(self):
                a1 = self._input_values[0]
                a2 = self._input_values[1]
                return [a1 * 0.5 + a2]

        m = ConcreteModel()

        m.gb = ExternalGreyBoxBlock(external_model=BasicGrayBox())
        m.gb_inactive = ExternalGreyBoxBlock(external_model=BasicGrayBox())
        m.gb_inactive.deactivate()
        m.a1 = Var(initialize=1)
        m.a1.fix()
        m.gb.inputs["a2"].unfix()
        m.gb.inputs["a3"].fix()
        m.a1_eq = Constraint(expr=m.a1 == m.gb.inputs["a1"])
        m.o1 = Var(initialize=1)
        m.o1_eq = Constraint(expr=m.o1 == m.gb.outputs["o1"])
        m.o1.fix()
        stats = collect_model_statistics(m)
        for k in stats:
            print(k)
        print(stats, len(stats))
        m.display()
        tab = " " * 4
        assert len(stats) == 13
        assert stats[0] == f"{tab}Activated Blocks: 1 (Deactivated: 0)"
        assert (
            stats[1] == f"{tab}Free Variables in Activated Constraints: 4 (External: 0)"
        )
        assert stats[2] == f"{tab*2}Free Variables with only lower bounds: 0"
        assert stats[3] == f"{tab*2}Free Variables with only upper bounds: 0"
        assert stats[4] == f"{tab*2}Free Variables with upper and lower bounds: 0"
        assert (
            stats[5]
            == f"{tab}Fixed Variables in Activated Constraints: 3 (External: 0)"
        )
        assert stats[6] == f"{tab}Activated Equality Constraints: 5 (Deactivated: 3)"
        assert stats[7] == f"{tab}Activated Inequality Constraints: 0 (Deactivated: 0)"
        assert stats[8] == f"{tab}Activated Objectives: 0 (Deactivated: 0)"
        assert stats[9] == f"{tab}GreyBox Statistics"
        assert stats[10] == f"{tab*2}Activated GreyBox models: 1 (Deactivated: 1)"
        assert stats[11] == f"{tab*2}Activated GreyBox Equalities: 3 (Deactivated: 3)"
        assert (
            stats[12]
            == f"{tab*2}Free Variables in Activated GreyBox Equalities: 4 (Fixed: 1)"
        )


@pytest.mark.unit
def test_write_report_section_all():
    stream = StringIO()

    write_report_section(
        stream=stream,
        lines_list=["a", "b", "c"],
        title="foo",
        line_if_empty="bar",
        end_line="baz",
        header="-",
        footer="=",
    )

    expected = """------------------------------------------------------------------------------------
foo

    a
    b
    c

baz
====================================================================================
"""
    assert stream.getvalue() == expected


@pytest.mark.unit
def test_write_report_section_no_lines():
    stream = StringIO()

    write_report_section(
        stream=stream,
        lines_list=[],
        title="foo",
        line_if_empty="bar",
        end_line="baz",
        header="-",
        footer="=",
    )

    expected = """------------------------------------------------------------------------------------
foo

    bar

baz
====================================================================================
"""
    assert stream.getvalue() == expected


@pytest.mark.unit
def test_write_report_section_lines_only():
    stream = StringIO()

    write_report_section(
        stream=stream,
        lines_list=["a", "b", "c"],
    )

    expected = """------------------------------------------------------------------------------------
    a
    b
    c

"""
    assert stream.getvalue() == expected
