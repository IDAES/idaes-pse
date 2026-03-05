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
Pydantic schema for the model variables returned by the ModelVariables action
in `runner_actions.py`.

The structure is:

```
{"<IDAES/pyomo.component.name>":
    {
        "t": "<IDAES/Pyomo type>",
        "sub": {
           "<nested.IDAES/pyomo.component.name">:
           {
                # ..etc..
                # ..until you reach a "leaf"
                #  parameter or variable, which you
                # can tell with special values for "t"
                # t = "V" => variable
                  "t": "V",
                  # List of one or more values. For a scalar variable,
                  # there is only one value and the 'index' is None/null
                  v:  [[index, value, fixed, stale, lower_bound, upper_bound], ...]
                # t = "P" => parameter
                  "t": "P",
                  # List of one or more values. For a scalar parameter,
                  # there is only one value and the 'index' is None/null
                  v: [[index, value], ...]

           }
        }
    }
}
```

"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, RootModel, field_validator

# JSON-compatible scalar values.
JsonScalar = str | int | float | bool | None
# Index values are scalar for scalar vars/params, or a flat tuple/list for indexed ones.
IndexValue = JsonScalar | list[JsonScalar] | tuple[JsonScalar, ...]


class ParameterLeaf(BaseModel):
    """Leaf node for parameter values."""

    t: Literal["P"]
    # [index, value]
    v: list[tuple[IndexValue, JsonScalar]]


class VariableLeaf(BaseModel):
    """Leaf node for variable values."""

    t: Literal["V"]
    # [index, value, fixed, stale, lower_bound, upper_bound]
    v: list[tuple[IndexValue, JsonScalar, bool, bool, JsonScalar, JsonScalar]]


class BlockNode(BaseModel):
    """Tree node for blocks/components with nested children."""

    t: str
    sub: dict[str, "ModelVarNode"]

    @field_validator("t")
    @classmethod
    def _block_type_must_not_use_leaf_codes(cls, value: str) -> str:
        if value in {"P", "V"}:
            raise ValueError("Block node type must be a class name, not leaf code")
        return value


ModelVarNode = BlockNode | ParameterLeaf | VariableLeaf


class ModelVarsSchema(RootModel[dict[str, ModelVarNode]]):
    """Top-level schema for model variables JSON."""


BlockNode.model_rebuild()
