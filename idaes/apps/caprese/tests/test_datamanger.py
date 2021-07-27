from collections import OrderedDict

import pyomo.environ as pyo
from idaes.apps.caprese.data_manager import (
        empty_dataframe_from_variables,
        add_variable_values_to_dataframe,
        )

from pyomo.core.base.componentuid import ComponentUID
from pyomo.common.collections import ComponentMap
import pandas as pd


def test_dataframes():
    m = pyo.ConcreteModel()
    m.var = pyo.Var(['a','b'], [0,1,2,3,4],
                    initialize={
                        ("a",0): 1.1,
                        ("a",1): 1.2,
                        ("a",2): 1.3,
                        ("a",3): 1.4,
                        ("a",4): 1.5,
                        ("b",0): 2.1,
                        ("b",1): 2.2,
                        ("b",2): 2.3,
                        ("b",3): 2.4,
                        ("b",4): 2.5,
                    })
    variables = [
        pyo.Reference(m.var["a",:]),
        pyo.Reference(m.var["b",:]),
    ]
    df = empty_dataframe_from_variables(variables)
    df = add_variable_values_to_dataframe(
        df,
        variables,
        0,
        time_subset=[0,1,2],
    )
    df = add_variable_values_to_dataframe(
        df,
        variables,
        1,
        time_subset=[3,4],
    )
    print(df)
    
    return df, m


def main():
    df, m = test_dataframes()
    
    return df, m

if __name__ == "__main__":
    df, m = main()
    
