from .phi_ideal_type01 import phi_ideal_expressions_type01
from .phi_ideal_type02 import phi_ideal_expressions_type02
from .phi_ideal_type03 import phi_ideal_expressions_type03
from .phi_residual_type01 import phi_residual_expressions_type01
from .phi_residual_type02 import phi_residual_expressions_type02
from .phi_residual_type03 import phi_residual_expressions_type03
from .phi_residual_type04 import phi_residual_expressions_type04
from .surface_tension_type01 import surface_tension_type01
from .sat_delta_approx import sat_delta_type01, sat_delta_type02

phi_residual_types = {
    0: None,  # custom
    1: phi_residual_expressions_type01,
    2: phi_residual_expressions_type02,
    3: phi_residual_expressions_type03,
    4: phi_residual_expressions_type04,
}

phi_ideal_types = {
    0: None,  # custom
    1: phi_ideal_expressions_type01,
    2: phi_ideal_expressions_type02,
    3: phi_ideal_expressions_type03,
}

delta_sat_types = {
    0: None,  # custom
    1: sat_delta_type01,
    2: sat_delta_type02,
}

surface_tension_types = {
    0: None,
    1: surface_tension_type01,
}
