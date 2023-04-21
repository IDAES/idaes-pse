import pyomo.environ as pyo


def surface_tension_type01(model, parameters):
    s = parameters["transport"]["surface_tension"]["s"]
    n = parameters["transport"]["surface_tension"]["n"]
    tc = parameters["transport"]["surface_tension"]["Tc"]
    return sum(s[i] * (1 - model.T_star / model.tau / tc) ** n[i] for i in s)
