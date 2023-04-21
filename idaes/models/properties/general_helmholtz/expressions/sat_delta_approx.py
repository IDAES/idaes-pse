import pyomo.environ as pyo


def sat_delta_type01(model, name, parameters):
    c = parameters["aux"][name]["c"]
    n = parameters["aux"][name]["n"]
    t = parameters["aux"][name]["t"]
    return c + sum(n[i] * (1 - 1 / model.tau) ** t[i] for i in n)


def sat_delta_type02(model, name, parameters):
    c = parameters["aux"][name]["c"]
    n = parameters["aux"][name]["n"]
    t = parameters["aux"][name]["t"]
    return c * pyo.exp(sum(n[i] * (1 - 1 / model.tau) ** t[i] for i in n))
