import pyomo.environ as pyo


def phi_residual_expressions_type01(model, parameters):
    last_term = parameters["eos"]["last_term_residual"]
    n = parameters["eos"]["n"]
    t = parameters["eos"]["t"]
    d = parameters["eos"]["d"]
    c = parameters["eos"]["c"]
    first_term = 1
    rng = []
    for last_term in last_term:
        rng.append(range(first_term, last_term + 1))
        first_term = last_term + 1
    return {
        "phir": sum(n[i] * model.delta ** d[i] * model.tau ** t[i] for i in rng[0])
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        ),
        "phir_d": sum(
            n[i] * d[i] * model.delta ** (d[i] - 1) * model.tau ** t[i] for i in rng[0]
        )
        + sum(
            n[i]
            * pyo.exp(-model.delta ** c[i])
            * model.delta ** (d[i] - 1)
            * model.tau ** t[i]
            * (d[i] - c[i] * model.delta ** c[i])
            for i in rng[1]
        ),
        "phir_dd": sum(
            n[i] * d[i] * (d[i] - 1) * model.delta ** (d[i] - 2) * model.tau ** t[i]
            for i in rng[0]
        )
        + sum(
            n[i]
            * pyo.exp(-model.delta ** c[i])
            * model.delta ** (d[i] - 2)
            * model.tau ** t[i]
            * (
                (d[i] - c[i] * model.delta ** c[i])
                * (d[i] - 1 - c[i] * model.delta ** c[i])
                - c[i] ** 2 * model.delta ** c[i]
            )
            for i in rng[1]
        ),
        "phir_t": sum(
            n[i] * t[i] * model.delta ** d[i] * model.tau ** (t[i] - 1) for i in rng[0]
        )
        + sum(
            n[i]
            * t[i]
            * model.delta ** d[i]
            * model.tau ** (t[i] - 1)
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        ),
        "phir_tt": sum(
            n[i] * t[i] * (t[i] - 1) * model.delta ** d[i] * model.tau ** (t[i] - 2)
            for i in rng[0]
        )
        + sum(
            n[i]
            * t[i]
            * (t[i] - 1)
            * model.delta ** d[i]
            * model.tau ** (t[i] - 2)
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        ),
        "phir_dt": sum(
            n[i] * t[i] * d[i] * model.delta ** (d[i] - 1) * model.tau ** (t[i] - 1)
            for i in rng[0]
        )
        + sum(
            n[i]
            * t[i]
            * model.delta ** (d[i] - 1)
            * model.tau ** (t[i] - 1)
            * (d[i] - c[i] * model.delta ** c[i])
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        ),
    }
