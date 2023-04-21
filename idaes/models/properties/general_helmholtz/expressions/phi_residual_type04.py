import pyomo.environ as pyo


def phi_residual_expressions_type04(model, parameters):
    last_term = parameters["eos"]["last_term_residual"]
    n = parameters["eos"]["n"]
    t = parameters["eos"]["t"]
    d = parameters["eos"]["d"]
    first_term = 1
    rng = []
    for last_term in last_term:
        rng.append(range(first_term, last_term + 1))
        first_term = last_term + 1
    return {
        "phir": sum(n[i] * model.delta ** d[i] * model.tau ** t[i] for i in rng[0])
        + sum(
            pyo.exp(-model.delta**k)
            * sum(n[i] * model.delta ** d[i] * model.tau ** t[i] for i in rng[k])
            for k in range(1, len(rng))
        ),
        "phir_d": sum(
            n[i] * d[i] * model.delta ** (d[i] - 1) * model.tau ** t[i] for i in rng[0]
        )
        + sum(
            pyo.exp(-model.delta**k)
            * sum(
                n[i]
                * (d[i] - k * model.delta**k)
                * model.delta ** (d[i] - 1)
                * model.tau ** t[i]
                for i in rng[k]
            )
            for k in range(1, len(rng))
        ),
        "phir_dd": sum(
            n[i] * d[i] * (d[i] - 1) * model.delta ** (d[i] - 2) * model.tau ** t[i]
            for i in rng[0]
        )
        + sum(
            pyo.exp(-model.delta**k)
            * sum(
                n[i]
                * (
                    -k * model.delta**k * (d[i] - k * model.delta**k)
                    + (d[i] - 1) * (d[i] - k * model.delta**k)
                    - k**2 * model.delta**k
                )
                * model.delta ** (d[i] - 2)
                * model.tau ** t[i]
                for i in rng[k]
            )
            for k in range(1, len(rng))
        ),
        "phir_t": sum(
            n[i] * t[i] * model.delta ** d[i] * model.tau ** (t[i] - 1) for i in rng[0]
        )
        + sum(
            pyo.exp(-model.delta**k)
            * sum(
                n[i] * t[i] * model.delta ** d[i] * model.tau ** (t[i] - 1)
                for i in rng[k]
            )
            for k in range(1, len(rng))
        ),
        "phir_tt": sum(
            n[i] * t[i] * (t[i] - 1) * model.delta ** d[i] * model.tau ** (t[i] - 2)
            for i in rng[0]
        )
        + sum(
            pyo.exp(-model.delta**k)
            * sum(
                n[i] * t[i] * (t[i] - 1) * model.delta ** d[i] * model.tau ** (t[i] - 2)
                for i in rng[k]
            )
            for k in range(1, len(rng))
        ),
        "phir_dt": sum(
            n[i] * d[i] * t[i] * model.delta ** (d[i] - 1) * model.tau ** (t[i] - 1)
            for i in rng[0]
        )
        + sum(
            pyo.exp(-model.delta**k)
            * sum(
                n[i]
                * t[i]
                * (d[i] - k * model.delta**k)
                * model.delta ** (d[i] - 1)
                * model.tau ** (t[i] - 1)
                for i in rng[k]
            )
            for k in range(1, len(rng))
        ),
    }
