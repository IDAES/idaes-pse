from idaes.apps.caprese import NmpcManager
from pyomo.environ import SolverFactory

# plant = make_model(...)
# controller = make_model(...)
# sample_time = 1.

inputs = [
        plant.gas_inlet.flow[0],
        plant.solid_inlet.flow[0],
        ]

measurements = [
        *list(controller.gas_phase.material_holdup[0,...]),
        *list(controller.gas_phase.energy_holdup[0,...]),
        *list(controller.solid_phase.material_holdup[0,...]),
        *list(controller.solid_phase.energy_holdup[0,...]),
        ]

nmpc = NmpcManager(
        plant=(plant, plant.time, inputs),
        controller=(controller, controller.time, measurements),
        )

solver = SolverFactory('ipopt')

#setpoint = [
#        (vardata, value),
#        ]
nmpc.controller.add_setpoint_objective(setpoint)
with nmpc.controller.solve_setpoint_context():
    solver.solve(controller)

#tracking_weights = [
#        (vardata, value),
#        ]
nmpc.controller.add_tracking_objective(tracking_weights)

nmpc.controller.constrain_control_inputs_piecewise_constant(sample_time)

c_t0 = nmpc.controller.time.first()
p_t0 = nmpc.plant.time.first()

c_t1 = nmpc.controller.time[2]
p_t1 = nmpc.plant.time[2]

# Initialize from initial conditions:
for var in nmpc.controller.component_objects(NmpcVar):
    # var here is an NmpcVar; a custom time-indexed Var
    # that contains useful attributes such as setpoint
    var[:].set_value(var[c_t0])

# Unfix controller inputs:
nmpc.controller.vars.input[:,c_t1:].unfix()

# Solve control problem:
solver.solve(nmpc.controller, tee=True)

c_ts = nmpc.controller.sample_points[1]
p_ts = nmpc.controller.sample_points[1]

# Inject control inputs into plant:
nmpc.plant.vars.input.set_value(nmpc.controller.vars.input[:,c_ts].value)
# (^ input would be a custom _NmpcVector object that supports setting
# its matrix (inputs x time) values with a vector (of inputs))
#
# The syntax 
# `nmpc.plant.vars.input[:,:].set_value(inputs)`
# would be nice, but then would need to decide which axis to
# broadcast `inputs` along...

# Simulate plant:
solver.solve(nmpc.plant, tee=True)
# (alternatively, initialize by element or simulate, then solve)

N = 20
for i in range(N):
    # Get measurements from plant and add noise:
    measured = list(nmpc.plant.vars.measured[:,p_ts])
    measured = nmpc.add_measurement_noise(measured)
    nmpc.plant.advance_time(sample_time)
    
    nmpc.controller.advance_time(sample_time)
    # Load measurements into controller
    nmpc.controller.vars.measured[:,t0].set_value(measured)
    # ^ Requires ability to set value of a slice with an iterable.

    solver.solve(nmpc.controller)

    # Get inputs from controller and add noise:
    inputs = list(nmpc.controller.vars.input[:,c_ts].value)
    inputs = nmpc.add_input_noise(inputs)

    # Inject control inputs into plant:
    nmpc.plant.vars.inputs.set_value(inputs)

    # Simulate plant:
    solver.solve(nmpc.plant, tee=True)
