from pyomo.core.base.var import IndexedVar

NMPC_ATTRS = [
        'setpoint',
        'variance',
        'weight',
        'initial',
        'nominal',
        'reference', # Reference value, not a pointer to a reference object
        ]

# If I subclass Var here without overriding __new__,
# instances of this class will be SimpleVar or IndexedVar.
# This makes Var kinda tricky to subclass...

class NmpcVar(IndexedVar):
    def __init__(self, *args, **kwargs):
        if not args:
            raise NotImplementedError(
                '%s component must be indexed by at least one set.'
                % self.__class__
                )
        #for attr in NMPC_ATTRS:
        #    setattr(self, attr, kwargs.pop(attr, None))
        self.setpoint = kwargs.pop('setpoint', None)
        kwargs.setdefault('ctype', NmpcVar)
        super(NmpcVar, self).__init__(*args, **kwargs)

class _NmpcVector(IndexedVar):
    # TODO: Make sure I never add this monstrosity to an
    # active block in my model.
    # nmpc.controller.vars.input[:,:].fix() is still valid though
    # Then nmpc.controller.INPUT_BLOCK[i].var[t] still exists?
    # Should INPUT_BLOCK be attached?
    # Only reason would be so I can call component_objects(InputVar)
    # ^ This is probably a good idea, rather than having
    # component_objects(InputVector)

    def _generate_referenced_vars(self):
        _slice = self.referent.duplicate()
        # This class should only be instantiated as a reference...
        _slice._call_stack.pop()
        _slice._len -= 1
        for var in _slice:
            yield var
        # TODO: assert some properties of the call stack here
        # popped item is a getitem, slice
        # last remaining item is a getattr, NmpcVar

    def set_setpoint(self, setpoint):
        referent_gen = self._generate_referenced_vars()
        try:
            for var, sp in zip(referent_gen, setpoint):
                var.setpoint = sp
        except TypeError:
            for var in self._generate_referenced_vars():
                var.setpoint = setpoint

    def get_setpoint(self):
        for var in self._generate_referenced_vars():
            yield var.setpoint

    @property
    def values(self):
        referent_gen = self._generate_referenced_vars()
        return list(list(var[t].value for t in var) for var in referent_gen)

    @values.setter
    def values(self, vals):
        referent_gen = self._generate_referenced_vars()
        try:
            for var, val in zip(referent_gen, vals):
                # var is time-indexed
                var[:].set_value(val)
        except TypeError:
            # A scalar value was provided. Set for all i, t
            self[...].set_value(vals)
    # Pyomo seems to already have a similar functionality in
    # set_values, that sets values via a dict. Would it be
    # easier to do all my setting/getting via dicts?

    # But none of this really helps me set value or fix from an iterable...
    # Two syntaxes I would like to support:
    # var[:,:].set_value/fix(vector)
    # ^ Maybe this should be done through the new class...
    # var.set_value/fix(vector), automatically interpreting as
    # a vector of vars, eached set/fixed for all time
    #
    # But I rarely want to set/fix the entire var. More likely, I want:
    # var[:,t0:].set_value(some_vector)
    # (This is initialization, say to the setpoint)
    # or
    # var[:,t0:ts].set_value(some_vector # len n x n_t)
    # var[:].set_value(iter(some_vector))
    # var[t0:].set_value(setpoint)
    #
    # var[:,t0].set_value/fix(vector)
    # Then this could potentially be interpreted by the slice?
    # ^ This would be easier to implement than setting a matrix
    # from a vector.

# NOTE: .value attribute on a slice yields a slice
#       value function on a slice yields a list of values
#       behaves well even if value is None.
#       value function on a vardata, however, doesn't behave
#       well if the value is None.

# API I'd like:
# var[:,t].set_value(var[:,t0].value)
# var[:,t].value = var[:,t0].value
# var[:,t] = var[:,t0].value
# (^ set value from an iterable)

# var[:,:].set_value(var[:,t0].value)
# This is handled by the _NmpcVector class

# var[0,:].set_value(0)
# var[:,t0].set_value(0)
# ^ These will work currently if var is an IndexedVar

class DiffVar(NmpcVar):
    _attr = 'differential'

class DerivVar(NmpcVar):
    _attr = 'derivative'

class AlgVar(NmpcVar):
    _attr = 'algebraic'

class InputVar(NmpcVar):
    _attr = 'input'

class FixedVar(NmpcVar):
    _attr = 'fixed'

class MeasuredVar(NmpcVar):
    _attr = 'measurement'
