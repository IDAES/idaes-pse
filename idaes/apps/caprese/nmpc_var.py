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
        for attr in NMPC_ATTRS:
            setattr(self, attr, kwargs.pop(attr, None))
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
        _slice._call_stack.pop()
        _slice._len -= 1
        for var in _slice:
            yield var

    # What about...
    def __iter__(self):
        _slice = self.referent.duplicate()
        _slice._call_stack.pop()
        _slice._len -= 1
        for var in _slice:
            yield var
    # Will self[:,t0] still behave properly?

    def set_setpoint(self, setpoint):
        try:
            for var, sp in zip(self, setpoint):
                var.setpoint = sp
        except TypeError:
            for var in self:
                var.setpoint = setpoint
        #for var, sp in zip(self, setpoint):
        #    var.setpoint = sp
        # What if setpoint is not iterable, e.g. 
        # nmpc.controller.vars.derivative.set_setpoint(0.)

    # But none of this really helps me set value or fix from an iterable...
    # Two syntaxes I would like to support:
    # var[:,:].set_value/fix(vector)
    # ^ Maybe this should be done through the new class...
    # var.set_value/fix(vector), automatically interpreting as
    # a vector of vars, eached set/fixed for all time
    #
    # var[:,t0].set_value/fix(vector)
    # Then this could potentially be interpreted by the slice?

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

# var[0,:].set_value(0)
# var[:,t0].set_value(0)
# ^ These will work currently if var is an IndexedVar

# Want to be able to treat var as either a vector or a 
# matrix.
# var[:,:].set_value(var[:,t0])
# How do I know which dimension to broadcast here?
# Should var[0] return an IndexedVar?
# Default should be to iterate over first index?
# The "Pyomothonic" way of doing this would probably
# be to have different components for different dimensions.
# The "tuplized-index" components should be specifically for
# operating on all vars in one line
#
# How should this work for setpoint attributes?
# var.setpoint is a list of varlist[i].setpoint
# var[:].setpoint is undefined, because var[:] is undefined.
# var[:,t0].setpoint is undefined because a vardata has no
# setpoint
# var.set_setpoint(var.setpoint for var in varlist)

class DifferentialVar(NmpcVar):
    pass

class DerivativeVar(NmpcVar):
    # FIXME Find better name
    pass

class AlgebraicVar(NmpcVar):
    pass

class InputVar(NmpcVar):
    pass

class FixedVar(NmpcVar):
    pass

class MeasuredVar(NmpcVar):
    pass
