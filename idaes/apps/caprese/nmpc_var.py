from pyomo.core.base.var import IndexedVar

# If I subclass Var here without overriding __new__,
# instances of this class will be SimpleVar or IndexedVar.
# This makes Var somewhat tricky to subclass...

class NmpcVar(IndexedVar):
    def __init__(self, *args, **kwargs):
        if not args:
            raise NotImplementedError(
                '%s component must be indexed by at least one set.'
                % self.__class__
                )
        self.setpoint = kwargs.pop('setpoint', None)
        self.weight = kwargs.pop('weight', None)
        self.variance = kwargs.pop('variance', None)
        self.nominal = kwargs.pop('nominal', None)
        self.noise_bounds = kwargs.pop('noise_bounds', None)
        kwargs.setdefault('ctype', type(self))
        super(NmpcVar, self).__init__(*args, **kwargs)

class _NmpcVector(IndexedVar):

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
