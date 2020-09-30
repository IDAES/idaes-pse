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

class NmpcVector(IndexedVar):
    pass
