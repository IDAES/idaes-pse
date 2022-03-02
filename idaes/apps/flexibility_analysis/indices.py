class _ConIndex(object):
    def __init__(self, con, bound):
        self._con = con
        self._bound = bound

    @property
    def con(self):
        return self._con

    @property
    def bound(self):
        return self._bound

    def __repr__(self):
        if self.bound is None:
            return str(self.con)
        else:
            return str((str(self.con), str(self.bound)))

    def __str__(self):
        return repr(self)

    def __eq__(self, other):
        if isinstance(other, _ConIndex):
            return self.con is other.con and self.bound is other.bound
        return False

    def __hash__(self):
        return hash((self.con, self.bound))


class _VarIndex(object):
    def __init__(self, var, bound):
        self._var = var
        self._bound = bound

    @property
    def var(self):
        return self._var

    @property
    def bound(self):
        return self._bound

    def __repr__(self):
        if self.bound is None:
            return str(self.var)
        else:
            return str((str(self.var), str(self.bound)))

    def __str__(self):
        return repr(self)

    def __eq__(self, other):
        if isinstance(other, _VarIndex):
            return self.var is other.var and self.bound is other.bound
        return False

    def __hash__(self):
        return hash((id(self.var), self.bound))
