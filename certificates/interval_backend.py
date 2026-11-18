"""Minimal, outward-rounded mpmath.iv adapter for replaying the Arb formulas.
This is an independent arithmetic backend, not an independent derivation of bounds.
Endpoint extraction uses exact binary mpf tuples, never float rounding.
"""
from mpmath import iv, mp, make_mpf


def point(t):
    return iv.mpf(make_mpf(t))


class arb:
    def __init__(self, x=0, radius=None):
        self.v = x.v if isinstance(x, arb) else iv.mpf(x)
        if radius is not None:
            r = abs(arb(radius)).v.b
            self.v = self.v + iv.mpf([-r, r])

    @classmethod
    def wrap(cls, x):
        z = cls.__new__(cls)
        z.v = x
        return z

    @staticmethod
    def pi(): return arb.wrap(+iv.pi)
    def lower(self): return arb.wrap(point(self.v._mpi_[0]))
    def upper(self): return arb.wrap(point(self.v._mpi_[1]))
    def mid(self): return (self.lower() + self.upper()) / 2
    def rad(self): return (self.upper() - self.lower()) / 2
    def abs_upper(self): return abs(self).upper()
    def abs_lower(self): return abs(self).lower()
    def contains(self, x):
        x = arb(x)
        return self.lower().v <= x.lower().v and x.upper().v <= self.upper().v
    def overlaps(self, x):
        x = arb(x)
        return self.lower().v <= x.upper().v and x.lower().v <= self.upper().v
    def sqrt(self): return arb.wrap(iv.sqrt(self.v))
    def exp(self): return arb.wrap(iv.exp(self.v))
    def cos(self): return arb.wrap(iv.cos(self.v))
    def sin(self): return arb.wrap(iv.sin(self.v))
    def __add__(self, x):
        if isinstance(x, acb): return acb(self) + x
        return arb.wrap(self.v + arb(x).v)
    __radd__ = __add__
    def __sub__(self, x):
        if isinstance(x, acb): return acb(self) - x
        return arb.wrap(self.v - arb(x).v)
    def __rsub__(self, x): return arb(x) - self
    def __mul__(self, x):
        if isinstance(x, acb): return x * self
        return arb.wrap(self.v * arb(x).v)
    __rmul__ = __mul__
    def __truediv__(self, x): return arb.wrap(self.v / arb(x).v)
    def __rtruediv__(self, x): return arb(x) / self
    def __pow__(self, n): return arb.wrap(self.v ** n)
    def __neg__(self): return arb.wrap(-self.v)
    def __abs__(self): return arb.wrap(abs(self.v))
    def __lt__(self, x): return bool(self.v < arb(x).v)
    def __le__(self, x): return bool(self.v <= arb(x).v)
    def __gt__(self, x): return bool(self.v > arb(x).v)
    def __ge__(self, x): return bool(self.v >= arb(x).v)
    def __float__(self): return float(make_mpf(self.mid().v._mpi_[0]))
    def str(self, *args, **kwargs): return str(self.v)
    def __repr__(self): return str(self.v)


class acb:
    def __init__(self, real=0, imag=0):
        if isinstance(real, acb):
            self.real, self.imag = real.real, real.imag
        else:
            self.real, self.imag = arb(real), arb(imag)
    def __add__(self, x):
        x = acb(x)
        return acb(self.real+x.real, self.imag+x.imag)
    __radd__ = __add__
    def __neg__(self): return acb(-self.real, -self.imag)
    def __sub__(self, x): return self + -acb(x)
    def __rsub__(self, x): return acb(x) - self
    def __mul__(self, x):
        x = acb(x)
        return acb(self.real*x.real-self.imag*x.imag, self.real*x.imag+self.imag*x.real)
    __rmul__ = __mul__
    def __truediv__(self, x):
        if not isinstance(x, acb): return acb(self.real/x, self.imag/x)
        den = x.real*x.real+x.imag*x.imag
        return self*x.conjugate()/den
    def exp(self):
        r = self.real.exp()
        return acb(r*self.imag.cos(), r*self.imag.sin())
    def conjugate(self): return acb(self.real, -self.imag)
    def abs_upper(self): return (self.real**2+self.imag**2).sqrt().upper()
    def __repr__(self): return f'({self.real}) + i({self.imag})'


class acb_mat:
    def __init__(self, rows, cols, data):
        self.rows, self.cols, self.data = rows, cols, data
    def __getitem__(self, k): return self.data[k[0]*self.cols+k[1]]
    def __mul__(self, other):
        assert self.cols == other.rows
        return acb_mat(self.rows, other.cols, [
            sum((self[i,k]*other[k,j] for k in range(self.cols)), acb(0))
            for i in range(self.rows) for j in range(other.cols)])


class Context:
    @property
    def prec(self): return iv.prec
    @prec.setter
    def prec(self, value):
        iv.prec = value
        mp.prec = value


ctx = Context()
