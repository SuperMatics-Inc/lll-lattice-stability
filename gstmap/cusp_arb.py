"""Rigorous (Arb) evaluation of the cusp bound at t2 = T:  for all t1, Sigma_tau(s*,0) <= BOUND * c^2 < 0 (see PROOFS.md, Lemma C)."""
import sys
from flint import arb, ctx
ctx.prec = 200
def cusp_bound(T):
    T = arb(T)
    pi = arb.pi()
    a = (arb(10) / 7).log()
    x = (-a).exp()                                            # = 7/10
    sqrt_ = lambda z: z.sqrt()
    u = 1 - sqrt_(a / (pi * T))                               # 1 - s*
    delta = 2 * (-pi * T * u * u).exp() / (1 - (-2 * pi * T * u).exp())
    Tn = 2 * (-pi * T).exp() / (1 - (-3 * pi * T).exp())
    etaA = max((2 * delta + 2 * (x + delta) * Tn + Tn * Tn).upper(), (2 * (x + 1) * Tn + Tn * Tn).upper())
    etaA = arb(etaA)
    b = x * x - 2 * x * delta - delta * delta
    bound = (arb(2) * x - 1 + etaA) ** 2 - b ** 2
    return dict(x=x, delta=delta, Tn=Tn, etaA=etaA, b=b, bound=bound, ok=bool(bound < 0))
if __name__ == "__main__":
    fails=0
    for T in sys.argv[1:] or ["3", "5/2", "4", "10"]:
        from fractions import Fraction
        f = Fraction(T)
        r = cusp_bound(arb(f.numerator) / f.denominator)
        print("t2 =", T, {k: (v if isinstance(v, bool) else v.str(8)) for k, v in r.items()})
        fails+=int(not r['ok'])
    sys.exit(bool(fails))
