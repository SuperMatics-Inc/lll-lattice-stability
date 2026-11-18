"""Find box certificates for Sigma > 0 on T^2 except at the origin, or for instability."""
import sys, time
from fractions import Fraction
from flint import ctx, arb
import gsteng as G
from unst import try_unstable


def A(x):
    x = Fraction(x)
    return arb(x.numerator) / x.denominator


def rmax_for(w):
    """near radius shrinks with the box width (high-degree Taylor coefficients lose accuracy ~ w * rho^j)"""
    w = float(w)
    r = 0.3 * (1e-3 / w) ** 0.4 if w > 1e-3 else 0.3
    return Fraction(round(max(r, 0.04) * 1000), 1000)


def N_for_rmax(rmax, c2):
    r = float(rmax) * max(1.0, float(c2)) ** 0.5
    return 24 if r <= 0.2 else (30 if r <= 0.26 else 36)


def try_stable(c1, c2, w, M=6, prec=192, rin=Fraction(1, 100), g=Fraction(21, 20), rmax=None, N=None, do_far=True, verbose=False, ngrid=1):
    c1 = Fraction(c1); c2 = Fraction(c2); w = Fraction(w)
    ctx.prec = prec
    rmax = rmax_for(w) if rmax is None else Fraction(rmax)
    N = N_for_rmax(rmax, c2) if N is None else N
    P0, Pb = G.make_params(c1, c2, w, w, M, ngrid)
    lam = G.lam_AD(P0, Pb)
    t0 = time.time()
    nd = G.NearData(P0, Pb, w, w, lam, N)
    okN, stN = G.near_cover(nd, A(rmax), A(rin), g=A(g), verbose=verbose)
    tN = time.time() - t0
    info = dict(rmax=str(rmax), N=N, near=okN, tN=round(tN, 2), near_info={k: (None if v is None else str(v)[:30]) for k, v in stN.items()})
    if not okN or not do_far:
        return okN, info
    KAP6 = G.build_kappa6(Pb)
    okF, stF = G.far_cover(P0, Pb, w, w, lam, A(rmax), KAP6, hmin_exp=10, verbose=verbose, h0e=3)
    info.update(far=okF, tF=round(time.time() - t0 - tN, 2), far_ncells=stF.get('ncells'), far_fail=stF.get('fail'),
                far_worst=None if stF.get('worst') is None else str(stF['worst'])[:25])
    return okN and okF, info


if __name__ == "__main__":
    a = sys.argv
    ok, info = try_stable(Fraction(a[1]), Fraction(a[2]), Fraction(a[3]), verbose=False)
    print("STABLE-CERT", ok, info)
    sys.exit(0 if ok else 1)
