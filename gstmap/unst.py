"""Unstable-box certificates: a point q0 and a rigorous upper bound Sigma_tau(q0) < 0 for every tau in a box."""
import sys, time
from fractions import Fraction
import numpy as np
from flint import ctx, arb
import gsteng as G
import gstsym as S

RHOS = [0.02, 0.03, 0.045, 0.065, 0.09, 0.12, 0.16, 0.2]


def float_candidates(c1, c2, nth=361, rhos=RHOS, keep=4):
    """float64 guess: for each rho the minimising angle; returns list of (Sigma/rho^4, rho, s, xi) most negative first."""
    tau = complex(float(c1), float(c2))
    th = np.linspace(0, np.pi, nth)
    out = []
    # Small-q candidates must not evaluate A^2-|B|^2 in float arithmetic.
    # Taylor coefficients select a direction/radius only; Arb still proves each claim.
    pol=S.sigma_taylor(tau,6)[3]
    c4=S.homog(pol,4).real; c6=S.homog(pol,6).real
    co,si=np.cos(th),np.sin(th)
    p4=sum(c4[k]*co**(4-k)*si**k for k in range(5))
    p6=sum(c6[k]*co**(6-k)*si**k for k in range(7))
    rho=1/512
    v=p4+rho*rho*p6
    from scipy.optimize import minimize_scalar
    i=int(v.argmin())
    def value(theta,degree,coeff):
        return sum(coeff[k]*np.cos(theta)**(degree-k)*np.sin(theta)**k for k in range(degree+1))
    angle=minimize_scalar(lambda theta:value(theta,4,c4),
                          bounds=(th[i]-np.pi/nth,th[i]+np.pi/nth),method='bounded').x
    val=value(angle,4,c4)+rho*rho*value(angle,6,c6)
    out.append((float(val),rho,float(rho*np.cos(angle)),float(rho*np.sin(angle))))
    for r in rhos:
        v = np.asarray(S.sigma(tau, r * np.cos(th), r * np.sin(th))) / r ** 4
        i = int(v.argmin())
        out.append((float(v[i]), r, float(r * np.cos(th[i])), float(r * np.sin(th[i]))))
    out.sort()
    return out[:keep]


def rat(x, bits=24):
    return Fraction(round(x * 2 ** bits), 2 ** bits)


def N_for(rho, c2):
    x = rho * max(1.0, float(c2)) ** 0.5
    return 12 if x <= 0.11 else (18 if x <= 0.2 else 26)


def try_unstable(c1, c2, w1, w2, M=6, prec=192, ncand=3, verbose=False):
    """returns (ok, info). ok => Sigma_tau(q0)<0 for all tau in box, rigorous."""
    ctx.prec = prec
    cands = float_candidates(c1, c2)
    P0 = Pb = lam = None
    cache = {}
    best = None
    for (v, rho, s, xi) in cands[:ncand]:
        if v >= 0:
            continue
        N = N_for(rho, c2)
        if N not in cache:
            if P0 is None:
                P0, Pb = G.make_params(c1, c2, w1, w2, M)
                lam = G.lam_AD(P0, Pb)
            cache[N] = G.NearData(P0, Pb, w1, w2, lam, N)
        nd = cache[N]
        sq = rat(s); xq = rat(xi)
        sa = arb(sq.numerator) / sq.denominator
        xa = arb(xq.numerator) / xq.denominator
        up = nd.sigma_upper_at(sa, xa)
        if verbose:
            print("  cand rho", rho, "float S/rho^4", v, "N", N, "upper bound", up.str(8))
        if best is None or up < best[0]:
            best = (up, float(sq), float(xq), N)
        if up < 0:
            return True, dict(q0=(sq, xq), N=N, upper=float(up.mid()), rho=rho)
    return False, dict(best=None if best is None else (float(best[0].mid()), best[1], best[2], best[3]))


if __name__ == "__main__":
    a = sys.argv
    c1, c2, w1, w2 = Fraction(a[1]), Fraction(a[2]), Fraction(a[3]), Fraction(a[4])
    t0 = time.time()
    ok, info = try_unstable(c1, c2, w1, w2, verbose=True)
    print("UNSTABLE-CERT", ok, info, f"({time.time()-t0:.1f}s)")
    sys.exit(0 if ok else 1)
