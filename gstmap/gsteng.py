r"""Enclose the lowest-Landau-level symbol over boxes of lattice shapes with Arb.

Write Sigma(s,xi) = A^2 - |B|^2, with A = 2K - lam and K,B the Gaussian sums of the GST Hamiltonian.
Starting from the single-shape formulas in ../certificates/cert_sigma.py, enclose each coefficient by
S(tau) in S(tau0) + sum_i dS/dt_i(box) * [-w_i, w_i].

P0 gives parameters at the exact rational centre tau0, up to rounding. Pb encloses them on the closed box.
"""
import time
from math import comb, factorial
from fractions import Fraction
from flint import arb, acb, acb_mat, ctx

T0 = time.time()
def hdr(msg):
    print(f"[{time.time()-T0:7.1f}s] {msg}", flush=True)

def qarb(x):
    x = Fraction(x)
    return arb(x.numerator) / x.denominator if x.denominator != 1 else arb(x.numerator)

def ub(x):
    return x.abs_upper()

_TERMS = {}
def terms(M):
    if M not in _TERMS:
        _TERMS[M] = [(m, n) for m in range(-M, M + 1) for n in range(-M, M + 1)]
    return _TERMS[M]

class Params:
    """lattice parameter tau = T1 + i T2 (arb balls) and derived constants."""
    def __init__(self, T1, T2, M):
        self.T1 = T1; self.T2 = T2; self.M = M
        self.PI = arb.pi(); self.TWOPI = 2 * self.PI
        self.C0 = T2.sqrt() / self.PI
        self.TAUABS = (T1 * T1 + T2 * T2).sqrt()
        self.KAPPA = acb(-self.TWOPI * T2, -self.TWOPI * T1)
        self.cache = {}

def make_params(c1, c2, w1, w2, M, ngrid=1):
    """(P0, Pb): exact-centre params and box params (box = [c1-w1,c1+w1] x [c2-w2,c2+w2]); all Fractions.
    ngrid>1: Pb.grid lists the ngrid^2 sub-box params; derivative jets over Pb are then the hull of the sub-box enclosures
    (a sharper but still valid enclosure of the derivative over the whole box)."""
    P0 = Params(qarb(c1), qarb(c2), M)
    Pb = Params(qarb(c1) + arb(0, qarb(w1)), qarb(c2) + arb(0, qarb(w2)), M)
    Pb.grid = None
    if ngrid > 1:
        Pb.grid = []
        for i in range(ngrid):
            for j in range(ngrid):
                a1 = c1 - w1 + (2 * i + 1) * w1 / ngrid
                a2 = c2 - w2 + (2 * j + 1) * w2 / ngrid
                Pb.grid.append(Params(qarb(a1) + arb(0, qarb(w1 / ngrid)), qarb(a2) + arb(0, qarb(w2 / ngrid)), M))
    return P0, Pb

# ----------------------------------------------------------------------------------------------------------------------
# jets of K (kind 'A') and B (kind 'B'); der in {0,1,2}: value, d/dt1, d/dt2
# ----------------------------------------------------------------------------------------------------------------------
def term_data(P, kind, m, n, s0, xi0):
    PI, TWOPI, T1, T2 = P.PI, P.TWOPI, P.T1, P.T2
    half_inv_t2 = 1 / (2 * T2)
    if kind == 'A':
        x = s0 + m
        E0 = acb(-PI * T2 * (x * x + n * n), TWOPI * (T1 * n * x + xi0 * n))
        Es = acb(-TWOPI * T2 * x, TWOPI * T1 * n)
        e2 = acb(-PI * T2, 0)
        ex = acb(0, TWOPI * n)
        # d/dt1: 2 pi i n x  (x = x0 + u)       d/dt2: -pi (x^2 + n^2) + 1/(2 t2)
        mult1 = (acb(0, TWOPI * n * x), acb(0, TWOPI * n), acb(0))
        mult2 = (acb(-PI * (x * x + n * n) + half_inv_t2), acb(-TWOPI * x), acb(-PI))
    else:
        x = s0 + m; y = n - s0
        E0 = acb(-PI * T2 * (x * x + y * y), TWOPI * (T1 * x * y + xi0 * (n - m)))
        Es = P.KAPPA * (x - y)
        e2 = P.KAPPA
        ex = acb(0, TWOPI * (n - m))
        # x = x0+u, y = y0-u:  xy = x0y0 + (y0-x0)u - u^2 ; x^2+y^2 = x0^2+y0^2 + 2(x0-y0)u + 2u^2
        mult1 = (acb(0, TWOPI * x * y), acb(0, TWOPI * (y - x)), acb(0, -TWOPI))
        mult2 = (acb(-PI * (x * x + y * y) + half_inv_t2), acb(-TWOPI * (x - y)), acb(-TWOPI))
    return P.C0 * E0.exp(), Es, e2, ex, {1: mult1, 2: mult2}

def jets_multi(P, kind, s0, xi0, order, ders=(0,)):
    """box-truncated Taylor coefficients (acb) a+b<=order, for each requested derivative index; one shared lattice pass."""
    TT = terms(P.M); TN = len(TT)
    Fs = {d: [None] * ((order + 1) * TN) for d in ders}
    X = [None] * (TN * (order + 1))
    for idx, (m, n) in enumerate(TT):
        f0, Es, e2, ex, mult = term_data(P, kind, m, n, s0, xi0)
        c = [acb(1), Es]
        e22 = 2 * e2
        for a in range(2, order + 1):
            c.append((Es * c[a - 1] + e22 * c[a - 2]) / a)
        for d in ders:
            Fd = Fs[d]
            if d == 0:
                for a in range(order + 1):
                    Fd[a * TN + idx] = f0 * c[a]
            else:
                q0, q1, q2 = mult[d]
                for a in range(order + 1):
                    v = q0 * c[a]
                    if a >= 1:
                        v = v + q1 * c[a - 1]
                    if a >= 2:
                        v = v + q2 * c[a - 2]
                    Fd[a * TN + idx] = f0 * v
        xp = acb(1)
        for b in range(order + 1):
            X[idx * (order + 1) + b] = xp
            xp = xp * ex / (b + 1)
    Xm = acb_mat(TN, order + 1, X)
    out = {}
    for d in ders:
        C = acb_mat(order + 1, TN, Fs[d]) * Xm
        out[d] = {(a, b): C[a, b] for a in range(order + 1) for b in range(order + 1 - a)}
    return out

# ----------------------------------------------------------------------------------------------------------------------
# rigorous lattice-tail bounds (shells t = max(|m|,|n|) > M)
# ----------------------------------------------------------------------------------------------------------------------
def bound_poly(a, b, Ebar, e2b, exb):
    s = arb(0)
    for j in range(a // 2 + 1):
        s += e2b ** j * Ebar ** (a - 2 * j) / (factorial(j) * factorial(a - 2 * j))
    return s * exb ** b / factorial(b)

def shell_term(P, a, b, t, der):
    t = arb(t)
    Ebar = P.TWOPI * P.TAUABS * (2 * t + 1); e2b = P.TWOPI * P.TAUABS; exb = 4 * P.PI * t
    base = 8 * t * P.C0 * (-P.PI * P.T2 * (t - arb(1) / 2) ** 2).exp()
    if der == 0:
        return base * bound_poly(a, b, Ebar, e2b, exb)
    Q0 = P.TWOPI * (t + 1) ** 2 + 1 / (2 * P.T2); Q1 = P.TWOPI * (2 * t + 1); Q2 = P.TWOPI
    s = Q0 * bound_poly(a, b, Ebar, e2b, exb)
    if a >= 1:
        s += Q1 * bound_poly(a - 1, b, Ebar, e2b, exb)
    if a >= 2:
        s += Q2 * bound_poly(a - 2, b, Ebar, e2b, exb)
    return base * s

def shell_tail(P, a, b, der, extra=30):
    tot = arb(0)
    M = P.M
    for t in range(M + 1, M + extra + 1):
        tot += shell_term(P, a, b, t, der)
    T = M + extra
    rho = (arb(T + 1) / T) ** (a + b + 3) * (-P.TWOPI * P.T2 * T).exp()
    assert rho < arb(1) / 2, "tail ratio"
    return (tot + shell_term(P, a, b, T + 1, der) / (1 - rho)).upper()

def tailJ(P, a, b, der=0):
    k = ('tail', a, b, der)
    if k not in P.cache:
        P.cache[k] = shell_tail(P, a, b, der)
    return P.cache[k]

def infl(z, r):
    return z + acb(arb(0, r), arb(0, r))

def _hull_acb(zs):
    re = zs[0].real; im = zs[0].imag
    for z in zs[1:]:
        re = re.union(z.real); im = im.union(z.imag)
    return acb(re, im)

def jets_full(P, kind, s0, xi0, order, ders=(0,)):
    grid = getattr(P, 'grid', None)
    if grid:
        outs = [jets_full(Q, kind, s0, xi0, order, ders) for Q in grid]
        return {d: {k: _hull_acb([o[d][k] for o in outs]) for k in outs[0][d]} for d in ders}
    J = jets_multi(P, kind, s0, xi0, order, ders)
    return {d: {k: infl(v, tailJ(P, k[0], k[1], d)) for k, v in J[d].items()} for d in ders}

# ----------------------------------------------------------------------------------------------------------------------
# polynomial helpers
# ----------------------------------------------------------------------------------------------------------------------
def pmul(P_, Q, deg):
    R = {}
    for (a, b), x in P_.items():
        for (c, d), y in Q.items():
            if a + b + c + d <= deg:
                k = (a + c, b + d)
                R[k] = R[k] + x * y if k in R else x * y
    return R

def real_of(v):
    return v.real if isinstance(v, acb) else v

def padd(P_, Q):
    R = dict(P_)
    for k, v in Q.items():
        R[k] = R[k] + v if k in R else v
    return R

def sigma_pt_poly(Apt, Bpt, deg):
    AA = pmul(Apt, Apt, deg)
    Bq = {k: acb(v) if not isinstance(v, acb) else v for k, v in Bpt.items()}
    Bb = {k: v.conjugate() for k, v in Bq.items()}
    BB = pmul(Bq, Bb, deg)
    S = {}
    for k in set(AA) | set(BB):
        S[k] = (acb(AA.get(k, arb(0))) - BB.get(k, acb(0))).real
    return S

# ----------------------------------------------------------------------------------------------------------------------
# lambda (= K(0,0)) value and tau-derivatives
# ----------------------------------------------------------------------------------------------------------------------
def lam_AD(P0, Pb):
    Jp = jets_full(P0, 'B', arb(0), arb(0), 0)[0][(0, 0)].real
    Jd = jets_full(Pb, 'B', arb(0), arb(0), 0, (1, 2))
    return Jp, Jd[1][(0, 0)].real, Jd[2][(0, 0)].real

# ----------------------------------------------------------------------------------------------------------------------
# mean-value (centred) enclosure of the Taylor coefficients S_ab of Sigma at (s0,xi0), a+b<=order (and products up to sdeg)
# ----------------------------------------------------------------------------------------------------------------------
def sigma_box(P0, Pb, w1, w2, s0, xi0, order, lam, sdeg=None, info=None):
    """returns dict (a,b)-> arb enclosing the (a,b) Taylor coefficient of Sigma(.;tau) at (s0,xi0) for ALL tau in the box."""
    if sdeg is None:
        sdeg = order
    I1 = arb(0, qarb(w1)); I2 = arb(0, qarb(w2))
    lam_pt, lam_d1, lam_d2 = lam
    JKp = jets_full(P0, 'A', s0, xi0, order)[0]
    JBp = jets_full(P0, 'B', s0, xi0, order)[0]
    JKd = jets_full(Pb, 'A', s0, xi0, order, (1, 2))
    JBd = jets_full(Pb, 'B', s0, xi0, order, (1, 2))
    Apt = {}; Ad = {1: {}, 2: {}}
    for k, v in JKp.items():
        Apt[k] = 2 * v.real - (lam_pt if k == (0, 0) else 0)
    for d, lamd in ((1, lam_d1), (2, lam_d2)):
        for k, v in JKd[d].items():
            Ad[d][k] = 2 * v.real - (lamd if k == (0, 0) else 0)
    Bpt = dict(JBp); Bd = {1: dict(JBd[1]), 2: dict(JBd[2])}
    Aball = {k: Apt[k] + Ad[1][k] * I1 + Ad[2][k] * I2 for k in Apt}
    Bball = {k: Bpt[k] + Bd[1][k] * I1 + Bd[2][k] * I2 for k in Bpt}
    Bbc = {k: v.conjugate() for k, v in Bball.items()}
    Spt = sigma_pt_poly(Apt, Bpt, sdeg)
    Sd = {}
    for d in (1, 2):
        t1 = pmul(Aball, Ad[d], sdeg)
        t2 = pmul(Bbc, Bd[d], sdeg)
        Sd[d] = {k: 2 * t1.get(k, arb(0)) - 2 * t2.get(k, acb(0)).real for k in set(t1) | set(t2)}
    S = {}
    for k, v in Spt.items():
        S[k] = v + Sd[1].get(k, arb(0)) * I1 + Sd[2].get(k, arb(0)) * I2
    if info is not None:
        info['Aball'] = Aball; info['Bball'] = Bball
        info['aff'] = {k: (Spt[k], Sd[1].get(k, arb(0)), Sd[2].get(k, arb(0))) for k in Spt}
    return S

def p4_box(P0, Pb, w1, w2, lam, info=None):
    """box enclosure of the degree-4 coefficient polynomial of Sigma via the closed form (Prop. 3.1 of the paper)
        P4 = 2 lam (A4 - Re B4) - (Im B2)^2,   A4 = 2 K4 (degree-4 part),
    with a centred mean-value form applied to the *closed form* (so that the large cancelling lattice sums A4, Re B4 are
    differenced before the tau-box enters).  dict (a,b) -> arb, a+b = 4."""
    I1 = arb(0, qarb(w1)); I2 = arb(0, qarb(w2))
    JKp = jets_full(P0, 'A', arb(0), arb(0), 4)[0]; JBp = jets_full(P0, 'B', arb(0), arb(0), 4)[0]
    JKd = jets_full(Pb, 'A', arb(0), arb(0), 4, (1, 2)); JBd = jets_full(Pb, 'B', arb(0), arb(0), 4, (1, 2))
    lam_pt, lam_d1, lam_d2 = lam
    lam_ball = lam_pt + lam_d1 * I1 + lam_d2 * I2
    def Dq(JK, JB, k):
        return 2 * JK[k].real - JB[k].real
    k4 = [(a, 4 - a) for a in range(5)]
    k2 = [(a, 2 - a) for a in range(3)]
    D_pt = {k: Dq(JKp, JBp, k) for k in k4}
    D_d = {i: {k: Dq(JKd[i], JBd[i], k) for k in k4} for i in (1, 2)}
    D_ball = {k: D_pt[k] + D_d[1][k] * I1 + D_d[2][k] * I2 for k in k4}
    M_pt = {k: JBp[k].imag for k in k2}
    M_d = {i: {k: JBd[i][k].imag for k in k2} for i in (1, 2)}
    M_ball = {k: M_pt[k] + M_d[1][k] * I1 + M_d[2][k] * I2 for k in k2}
    Fpt = {k: 2 * lam_pt * D_pt[k] for k in k4}
    sq = pmul(M_pt, M_pt, 4)
    for k in k4:
        Fpt[k] = Fpt[k] - sq.get(k, arb(0))
    Fd = {}
    lam_d = {1: lam_d1, 2: lam_d2}
    for i in (1, 2):
        cr = pmul(M_ball, M_d[i], 4)
        Fd[i] = {k: 2 * lam_d[i] * D_ball[k] + 2 * lam_ball * D_d[i][k] - 2 * cr.get(k, arb(0)) for k in k4}
    if info is not None:
        info['aff4'] = {k: (Fpt[k], Fd[1][k], Fd[2][k]) for k in k4}
    return {k: Fpt[k] + Fd[1][k] * I1 + Fd[2][k] * I2 for k in k4}

# ----------------------------------------------------------------------------------------------------------------------
# global sup bound of 6th derivatives over the full torus for tau in the box (kappa6), as in the paper
# ----------------------------------------------------------------------------------------------------------------------
DEG6 = 6
def global_U(P, kind, a, b):
    tot = arb(0)
    half = arb(1) / 2
    PI, TWOPI, T1, T2 = P.PI, P.TWOPI, P.T1, P.T2
    for (m, n) in terms(P.M):
        am, an = abs(m), abs(n)
        dm = max(am - half, arb(0)) if am > 0 else arb(0)
        if kind == 'A':
            ex_ = (-PI * T2 * (dm * dm + an * an)).exp()
            Ebar = TWOPI * T2 * (am + half) + TWOPI * abs(T1) * an
            e2b = PI * T2; exb = TWOPI * an
        else:
            dn = max(an - half, arb(0)) if an > 0 else arb(0)
            ex_ = (-PI * T2 * (dm * dm + dn * dn)).exp()
            Ebar = TWOPI * P.TAUABS * (abs(m - n) + 1); e2b = TWOPI * P.TAUABS; exb = TWOPI * abs(n - m)
        tot += P.C0 * ex_ * bound_poly(a, b, Ebar, e2b, exb)
    tot += tailJ(P, a, b, 0)
    return (factorial(a) * factorial(b) * tot).upper()

def build_kappa6(P):
    UK = {}; UB = {}
    for a in range(DEG6 + 1):
        for b in range(DEG6 + 1 - a):
            UK[(a, b)] = global_U(P, 'A', a, b)
            UB[(a, b)] = global_U(P, 'B', a, b)
    UA = {k: 2 * v for k, v in UK.items()}
    UA[(0, 0)] = UA[(0, 0)] + UK[(0, 0)]
    kap = arb(0)
    for a in range(DEG6 + 1):
        b = DEG6 - a
        s = arb(0)
        for i in range(a + 1):
            for j in range(b + 1):
                s += comb(a, i) * comb(b, j) * (UA[(i, j)] * UA[(a - i, b - j)] + UB[(i, j)] * UB[(a - i, b - j)])
        kap += s / (factorial(a) * factorial(b))
    return kap.upper()

# ----------------------------------------------------------------------------------------------------------------------
# far-region cells
# ----------------------------------------------------------------------------------------------------------------------
ORD = 5
def cell_low(P0, Pb, w1, w2, lam, s0, xi0, h, KAP6):
    S = sigma_box(P0, Pb, w1, w2, s0, xi0, ORD, lam)
    low = S[(0, 0)]
    hp = arb(1)
    for d in range(1, ORD + 1):
        hp = hp * h
        sd = arb(0)
        for a in range(d + 1):
            sd += ub(S[(a, d - a)])
        low = low - sd * hp
    return low - KAP6 * h ** 6

def far_cover(P0, Pb, w1, w2, lam, rmax, KAP6, hmin_exp=9, verbose=True, h0e=3):
    """adaptive quadtree over [0,1/2]x[-1/2,1/2] minus the disc |q|<=rmax; returns (ok, stats)."""
    hmin = arb(1) / (2 ** hmin_exp)
    half = arb(1) / 2
    h0 = arb(1) / (2 ** h0e)
    stack = []
    nx = 2 ** (h0e - 2)    # cells along s in [0,1/2] (each of width 2*h0)
    for i in range(nx):
        for j in range(2 * nx):
            stack.append((h0 + 2 * h0 * i, -half + h0 + 2 * h0 * j, h0))
    passed = {}; skipped = 0; ncells = 0; worst = None
    rm2 = rmax * rmax
    while stack:
        s0, xi0, h = stack.pop()
        fs = s0 + h; fx = abs(xi0) + h
        if (fs * fs + fx * fx) <= rm2:
            skipped += 1
            continue
        low = cell_low(P0, Pb, w1, w2, lam, s0, xi0, h, KAP6)
        ncells += 1
        if low > 0:
            key = float(h.mid())
            passed[key] = passed.get(key, 0) + 1
            lw = low.lower()
            worst = lw if worst is None or lw < worst else worst
        else:
            if h <= hmin:
                if verbose:
                    print("FAIL at cell", float(s0.mid()), float(xi0.mid()), float(h.mid()), "lower bound", low)
                return False, dict(fail=(float(s0.mid()), float(xi0.mid()), float(h.mid())), ncells=ncells)
            h2 = h / 2
            for ds in (-1, 1):
                for dx in (-1, 1):
                    stack.append((s0 + ds * h2, xi0 + dx * h2, h2))
    return True, dict(ncells=ncells, passed=dict(sorted(passed.items())), skipped=skipped, worst=worst)

# ----------------------------------------------------------------------------------------------------------------------
# near region: polar rings with angular cells
# ----------------------------------------------------------------------------------------------------------------------
def taylor_tail_eps(P, rho, N):
    """rigorous bounds epsK, epsB for the Taylor remainder (total degree > N) of K and B on |s|,|xi| <= rho at tau in the
    (ball) parameter P, over all lattice terms: box terms via exact degree-weighted exponential tails, shells t>M whole."""
    PI, TWOPI, T1, T2 = P.PI, P.TWOPI, P.T1, P.T2
    epsK = arb(0); epsB = arb(0)
    rho_u = rho.upper()
    def wtail(alpha, beta):
        # sum_{d>N} p_d rho^d  where sum p_d z^d = exp(alpha z + beta z^2); monotone in alpha,beta: use exact upper bounds
        alpha = alpha.upper(); beta = beta.upper()
        p = [arb(1), alpha]
        for d in range(2, N + 1):
            p.append((alpha * p[d - 1] + 2 * beta * p[d - 2]) / d)
        part = arb(0); rp = arb(1)
        for d in range(N + 1):
            part += p[d] * rp
            rp = rp * rho_u
        y = alpha * rho_u + beta * rho_u * rho_u
        return (y.exp() - part).upper()
    for (m, n) in terms(P.M):
        an = abs(n)
        f0 = P.C0 * (-PI * T2 * (m * m + n * n)).exp()
        e1 = ub(acb(-TWOPI * T2 * m, TWOPI * T1 * n)); ex = TWOPI * an; e2 = PI * T2
        epsK += f0 * wtail(e1 + ex, e2)
        e1 = TWOPI * P.TAUABS * abs(m - n); ex = TWOPI * abs(m - n); e2 = TWOPI * P.TAUABS
        epsB += f0 * wtail(e1 + ex, e2)
    def shell_eps(t):
        t = arb(t)
        alpha = (TWOPI * P.TAUABS * 2 + 4 * PI) * t + TWOPI * P.TAUABS
        y = alpha * rho + TWOPI * P.TAUABS * rho * rho
        return 8 * t * P.C0 * (-PI * T2 * (t - arb(1) / 2) ** 2).exp() * y.exp()
    sh = arb(0)
    M = P.M
    for t in range(M + 1, M + 31):
        sh += shell_eps(t)
    T = M + 30
    rr = (arb(T + 1) / T) * (-TWOPI * T2 * T).exp() * ((TWOPI * P.TAUABS * 2 + 4 * PI) * rho).exp()
    assert rr < arb(1) / 2, "near tail ratio"
    sh = (sh + shell_eps(T + 1) / (1 - rr)).upper()
    return (epsK + sh).upper(), (epsB + sh).upper()

class NearData:
    def __init__(self, P0, Pb, w1, w2, lam, N):
        self.N = N
        info = {}
        self.S = sigma_box(P0, Pb, w1, w2, arb(0), arb(0), N, lam, sdeg=2 * N, info=info)
        self.Pb = Pb
        # low-order coefficients are exactly 0 (analytic, Prop 3.1 of the paper); enclosures only checked to contain 0
        for (a, b), v in self.S.items():
            if a + b <= 3:
                assert arb(0).overlaps(v), f"low-order coefficient {(a, b)} excludes 0: {v}"
        # degree-4 coefficients: replaced by the (much sharper) closed-form enclosure; consistency with the generic one checked
        info4 = {}
        self.Q4 = p4_box(P0, Pb, w1, w2, lam, info=info4)
        for k, v in self.Q4.items():
            assert v.overlaps(self.S[k]), f"closed-form P4 inconsistent with generic S_4 at {k}: {v} vs {self.S[k]}"
        self.byj = {}
        self.aff = {}
        for (a, b), v in self.S.items():
            j = a + b
            if j >= 5:
                self.byj.setdefault(j, []).append((a, b, v))
                self.aff.setdefault(j, []).append((a, b) + info['aff'][(a, b)])
        self.byj[4] = [(a, b, v) for (a, b), v in self.Q4.items()]
        self.aff[4] = [(a, b) + info4['aff4'][(a, b)] for (a, b) in self.Q4]
        self.I1 = arb(0, qarb(w1)); self.I2 = arb(0, qarb(w2))
        self.jmax = 2 * N
        self.Acoef = info['Aball']
        self.Bcoef = info['Bball']
        self.eps_cache = {}

    def etot(self, rho):
        key = float(rho.mid())
        if key in self.eps_cache:
            return self.eps_cache[key]
        epsK, epsB = taylor_tail_eps(self.Pb, rho, self.N)
        epsA = 2 * epsK
        Abar = arb(0); Bbar = arb(0)
        for (a, b), v in self.Acoef.items():
            Abar += ub(v) * rho ** (a + b)
        for (a, b), v in self.Bcoef.items():
            Bbar += ub(v) * rho ** (a + b)
        E = (epsA * (2 * Abar + epsA) + epsB * (2 * Bbar + epsB)).upper()
        self.eps_cache[key] = E
        return E

    def lo_j(self, cth, sth):
        """lower bounds lo[j] (arb) of S_j(Theta;tau) = sum_{a+b=j} S_ab(tau) cos^a sin^b over the angle ball and ALL tau in the box,
        using the centred form  c_j(Theta) + D1*G1_j(Theta) + D2*G2_j(Theta)  (same D_i in all monomials, so the first-order
        dependence of the box on the direction is not lost)."""
        cp = [arb(1)]; sp = [arb(1)]
        for i in range(self.jmax):
            cp.append(cp[-1] * cth); sp.append(sp[-1] * sth)
        lo = {}
        for j, lst in self.aff.items():
            c = arb(0); g1 = arb(0); g2 = arb(0)
            for (a, b, v0, d1, d2) in lst:
                m = cp[a] * sp[b]
                c += v0 * m; g1 += d1 * m; g2 += d2 * m
            lo[j] = (c + self.I1 * g1 + self.I2 * g2).lower()
        return lo

    def ring_margin(self, lo, ra, rb, Eb):
        """margin of  Sigma >= rho^4 * F  on ra<=rho<=rb for angle cell with bounds lo"""
        F = lo[4]
        for j in range(5, self.jmax + 1):
            lj = lo[j]
            if lj > 0:
                F += ra ** (j - 4) * lj
            else:
                F += rb ** (j - 4) * lj
        return F - Eb / ra ** 4

    def sigma_upper_at(self, s, xi):
        """rigorous upper bound of Sigma_tau(s,xi) valid for every tau in the box; (s,xi) exact arb points (rational)."""
        rho = (s * s + xi * xi).sqrt().upper()
        sp = [arb(1)]; xp = [arb(1)]
        for i in range(self.jmax):
            sp.append(sp[-1] * s); xp.append(xp[-1] * xi)
        c = arb(0); g1 = arb(0); g2 = arb(0)
        for j, lst in self.aff.items():
            for (a, b, v0, d1, d2) in lst:
                m = sp[a] * xp[b]
                c += v0 * m; g1 += d1 * m; g2 += d2 * m
        return (c + self.I1 * g1 + self.I2 * g2 + self.etot(rho)).upper()

def near_cover(nd, rmax, rin, g=arb(5) / 4, depth0=7, depth_max=16, verbose=True):
    """certify Sigma>0 on 0<|q|<=rmax via angular cells [0,pi) x rings; inner disc (0,rin] handled in one piece."""
    # ring boundaries: rin, rin*g, ..., >= rmax
    rings = []
    ra = rin
    while True:
        rb = ra * g
        if rb >= rmax:
            rb = rmax
            rings.append((ra, rb))
            break
        rings.append((ra, rb))
        ra = rb
    Eb = {}
    def Efor(r):
        k = float(r.mid())
        if k not in Eb:
            Eb[k] = nd.etot(r)
        return Eb[k]
    Ein = Efor(rin)
    pi = arb.pi()
    stack = [(i, depth0) for i in range(2 ** depth0)]
    ncell = 0; minmargin = None; fail = None
    while stack:
        i, dep = stack.pop()
        thm = pi * (2 * i + 1) / 2 ** (dep + 1)
        th = thm + arb(0, pi.upper() / 2 ** (dep + 1))
        lo = nd.lo_j(th.cos(), th.sin())
        ncell += 1
        # inner disc
        F = lo[4]
        for j in range(5, nd.jmax + 1):
            lj = lo[j]
            if lj < 0:
                F += rin ** (j - 4) * lj
        Fin = F - Ein / rin ** 4
        ok = bool(Fin > 0)
        mm = Fin
        if ok:
            for (ra, rb) in rings:
                m_ = nd.ring_margin(lo, ra, rb, Efor(rb))
                if not (m_ > 0):
                    ok = False
                    break
                mm = m_ if m_ < mm else mm
        if ok:
            if minmargin is None or mm.lower() < minmargin:
                minmargin = mm.lower()
        else:
            if dep >= depth_max:
                fail = (i, dep)
                break
            stack.append((2 * i, dep + 1)); stack.append((2 * i + 1, dep + 1))
    return fail is None, dict(ncell=ncell, minmargin=minmargin, fail=fail, nrings=len(rings))
