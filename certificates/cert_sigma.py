r"""Use Arb to verify Sigma(s,xi) = A^2 - |B|^2 > 0 on the Brillouin torus T^2 \ {0}
for rational lattice parameters t1 and t2, and Sigma = P4 + O(|k|^6) with P4 > 0 near the origin.

Closed Gaussian lattice sums (derived in the report / paper; checked against the LLL kernel to 1e-16 in closed_form.py):
    G(x,y) = c exp(-pi t2 (x^2+y^2) + 2 pi i t1 x y),  c = sqrt(t2)/pi,
    K(s,xi) = sum_{(m,n) in Z^2} G(m+s, n)     e^{2 pi i xi n}                    (real),   A = 2K - lam,  lam = K(0,0)
    B(s,xi) = sum_{(m,n) in Z^2} G(m+s, n-s)   e^{2 pi i xi (n-m)}.
Sigma is 1-periodic in s and xi (|B| periodic) and even: Sigma(-s,-xi) = Sigma(s,xi)  (B(-s,-xi)=B(s,xi) by relabelling).
Strategy
  (I)  near region  rho := max(|s|,|xi|) <= r1 :  Taylor polynomial of A,B to degree N at 0 (coefficients enclosed incl. lattice tails),
       Section 4 of the report proves that the Sigma coefficients of degree <= 3 vanish exactly. We check that
       the intervals contain zero. Minimizing on the boundary of the sup-norm unit square gives P4 > p4lo * rho^4.
       The inequality p4lo - sum_{j>=5} c_j r1^{j-4} - (Taylor tail)/r1^4 > 0 suffices for every rho <= r1 by monotonicity.
  (II) far region: adaptive subdivision of [0,1/2] x [-1/2,1/2] (Sigma even => half torus), degree-5 jets of A,B at cell centres,
       degree-6 remainder from a global sup bound on the 6th derivatives of A and B, test
       Sigma(c) - sum_{1<=d<=5} h^d sum_{a+b=d}|S_ab| - kappa6 h^6 > 0 on each cell; cells inside the near box [-r1,r1]^2 are skipped.
Usage: python cert_sigma.py T1 T2 [r1=1/20] [N=20] [M=8] [prec=192] [hmin_exp=9]      (T1,T2 rationals like -2/5 19/20)
"""
import sys, time
from math import comb, factorial
from flint import arb, acb, acb_mat, ctx

def parse_q(s):
    a, _, b = s.partition('/')
    return arb(int(a)) / int(b) if b else arb(int(a))

N = int(sys.argv[4]) if len(sys.argv) > 4 else 20
M = int(sys.argv[5]) if len(sys.argv) > 5 else 8
PREC = int(sys.argv[6]) if len(sys.argv) > 6 else 192
HMIN_EXP = int(sys.argv[7]) if len(sys.argv) > 7 else 9
ctx.prec = PREC
T1 = parse_q(sys.argv[1]); T2 = parse_q(sys.argv[2]); R1 = parse_q(sys.argv[3]) if len(sys.argv) > 3 else arb(1) / 20

PI = arb.pi()
C0 = T2.sqrt() / PI
TAUABS = (T1 * T1 + T2 * T2).sqrt()
TWOPI = 2 * PI
TERMS = [(m, n) for m in range(-M, M + 1) for n in range(-M, M + 1)]
TN = len(TERMS)
ORD = 5          # jet order in cells
DEG6 = 6         # remainder order
T0 = time.time()

def hdr(msg):
    print(f"[{time.time()-T0:7.1f}s] {msg}", flush=True)

def ub(x):       # rigorous upper bound (arb) of |x| for arb/acb
    return x.abs_upper()

# ------------------------------------------------------------------------------------------------------------------
# exact (ball) term data
# ------------------------------------------------------------------------------------------------------------------
KAPPA = acb(-TWOPI * T2, -TWOPI * T1)          # -2 pi (t2 + i t1)

def term_data(kind, m, n, s0, xi0):
    if kind == 'A':
        x = s0 + m
        E0 = acb(-PI * T2 * (x * x + n * n), TWOPI * (T1 * n * x + xi0 * n))
        Es = acb(-TWOPI * T2 * x, TWOPI * T1 * n)
        e2 = acb(-PI * T2, 0)
        ex = acb(0, TWOPI * n)
    else:
        x = s0 + m; y = n - s0
        E0 = acb(-PI * T2 * (x * x + y * y), TWOPI * (T1 * x * y + xi0 * (n - m)))
        Es = KAPPA * (x - y)
        e2 = KAPPA
        ex = acb(0, TWOPI * (n - m))
    return C0 * E0.exp(), Es, e2, ex

def jets_raw(kind, s0, xi0, order):
    """box-truncated Taylor coefficients C[(a,b)] (acb), a+b<=order, of K (kind 'A') or B (kind 'B') at (s0,xi0)."""
    F = [None] * ((order + 1) * TN); X = [None] * (TN * (order + 1))
    for idx, (m, n) in enumerate(TERMS):
        f0, Es, e2, ex = term_data(kind, m, n, s0, xi0)
        c = [acb(1), Es]
        e22 = 2 * e2
        for a in range(2, order + 1):
            c.append((Es * c[a - 1] + e22 * c[a - 2]) / a)
        for a in range(order + 1):
            F[a * TN + idx] = f0 * c[a]
        xp = acb(1)
        for b in range(order + 1):
            X[idx * (order + 1) + b] = xp
            xp = xp * ex / (b + 1)
    C = acb_mat(order + 1, TN, F) * acb_mat(TN, order + 1, X)
    return {(a, b): C[a, b] for a in range(order + 1) for b in range(order + 1 - a)}

# ------------------------------------------------------------------------------------------------------------------
# rigorous lattice-tail bounds (terms with max(|m|,|n|) = t > M; 8t lattice points on that shell)
#   |f_mn| <= c exp(-pi t2 (t-1/2)^2),  |E_s| <= 2 pi|tau| (2t+1),  |e2| <= 2 pi|tau|,  |e_xi| <= 4 pi t   (|s|<=1/2)
#   (d/ds)^a (d/dxi)^b f = f * a! [t^a] exp(E_s t + e2 t^2) * e_xi^b   =>  |.| <= |f| a! b! bound(a,b,t)
# ------------------------------------------------------------------------------------------------------------------
def bound_poly(a, b, Ebar, e2b, exb):
    s = arb(0)
    for j in range(a // 2 + 1):
        s += e2b ** j * Ebar ** (a - 2 * j) / (factorial(j) * factorial(a - 2 * j))
    return s * exb ** b / factorial(b)

def shell_term(a, b, t):
    t = arb(t)
    Ebar = TWOPI * TAUABS * (2 * t + 1); e2b = TWOPI * TAUABS; exb = 4 * PI * t
    return 8 * t * C0 * (-PI * T2 * (t - arb(1) / 2) ** 2).exp() * bound_poly(a, b, Ebar, e2b, exb)

def shell_tail(a, b, extra=30):
    """sum_{t>M} 8t c e^{-pi t2 (t-1/2)^2} bound(a,b,t) (bound = coefficient bound without a!b!); rigorous upper bound."""
    tot = arb(0)
    for t in range(M + 1, M + extra + 1):
        tot += shell_term(a, b, t)
    T = M + extra
    rho = (arb(T + 1) / T) ** (a + b + 1) * (-TWOPI * T2 * T).exp()     # ratio bound for t>=T (poly deg <= a+b+1, gaussian factor)
    assert rho < arb(1) / 2, "tail ratio"
    return (tot + shell_term(a, b, T + 1) / (1 - rho)).upper()

# jets of the full-lattice sums: box sum + [-tail,tail]-ball.  tailJ(a,b) bounds coefficient (a,b) of each of K and B outside the box.
TAILJ = {}
def tailJ(a, b):
    if (a, b) not in TAILJ:
        TAILJ[(a, b)] = shell_tail(a, b)
    return TAILJ[(a, b)]

def infl(z, r):
    return z + acb(arb(0, r), arb(0, r))

def jets_full(kind, s0, xi0, order):
    J = jets_raw(kind, s0, xi0, order)
    return {k: infl(v, tailJ(*k)) for k, v in J.items()}

def lam_enclosure():
    J = jets_full('B', arb(0), arb(0), 0)
    lam = J[(0, 0)].real
    return lam

LAM = lam_enclosure()

def Acoef(Jk):
    """Taylor coefficients of A = 2K - lam from those of K (real parts; K real for real (s,xi))."""
    out = {}
    for k, v in Jk.items():
        r = 2 * v.real
        if k == (0, 0):
            r = r - LAM
        out[k] = r
    return out

def Bcoef(Jb):
    return dict(Jb)

# ------------------------------------------------------------------------------------------------------------------
# global sup bounds U_K(a,b), U_B(a,b) of |d^a_s d^b_xi| over the whole domain |s|<=1/2 (all xi)
# ------------------------------------------------------------------------------------------------------------------
def global_U(kind, a, b):
    tot = arb(0)
    half = arb(1) / 2
    for (m, n) in TERMS:
        am, an = abs(m), abs(n)
        dm = max(am - half, arb(0)) if am > 0 else arb(0)
        if kind == 'A':
            ex_ = (-PI * T2 * (dm * dm + an * an)).exp()
            Ebar = TWOPI * T2 * (am + half) + TWOPI * abs(T1) * an
            e2b = PI * T2; exb = TWOPI * an
        else:
            dn = max(an - half, arb(0)) if an > 0 else arb(0)
            ex_ = (-PI * T2 * (dm * dm + dn * dn)).exp()
            Ebar = TWOPI * TAUABS * (abs(m - n) + 1); e2b = TWOPI * TAUABS; exb = TWOPI * abs(n - m)
        tot += C0 * ex_ * bound_poly(a, b, Ebar, e2b, exb)
    tot += tailJ(a, b) if (a + b) <= 6 else 0
    return (factorial(a) * factorial(b) * tot).upper()

def build_kappa6():
    UK = {}; UB = {}
    for a in range(DEG6 + 1):
        for b in range(DEG6 + 1 - a):
            UK[(a, b)] = global_U('A', a, b)
            UB[(a, b)] = global_U('B', a, b)
    UA = {k: 2 * v for k, v in UK.items()}
    UA[(0, 0)] = UA[(0, 0)] + UK[(0, 0)]          # |lam| <= sup|K| as well
    kap = arb(0)
    for a in range(DEG6 + 1):
        b = DEG6 - a
        s = arb(0)
        for i in range(a + 1):
            for j in range(b + 1):
                s += comb(a, i) * comb(b, j) * (UA[(i, j)] * UA[(a - i, b - j)] + UB[(i, j)] * UB[(a - i, b - j)])
        kap += s / (factorial(a) * factorial(b))
    return kap.upper(), UA, UB

# ------------------------------------------------------------------------------------------------------------------
# polynomial helpers (dict (a,b) -> ball)
# ------------------------------------------------------------------------------------------------------------------
def pmul(P, Q, deg):
    R = {}
    for (a, b), x in P.items():
        for (c, d), y in Q.items():
            if a + b + c + d <= deg:
                k = (a + c, b + d)
                R[k] = R[k] + x * y if k in R else x * y
    return R

def sigma_poly(Ac, Bc, deg):
    """real polynomial coefficients of A^2 - B conj(B), degree <= deg."""
    AA = pmul(Ac, Ac, deg)
    Bq = {k: acb(v) if not isinstance(v, acb) else v for k, v in Bc.items()}
    Bb = {k: v.conjugate() for k, v in Bq.items()}
    BB = pmul(Bq, Bb, deg)
    S = {}
    for k in set(AA) | set(BB):
        x = AA.get(k, arb(0))
        y = BB.get(k, acb(0))
        S[k] = (acb(x) - y).real
    return S

# ------------------------------------------------------------------------------------------------------------------
# Part I: near region
# ------------------------------------------------------------------------------------------------------------------
def Rtail(y, K0):
    """sum_{k>=K0} y^k/k!  <=  y^K0/K0! e^y"""
    return y ** K0 / factorial(K0) * y.exp()

def near_certificate(r1, Ndeg):
    hdr(f"[near] r1={r1.mid()} N={Ndeg}")
    JK = jets_full('A', arb(0), arb(0), Ndeg)
    JB = jets_full('B', arb(0), arb(0), Ndeg)
    Ac = Acoef(JK); Bc = Bcoef(JB)
    hdr("[near] coefficients done; A(0)-lam, B(0)-lam enclosures:")
    print("    A00 =", Ac[(0, 0)].str(20, radius=True), "  Re B00 =", Bc[(0, 0)].real.str(20, radius=True), "  lam =", LAM.str(20, radius=True))
    S = sigma_poly(Ac, Bc, 2 * Ndeg)
    # Section 4 proves these coefficients vanish; check that their intervals contain zero.
    worst = 0.0
    for (a, b), v in S.items():
        if a + b <= 3:
            assert arb(0).overlaps(v), f"low order coeff {(a,b)} not containing 0: {v}"
            worst = max(worst, float(ub(v).upper()))
    print("    max |coef| enclosure over degree<=3 coefficients of Sigma: %.3e" % worst)
    for k in [(2, 0), (1, 1), (0, 2)]:
        print(f"    S{k} =", S[k].str(5, radius=True))
    # ---- P4 and its minimum on the sup-norm unit-square boundary
    P4 = {k: v for k, v in S.items() if k[0] + k[1] == 4}
    NSUB = 4096
    p4lo = None
    for edge in (0, 1):
        for i in range(NSUB):
            u = arb(-1) + arb(2 * i + 1) / NSUB                 # midpoint
            u = arb(u.mid(), arb(1) / NSUB)                      # ball of radius 1/NSUB
            if edge == 0:
                s_, x_ = arb(1), u
            else:
                s_, x_ = u, arb(1)
            val = arb(0)
            for (a, b), c in P4.items():
                val += c * s_ ** a * x_ ** b
            lo = val.lower()
            p4lo = lo if p4lo is None else min(p4lo, lo)
    print("    P4 min over sup-norm unit square boundary >=", p4lo.str(8))
    assert p4lo > 0, "P4 not positive"
    # ---- coefficient abs sums for degree 5..2N
    csum = {}
    for (a, b), v in S.items():
        d = a + b
        if d >= 5:
            csum[d] = csum.get(d, arb(0)) + ub(v)
    higher = arb(0)
    for d, c in csum.items():
        higher += c * r1 ** (d - 4)
    print("    sum_{j>=5} c_j r1^{j-4} <=", higher.upper().str(8))
    # ---- Taylor tail eps (beyond degree N) of K, B on |s|,|xi|<=r1, all lattice terms
    K0 = (Ndeg + 2) // 2                                         # ceil((N+1)/2)
    epsK = arb(0); epsB = arb(0)
    for (m, n) in TERMS:
        am = abs(m); an = abs(n)
        f0 = C0 * (-PI * T2 * (m * m + n * n)).exp()
        # kind A
        e1 = ub(acb(-TWOPI * T2 * m, TWOPI * T1 * n)); ex = TWOPI * an; e2 = PI * T2
        y = (e1 + ex) * r1 + e2 * r1 * r1
        epsK += f0 * Rtail(y, K0)
        # kind B
        e1 = TWOPI * TAUABS * abs(m - n); ex = TWOPI * abs(m - n); e2 = TWOPI * TAUABS
        y = (e1 + ex) * r1 + e2 * r1 * r1
        epsB += f0 * Rtail(y, K0)
    # shell terms t>M
    def shell_eps(t):
        t = arb(t)
        alpha = (TWOPI * TAUABS * 2 + 4 * PI) * t + TWOPI * TAUABS
        y = alpha * r1 + TWOPI * TAUABS * r1 * r1
        return 8 * t * C0 * (-PI * T2 * (t - arb(1) / 2) ** 2).exp() * Rtail(y, K0)
    sh = arb(0)
    for t in range(M + 1, M + 31):
        sh += shell_eps(t)
    T = M + 30
    alpha1 = (TWOPI * TAUABS * 2 + 4 * PI) * r1
    rho = (arb(T + 1) / T) ** (K0 + 1) * (-TWOPI * T2 * T).exp() * alpha1.exp()
    assert rho < arb(1) / 2, "near tail ratio"
    sh = (sh + shell_eps(T + 1) / (1 - rho)).upper()
    epsK = epsK + sh; epsB = epsB + sh
    epsA = 2 * epsK
    eps = max(epsA.upper(), epsB.upper())
    print("    Taylor-tail eps (A, B) at r1 <=", epsA.upper().str(5), epsB.upper().str(5))
    Abar = arb(0); Bbar = arb(0)
    for (a, b), v in Ac.items():
        Abar += ub(v) * r1 ** (a + b)
    for (a, b), v in Bc.items():
        Bbar += ub(v) * r1 ** (a + b)
    Etot = (2 * (Abar + Bbar) * eps + 2 * eps * eps).upper()
    margin = p4lo - higher - Etot / r1 ** 4
    print("    Abar, Bbar =", Abar.upper().str(6), Bbar.upper().str(6), " Etot/r1^4 =", (Etot / r1 ** 4).upper().str(5))
    print("    NEAR MARGIN  p4lo - sum_{j>=5} c_j r1^{j-4} - Etot/r1^4 =", margin.str(8))
    ok = bool(margin > 0)
    print("CERT near region: Sigma >= rho^4 * margin > 0 for 0 < rho <= r1 :", ok)
    return ok, p4lo, margin

# ------------------------------------------------------------------------------------------------------------------
# Part II: far region
# ------------------------------------------------------------------------------------------------------------------
def cell_test(s0, xi0, h, KAP6):
    JK = jets_full('A', s0, xi0, ORD)
    JB = jets_full('B', s0, xi0, ORD)
    Ac = Acoef(JK); Bc = Bcoef(JB)
    S = sigma_poly(Ac, Bc, ORD)
    low = S[(0, 0)]
    hp = arb(1)
    for d in range(1, ORD + 1):
        hp = hp * h
        sd = arb(0)
        for a in range(d + 1):
            sd += ub(S[(a, d - a)])
        low = low - sd * hp
    low = low - KAP6 * h ** 6
    return low

def far_certificate(r1, KAP6):
    hmin = arb(1) / (2 ** HMIN_EXP)
    half = arb(1) / 2
    stack = []
    h0 = arb(1) / 8
    for i in range(2):
        for j in range(4):
            stack.append((h0 + 2 * h0 * i, -half + h0 + 2 * h0 * j, h0))
    passed = {}; skipped = 0; ncells = 0
    worst = None
    while stack:
        s0, xi0, h = stack.pop()
        if (s0 + h <= r1) and (abs(xi0) + h <= r1):
            skipped += 1
            continue
        low = cell_test(s0, xi0, h, KAP6)
        ncells += 1
        if low > 0:
            key = float(h.mid())
            passed[key] = passed.get(key, 0) + 1
            lw = low.lower()
            worst = lw if worst is None or lw < worst else worst
        else:
            if h <= hmin:
                print("FAIL at cell", float(s0.mid()), float(xi0.mid()), float(h.mid()), "lower bound", low)
                return False, passed, skipped
            h2 = h / 2
            for ds in (-1, 1):
                for dx in (-1, 1):
                    stack.append((s0 + ds * h2, xi0 + dx * h2, h2))
        if ncells % 200 == 0:
            hdr(f"[far] tested {ncells} cells, stack {len(stack)}, passed {sum(passed.values())}")
    print("CERT far region: every cell passes:", True, " cells tested", ncells, " passed", sum(passed.values()), " skipped (inside near box)", skipped)
    print("    passed cells by half-width h:", dict(sorted(passed.items())))
    print("    smallest cell lower bound over passed cells (absolute):", worst.str(6) if worst is not None else None)
    return True, passed, skipped

if __name__ == "__main__":
    hdr(f"tau = {T1.mid()} + {T2.mid()} i   M={M} prec={PREC} r1={R1.mid()} N={N}  (rational inputs enclosed by balls)")
    print("lam =", LAM.str(25, radius=True))
    KAP6, UA, UB = build_kappa6()
    hdr(f"kappa6 (global sup of 6th derivative sum) <= {KAP6.str(6)}")
    okN, p4lo, marginN = near_certificate(R1, N)
    okF, passed, skipped = far_certificate(R1, KAP6)
    ok = okN and okF
    print("\nSUMMARY tau=%s+%si: near=%s far=%s  =>  CERT Sigma>0 on T^2\\{0}: %s" % (sys.argv[1], sys.argv[2], okN, okF, ok))
    sys.exit(0 if ok else 1)
