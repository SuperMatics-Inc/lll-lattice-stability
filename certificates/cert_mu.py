# Arb certification of the GST hexagonal-lattice facts that the paper only reads off Fig. 1.
# All "CERT" lines are rigorous ball-arithmetic statements (python-flint / Arb), prec=320 bits.
import sys, time
from flint import arb, ctx
ctx.prec = 320
PI = arb.pi(); SQ3 = arb(3).sqrt()
Q = (-PI*SQ3/2).exp()                     # q = exp(-pi*sqrt3/2)
C0 = (SQ3/2).sqrt()/PI                    # 1/(gamma sqrt(pi)), gamma^2 = 2 pi/sqrt3
K = 14

def _cos(w, x, r):                        # cos(w x + r pi/2) on ball x
    return (w*x + PI*r/2).cos()

def ell_der(r, x):
    """enclosure of ell^{(r)}(xi), xi in ball x; ell = 1 + 2 sum_{k>=1} q^{k^2} cos(2 pi k xi)"""
    tot = arb(1) if r == 0 else arb(0)
    for k in range(1, K+1):
        w = 2*PI*k
        tot += 2*Q**(k*k)*w**r*_cos(w, x, r)
    tail = 4*Q**((K+1)**2)*(2*PI*(K+1))**r
    return tot + arb(0, tail.upper())

def h_der(r, x):
    """h = 2 sum_{k>=0} q^{(2k+1)^2} cos(2 (2k+1) pi xi)"""
    tot = arb(0)
    for k in range(0, K+1):
        n = 2*k+1; w = 2*PI*n
        tot += 2*Q**(n*n)*w**r*_cos(w, x, r)
    n = 2*K+3
    tail = 4*Q**(n*n)*(2*PI*n)**r
    return tot + arb(0, tail.upper())

from math import comb
L0 = ell_der(0, arb(0)); H0 = h_der(0, arb(0))
LAM = C0*(L0**2 - 2*H0**2)

def jets_ab(x, nmax):
    """derivatives a^{(j)}, b^{(j)}, j=0..nmax on ball x"""
    l = [ell_der(j, x) for j in range(nmax+1)]
    h = [h_der(j, x) for j in range(nmax+1)]
    a = []; b = []
    for j in range(nmax+1):
        aj = 2*C0*(l[j]*L0 - 2*h[j]*H0)
        if j == 0: aj = aj - LAM
        pj = sum(comb(j, i)*(l[i]*l[j-i] - 2*h[i]*h[j-i]) for i in range(j+1))
        a.append(aj); b.append(C0*pj)
    return a, b

NT = 6
def jets_s(x):
    """enclosures of s^{(k)}(xi), k=0..3, xi in ball x, via Taylor form about the midpoint
    (point jets to order NT-1, order-NT jet over the whole ball) to avoid the a^2-b^2 cancellation."""
    m = arb(x.mid()); r = x.rad()
    d = arb(0, r.upper())
    ap, bp = jets_ab(m, NT-1)
    Sp = [sum(comb(j, i)*(ap[i]*ap[j-i] - bp[i]*bp[j-i]) for i in range(j+1)) for j in range(NT)]
    aI, bI = jets_ab(x, NT)
    SN = sum(comb(NT, i)*(aI[i]*aI[NT-i] - bI[i]*bI[NT-i]) for i in range(NT+1))
    from math import factorial
    out = []
    dpow = [arb(1)]
    for _ in range(NT): dpow.append(dpow[-1]*d)
    for k in range(4):
        e = SN*dpow[NT-k]/factorial(NT-k)
        for j in range(NT-k):
            e += Sp[k+j]*dpow[j]/factorial(j)
        out.append(e)
    a0, b0 = jets_ab(x, 0)
    return out, a0, b0

def mu_derivs(s):
    s0, s1, s2, s3 = s
    r = s0.sqrt()
    P = 2*s0*s2 - s1*s1
    mu2 = P/(4*s0*r)
    mu3 = s3/(2*r) - 3*s1*P/(8*s0*s0*r)
    return P, mu2, mu3

def lo(x): return x.lower()
def up(x): return x.upper()

# Every interval hull below must contain the complete endpoint balls, not
# merely their midpoints. Counts are printed only after a successful run.
CHECKS = {'region1': 0, 'region2': 0, 'cover': 0, 'root': 0,
          'cover_exact_6_5': 0, 'near_exact_1707_1000': 0}

# ---------------- Region II: xi in [W0, 1/2] ----------------
W0 = arb(1)/20
def region2():
    todo = []
    n0 = 45
    for i in range(n0):
        todo.append((W0 + (arb(1)/2 - W0)*i/n0, W0 + (arb(1)/2 - W0)*(i+1)/n0, 0))
    out = []
    while todo:
        lo_, hi_, dep = todo.pop(0)
        mid = (lo_+hi_)/2; rad = (hi_-lo_)/2
        x = arb(mid.mid(), rad.upper())
        assert x.contains(lo_) and x.contains(hi_), 'region-II endpoint containment'
        CHECKS['region2'] += 1
        s, a, b = jets_s(x)
        cls = None
        if lo(s[0]) > 0:
            P, mu2, mu3 = mu_derivs(s)
            if lo(P) > 0: cls = 'A'
            elif up(P) < 0: cls = 'B'
            elif ((lo(s[3]) > 0) or (up(s[3]) < 0)) and mu3.abs_lower() > 0:
                cls = 'C'      # P' = 2 s s''' has a strict sign, so P (hence mu'') has at most one zero here
        if cls is None:
            if dep > 40: raise SystemExit('depth exceeded at %s' % x)
            todo.insert(0, (mid, hi_, dep+1)); todo.insert(0, (lo_, mid, dep+1)); continue
        out.append((lo_, hi_, cls, mu2, mu3, s, a, b))
    out.sort(key=lambda t: float(t[0].mid()))
    return out

# ---------------- Region I: xi in [0, W0] ----------------
def region1(n=10):
    res = []
    for i in range(n):
        lo_ = W0*i/n; hi_ = W0*(i+1)/n
        x = arb(((lo_+hi_)/2).mid(), ((hi_-lo_)/2).upper())
        x0 = arb((hi_/2).mid(), (hi_/2).upper())      # hull of [0,hi]
        assert x.contains(lo_) and x.contains(hi_), 'region-I endpoint containment'
        assert x0.contains(arb(0)) and x0.contains(hi_), 'region-I auxiliary containment'
        CHECKS['region1'] += 1
        u = []; v = []
        for j in range(3):
            Wj = arb(1)/((j+1)*(j+2))
            u.append(-Wj*ell_der(j+2, x0)); v.append(-Wj*h_der(j+2, x0))
        a, b = jets_ab(x, 2)
        U = [C0*(2*v[0]**2 - u[0]**2),
             C0*(4*v[0]*v[1] - 2*u[0]*u[1]),
             C0*(4*v[1]**2 + 4*v[0]*v[2] - 2*u[1]**2 - 2*u[0]*u[2])]
        V = [a[j]+b[j] for j in range(3)]
        R = [U[0]*V[0], U[1]*V[0]+U[0]*V[1], U[2]*V[0]+2*U[1]*V[1]+U[0]*V[2]]
        if not (lo(U[0]) > 0 and lo(V[0]) > 0): raise SystemExit('U or V not positive in region I')
        m = R[0].sqrt()
        m1 = R[1]/(2*m); m2 = R[2]/(2*m) - R[1]**2/(4*m**3)
        mu2 = 2*m + 4*x*m1 + x*x*m2
        res.append((lo_, hi_, mu2, U[0], V[0], m))
    return res

def cover(lo_, hi_, target, dep=0):
    """adaptive bisection: certify |mu''|+|mu'''| >= target on [lo_,hi_] (arb endpoints)"""
    mid = (lo_+hi_)/2; rad = (hi_-lo_)/2
    x = arb(mid.mid(), rad.upper())
    assert x.contains(lo_) and x.contains(hi_), 'cover endpoint containment'
    CHECKS['cover'] += 1
    s, a, b = jets_s(x)
    if lo(s[0]) > 0:
        P, m2, m3 = mu_derivs(s)
        if (m2.abs_lower() + m3.abs_lower()) >= target:
            assert (m2.abs_lower() + m3.abs_lower()) >= arb(6)/5, 'exact 6/5 threshold'
            CHECKS['cover_exact_6_5'] += 1
            return True, 1
    if dep > 25: return False, 1
    ok1, n1 = cover(lo_, mid, target, dep+1)
    if not ok1: return False, n1
    ok2, n2 = cover(mid, hi_, target, dep+1)
    return ok2, n1+n2

if __name__ == '__main__':
    t0 = time.time()
    print('q =', Q.str(30), ' c0 =', C0.str(30), ' lambda =', LAM.str(30))
    print('CERT q<1/15:', Q < arb(1)/15, ' (q<0.1 needed by paper):', Q < arb(1)/10)
    # sign lemma l''(0)^2 < 2 h''(0)^2
    l2 = ell_der(2, arb(0)); h2 = h_der(2, arb(0))
    print("l''(0) =", l2.str(25), " h''(0) =", h2.str(25), " ratio l''/h'' =", (l2/h2).str(20), " sqrt2 =", arb(2).sqrt().str(10))
    print("CERT 2 h''(0)^2 - l''(0)^2 > 0:", (2*h2*h2 - l2*l2) > 0, " value", (2*h2*h2-l2*l2).str(15))
    C = (C0*(2*(h2/2)**2 - (l2/2)**2)*(2*LAM)).sqrt()   # m(0)^2 = U(0)*(a+b)(0), a(0)=b(0)=lam
    print("C = m(0) = lim mu/xi^2 =", C.str(20), " mu''(0) = 2C =", (2*C).str(20))
    r1 = region1(10)
    assert all(t[2] > arb(1707)/1000 for t in r1), 'exact 1707/1000 threshold'
    CHECKS['near_exact_1707_1000'] = len(r1)
    print('Region I [0,%s] : mu\'\' lower bound over 10 subintervals: %s' % (W0.str(3), min(float(lo(t[2])) for t in r1)))
    print('  CERT mu\'\'>0 on [0,1/20]:', all(lo(t[2]) > 0 for t in r1), ' min U lb:', min(float(lo(t[3])) for t in r1), ' min m lb:', min(float(lo(t[5])) for t in r1))
    r2 = region2()
    print('Region II [1/20,1/2]: %d subintervals' % len(r2))
    seq = ''.join(t[2] for t in r2)
    import re
    print('class sequence compressed:', re.sub(r'(.)\1+', lambda m: m.group(1)+'x'+str(len(m.group(0))), seq))
    # structure: A* C+ B*
    assert re.fullmatch(r'A+C+B+', seq), 'unexpected class sequence'
    cidx = [i for i, t in enumerate(r2) if t[2] == 'C']
    assert all(lo(-r2[i][5][3]) > 0 for i in cidx) or all(lo(r2[i][5][3]) > 0 for i in cidx)
    print('CERT exactly one zero of mu\'\' in (1/20,1/2): A* C+ B* with s\'\'\' of one strict sign on the C-block; mu\'\'>0 on [0,1/20]; mu\'\'(1/2)<0 (class B at end)')
    zl, zr = r2[cidx[0]][0], r2[cidx[-1]][1]
    print('  zero xi1 in', arb((zl+zr)/2).mid().str(8), '+-', ((zr-zl)/2).str(3))
    # certified lower bound c for |mu''|+|mu'''| on [1/20,1/2]; region I has mu'' >= 1.7 there
    c_target = arb(6)/5
    ok, ncell = cover(W0, arb(1)/2, c_target)
    c_lb1 = min(float(lo(t[2])) for t in r1)
    print('CERT |mu\'\'|+|mu\'\'\'| >= %s on [1/20,1/2] (%s, %d cells) and mu\'\' >= %.3f on [0,1/20]; hence on [0,1] by mu(1-xi)=mu(xi) and evenness' % (c_target.str(3), ok, ncell, c_lb1))
    assert ok and all(t[2] >= c_target for t in r1)
    # tight xi1 by bisection with rigorous point evaluation of P
    lo_, hi_ = zl, zr
    for _ in range(90):
        mid = (lo_+hi_)/2
        xm = arb(mid.mid())
        s, a, b = jets_s(xm)
        P, mu2, mu3 = mu_derivs(s)
        if lo(P) > 0: lo_ = mid.mid()
        elif up(P) < 0: hi_ = mid.mid()
        else: break
    x1 = arb(((lo_+hi_)/2).mid(), ((hi_-lo_)/2).upper())
    assert x1.contains(arb(lo_)) and x1.contains(arb(hi_)), 'root endpoint containment'
    CHECKS['root'] += 1
    s, a, b = jets_s(x1)
    P, mu2, mu3 = mu_derivs(s)
    print('CERT xi1 =', x1.str(25), '  mu\'\'\'(xi1) =', mu3.str(15), ' mu\'\'(xi1) enclosure', mu2.str(5))
    # mu'', mu''' at 1/2 and mu''(0)
    xh = arb(1)/2
    s, a, b = jets_s(xh); P, mu2, mu3 = mu_derivs(s)
    print('mu(1/2)=', s[0].sqrt().str(15), " mu''(1/2) =", mu2.str(15), " mu'''(1/2) =", mu3.str(10))
    # mu'''  zeros  (not needed) : report min |mu'''| over region II where mu'' small? just print sample
    # a>b>0 certification: s>0 verified on region II, U>0 on region I; also a>0,b>0:
    amin = min(float(lo(t[6][0])) for t in r2); bmin = min(float(lo(t[7][0])) for t in r2)
    print('CERT a>0 and b>0 on [1/20,1/2]: min lb a = %.6f, min lb b = %.6f' % (amin, bmin))
    print('done in %.1fs' % (time.time()-t0))
    # sanity (not a proof): central differences of mu from point evaluations vs the analytic jets
    h = arb(1)/10**8
    for xv in (arb(3)/10, arb(1)/5, arb(9)/20):
        f = [jets_s(xv + k*h)[0][0].sqrt() for k in (-1, 0, 1)]
        mu2fd = (f[0] - 2*f[1] + f[2])/(h*h)
        s, a, b = jets_s(xv); P, mu2, mu3 = mu_derivs(s)
        print('check xi=%s: mu\'\' analytic %s  central-diff %s' % (xv.mid().str(3), mu2.str(10), mu2fd.str(10)))
    print('CERT containment and exact-threshold self-checks passed:', CHECKS)
