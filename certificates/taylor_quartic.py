"""NUMERICAL (mpmath, 40 digits, NOT interval-certified): exact Taylor coefficients at k=0 of the LLL symbol
Sigma(s,xi)=A^2-|B|^2 from the closed Gaussian lattice sums
   A(s,xi) = 2*Khat1 - lam,  Khat1 = c sum_{P,d} exp(-pi t2[(P+s)^2+d^2] - 2 pi i d (xi + t1 (P+s)))   (real)
   B(s,xi) = c sum_{m,n}    exp(-pi t2[(m+s)^2+(n-s)^2] + 2 pi i t1 (m+s)(n-s) + 2 pi i xi (n-m))
Checks: orders 0..3 of Sigma vanish; A2 = Re B2 (quadratic term of A-|B| vanishes);
        quartic form P4 and the shear form Q_tau(theta); ratio P4/(|q|^4 Q).
Run:  uv run --no-project --with numpy --with scipy --with mpmath python taylor_quartic.py
"""
import mpmath as mp, numpy as np, sys
from pocket_scan import Q as Qshear
mp.mp.dps = 40
NORD = 8

def exp_coeffs(es, ess, exi, n_s, n_x):
    """Taylor coeffs of exp(es*s + ess*s^2) (list in s) and exp(exi*xi) (list in xi)"""
    cs = []
    for m in range(n_s + 1):
        v = mp.mpc(0)
        for j in range(m // 2 + 1):
            v += ess ** j * es ** (m - 2 * j) / (mp.factorial(j) * mp.factorial(m - 2 * j))
        cs.append(v)
    cx = [exi ** n / mp.factorial(n) for n in range(n_x + 1)]
    return cs, cx

def series(tau, which, M):
    t1, t2 = mp.mpf(tau.real), mp.mpf(tau.imag)
    c = mp.sqrt(t2) / mp.pi
    pi = mp.pi
    co = {(m, n): mp.mpc(0) for m in range(NORD + 1) for n in range(NORD + 1)}
    for P in range(-M, M + 1):
        for d in range(-M, M + 1):
            if which == 'A':
                # exponent in s: -pi t2 (P+s)^2 - 2 pi i d t1 (P+s) ; xi: -2 pi i d xi
                e0 = -pi * t2 * (P * P + d * d) - 2j * pi * d * t1 * P
                es = -2 * pi * t2 * P - 2j * pi * d * t1
                ess = -pi * t2
                exi = -2j * pi * d
            else:
                m_, n_ = P, d  # (m,n)
                e0 = -pi * t2 * (m_ * m_ + n_ * n_) + 2j * pi * t1 * m_ * n_
                es = -2 * pi * t2 * (m_ - n_) + 2j * pi * t1 * (n_ - m_)
                ess = -2 * pi * t2 - 2j * pi * t1
                exi = 2j * pi * (n_ - m_)
            pref = c * mp.exp(e0)
            cs, cx = exp_coeffs(es, ess, exi, NORD, NORD)
            for m in range(NORD + 1):
                for n in range(NORD + 1 - m):
                    co[(m, n)] += pref * cs[m] * cx[n]
    return co

def pmul(a, b):
    out = {}
    for (i, j), u in a.items():
        for (k, l), v in b.items():
            if i + k + j + l <= NORD:
                out[(i + k, j + l)] = out.get((i + k, j + l), 0) + u * v
    return out

def run(tau, M=9):
    A = series(tau, 'A', M); B = series(tau, 'B', M)
    lam = A[(0, 0)].real
    A = {k: 2 * v for k, v in A.items()}; A[(0, 0)] -= lam
    Bc = {k: mp.conj(v) for k, v in B.items()}  # conj of B(s,xi) for real (s,xi): conj coefficients
    S = pmul(A, A)
    BB = pmul(B, Bc)
    S = {k: (S.get(k, 0) - BB.get(k, 0)) for k in set(S) | set(BB)}
    return lam, A, B, S

def poly(co, order, q):
    return sum(co.get((m, order - m), 0) * q[0] ** m * q[1] ** (order - m) for m in range(order + 1))

if __name__ == "__main__":
    for tau in [complex(-0.5, np.sqrt(3) / 2), -0.4 + 0.95j, 0.4 + 1.0j, 1j, 0.5 + 1.04j]:
        tau = complex(tau)
        lam, A, B, S = run(tau)
        print(f"\ntau={tau}  lam={mp.nstr(lam, 15)}")
        for o in range(0, 4):
            mx = max(abs(S.get((m, o - m), 0)) for m in range(o + 1))
            print(f"  max|Sigma coeff order {o}| = {mp.nstr(mx, 3)}")
        for (m, n) in [(2, 0), (1, 1), (0, 2)]:
            print(f"  (A2 - Re B2)[{m},{n}] = {mp.nstr(A[(m, n)].real - B[(m, n)].real, 3)};  |A1|={mp.nstr(abs(A[(1,0)]),2)} |B1|={mp.nstr(abs(B[(1,0)]),2)}")
        t1, t2 = tau.real, tau.imag
        th = np.linspace(0, np.pi, 7)[:-1]
        print("  theta   P4(q)/|q|^4   Q_tau(theta) and its mirror Q_tau(pi-theta), ratios")
        for t in th:
            qx, qy = np.cos(t), np.sin(t)
            s = qx; xi = qy * t2 - qx * t1  # q=(s,(xi+s*t1)/t2)
            p4 = poly(S, 4, (mp.mpf(s), mp.mpf(xi))).real
            qq = Qshear(tau, t); qm = Qshear(tau, np.pi - t)
            print(f"  {t:6.3f}   {mp.nstr(p4, 12):>16}   Q={qq:12.8f} ratio={mp.nstr(p4 / qq, 10) if qq != 0 else 'nan'}   Qmirror={qm:12.8f} ratio={mp.nstr(p4 / qm, 10) if qm != 0 else 'nan'}")
        sys.stdout.flush()
