"""NUMERICAL (mpmath 40 digits): decompose the exact quartic coefficient of Sigma at k=0.
   D = A - |B| = c4 + O(k^6),  c4 = A4 - Re B4 - (Im B2)^2/(2 lam),  P4 = 2 lam c4   (A2 = Re B2, B1 = B3 = 0).
Reports: P4 vs 2*pi*lam*t2^2*Q_mirror(theta), the difference, Im B2(q) form, and a least-squares test of
   P4 - 2 pi lam t2^2 Q  ~  beta * (Im B2)^2.
Run: uv run --no-project --with numpy --with scipy --with mpmath python quartic_decomp.py
"""
import mpmath as mp, numpy as np, sys
from taylor_quartic import run, poly
from pocket_scan import Q as Qshear

def qform(co, order, t, t1, t2):
    qx, qy = np.cos(t), np.sin(t)
    s = mp.mpf(qx); xi = mp.mpf(qy * t2 - qx * t1)
    return poly(co, order, (s, xi))

for tau in [-0.4 + 0.95j, 0.4 + 1.0j, 0.5 + 1.04j, -0.45 + 0.9j, 0.3 + 1.0j]:
    tau = complex(tau); t1, t2 = tau.real, tau.imag
    lam, A, B, S = run(tau)
    lam = float(lam)
    print(f"\ntau={tau} lam={lam:.12f} 2*pi*lam*t2^2={2*np.pi*lam*t2**2:.10f}")
    ths = np.linspace(0, np.pi, 25)[:-1]
    rows = []; P4s = []; Qs = []
    for t in ths:
        p4 = float(qform(S, 4, t, t1, t2).real)
        a2 = float(qform(A, 2, t, t1, t2).real); rb2 = float(qform(B, 2, t, t1, t2).real); ib2 = float(qform(B, 2, t, t1, t2).imag)
        a4 = float(qform(A, 4, t, t1, t2).real); rb4 = float(qform(B, 4, t, t1, t2).real)
        c4 = a4 - rb4 - ib2 ** 2 / (2 * lam)
        qm = Qshear(tau, np.pi - t)
        P4s.append(p4); Qs.append(qm)
        rows.append((t, p4, 2 * lam * c4, a2 - rb2, ib2, qm))
    print("   theta      P4         2*lam*c4    A2-ReB2     ImB2(q)      Qmirror   P4/Qm")
    for r in rows[::3]:
        print(f"  {r[0]:7.4f} {r[1]:11.7f} {r[2]:11.7f} {r[3]:10.2e} {r[4]:11.7f} {r[5]:10.6f} {r[1]/r[5] if abs(r[5])>1e-9 else float('nan'):10.6f}")
    P4s = np.array(P4s); Qs = np.array(Qs)
    ib = np.array([r[4] for r in rows])
    # fit P4 = a*Q + b*ImB2^2 + c
    X = np.vstack([Qs, ib ** 2, np.ones_like(Qs)]).T
    coef, res, *_ = np.linalg.lstsq(X, P4s, rcond=None)
    print(f"  lstsq P4 ~ a Q + b (ImB2)^2 + c: a={coef[0]:.8f} b={coef[1]:.6f} c={coef[2]:.3e}  max resid={np.abs(X@coef-P4s).max():.3e}   (2 pi lam t2^2 = {2*np.pi*lam*t2**2:.8f})")
    X2 = np.vstack([Qs, np.ones_like(Qs)]).T
    c2, *_ = np.linalg.lstsq(X2, P4s, rcond=None)
    print(f"  lstsq P4 ~ a Q + c: a={c2[0]:.8f} c={c2[1]:.3e} max resid={np.abs(X2@c2-P4s).max():.3e}")
    sys.stdout.flush()
