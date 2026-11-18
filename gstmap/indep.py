"""Sample the interiors of the Arb certificate boxes with a separate double-precision calculation.

The evaluator implements its own Gaussian sums with NumPy and a different truncation, using the S transformation
Sigma_tau(s,xi) = Sigma_{-1/tau}(xi,-s) from Proposition 1. It supplies a numerical cross-check.
"""
import math
import numpy as np


def _lat(t1, t2, s, xi, M):
    """K (real), B (complex), lam for arrays s, xi (same shape)."""
    c = math.sqrt(t2) / math.pi
    K = np.zeros_like(s)
    B = np.zeros(s.shape, dtype=complex)
    lam = 0.0
    rng = range(-M, M + 1)
    for m in rng:
        for n in rng:
            lam += c * math.exp(-math.pi * t2 * (m * m + n * n)) * math.cos(2 * math.pi * t1 * m * n)
            x = m + s
            K += c * np.exp(-math.pi * t2 * (x * x + n * n)) * np.cos(2 * math.pi * (t1 * x * n + xi * n))
            y = n - s
            B += c * np.exp(-math.pi * t2 * (x * x + y * y) + 2j * math.pi * (t1 * x * y + xi * (n - m)))
    return K, B, lam


def sigma_direct(t1, t2, s, xi, M=8):
    s = np.asarray(s, float)
    xi = np.asarray(xi, float)
    K, B, lam = _lat(t1, t2, s, xi, M)
    A = 2 * K - lam
    return A * A - (B.real ** 2 + B.imag ** 2)


def sigma_dual(t1, t2, s, xi, M=8):
    """Sigma_tau(s,xi) evaluated through tau' = -1/tau at (xi,-s)."""
    n2 = t1 * t1 + t2 * t2
    t1p, t2p = -t1 / n2, t2 / n2
    return sigma_direct(t1p, t2p, np.asarray(xi, float), -np.asarray(s, float), M)


def torus_min_ratio(t1, t2, nrad=60, nth=72, ngrid=121, rmin=0.01, rmax=0.5, use_dual=True):
    """min over Bloch torus of Sigma/rho^4 (|q|>=rmin) via polar grid on |q|<=rmax and Cartesian grid outside."""
    f = sigma_dual if use_dual else sigma_direct
    best = (1e99, None)
    rs = rmin * (rmax / rmin) ** (np.arange(nrad) / (nrad - 1))
    th = np.linspace(0, 2 * math.pi, nth, endpoint=False)
    R, TH = np.meshgrid(rs, th, indexing="ij")
    S, X = R * np.cos(TH), R * np.sin(TH)
    val = f(t1, t2, S, X) / R ** 4
    i = np.unravel_index(np.argmin(val), val.shape)
    best = (float(val[i]), (float(S[i]), float(X[i])))
    g = np.linspace(-0.5, 0.5, ngrid)
    S2, X2 = np.meshgrid(g, g, indexing="ij")
    R2 = np.hypot(S2, X2)
    mask = R2 >= rmax
    v2 = f(t1, t2, S2, X2)
    # outside the polar disc: report Sigma itself (positive) normalised by rho^4 for comparability
    r4 = np.where(mask, R2, 1.0) ** 4
    ratio = np.where(mask, v2 / r4, 1e99)
    j = np.unravel_index(np.argmin(ratio), ratio.shape)
    if ratio[j] < best[0]:
        best = (float(ratio[j]), (float(S2[j]), float(X2[j])))
    return best


def box_replay(c1, c2, w, kind, n_rand=12, seed=1, q0=None, **kw):
    """Replay a certificate box at corners, centre and random interior points.
    kind 'S': torus_min_ratio must be > 0 everywhere;   kind 'U': Sigma(q0) < 0 (q0 given) at all points."""
    rng = np.random.default_rng(seed)
    pts = [(c1 + a * w, c2 + b * w) for a in (-1, 0, 1) for b in (-1, 0, 1)]
    pts += [(c1 + w * (2 * rng.random() - 1), c2 + w * (2 * rng.random() - 1)) for _ in range(n_rand)]
    out = []
    for (t1, t2) in pts:
        if kind == "S":
            out.append(torus_min_ratio(t1, t2, **kw)[0])
        else:
            out.append(float(sigma_dual(t1, t2, np.array([q0[0]]), np.array([q0[1]]))[0]))
    return out
