"""Exploration (float64, NUMERICAL): GST LLL Bloch symbol Sigma_tau(s,xi)=A^2-|B|^2 and its Taylor/quartic data.
Definitions are those of the sci-gstlattice manuscript: G(x,y)=c exp(-pi t2 (x^2+y^2)+2 pi i t1 x y), c=sqrt(t2)/pi,
K(s,xi)=sum G(m+s,n)e^{2 pi i xi n}, B(s,xi)=sum G(m+s,n-s)e^{2 pi i xi (n-m)}, A=2K-lam, Sigma=A^2-|B|^2.
"""
import numpy as np
from math import factorial

def lattice(M):
    m, n = np.meshgrid(np.arange(-M, M + 1), np.arange(-M, M + 1), indexing='ij')
    return m.ravel().astype(float), n.ravel().astype(float)

def lam_of(tau, M=10):
    t1, t2 = tau.real, tau.imag
    m, n = lattice(M)
    c = np.sqrt(t2) / np.pi
    return (c * np.exp(-np.pi * t2 * (m * m + n * n) + 2j * np.pi * t1 * m * n)).sum().real

def KB(tau, s, xi, M=10):
    """K, B arrays for points (s,xi) (1-D arrays of equal length)."""
    t1, t2 = tau.real, tau.imag
    c = np.sqrt(t2) / np.pi
    m, n = lattice(M)
    s = np.asarray(s, float)[:, None]; xi = np.asarray(xi, float)[:, None]
    GA = np.exp(-np.pi * t2 * ((m + s) ** 2 + n ** 2) + 2j * np.pi * t1 * (m + s) * n + 2j * np.pi * xi * n)
    GB = np.exp(-np.pi * t2 * ((m + s) ** 2 + (n - s) ** 2) + 2j * np.pi * t1 * (m + s) * (n - s) + 2j * np.pi * xi * (n - m))
    return c * GA.sum(1), c * GB.sum(1)

def sigma(tau, s, xi, M=10):
    lam = lam_of(tau, M)
    K, B = KB(tau, s, xi, M)
    A = 2 * K.real - lam
    return A * A - np.abs(B) ** 2

# ---- Taylor coefficients at the origin via moments ----
def _poly_exp(e_s, e_x, e_ss, N):
    """Taylor coefficients (shape (npts,N+1,N+1), index [a,b] of s^a xi^b) of exp(e_s*s+e_x*xi+e_ss*s^2)."""
    npts = len(e_s)
    # build E polynomial, then exp by series sum E^j/j!
    E = np.zeros((npts, N + 1, N + 1), complex)
    E[:, 1, 0] = e_s; E[:, 0, 1] = e_x; E[:, 2, 0] = e_ss
    out = np.zeros_like(E); out[:, 0, 0] = 1
    term = out.copy()
    for j in range(1, N + 1):
        term = _pmul(term, E, N) / j
        out += term
    return out

def _pmul(P, Q, N):
    R = np.zeros_like(P)
    for a in range(N + 1):
        for b in range(N + 1 - a):
            p = P[:, a, b]
            if not np.any(p):
                continue
            R[:, a:, b:] += p[:, None, None] * Q[:, :N + 1 - a, :N + 1 - b]
    return R

def taylor_AB(tau, N=8, M=10):
    t1, t2 = tau.real, tau.imag
    tb = np.conj(tau)
    c = np.sqrt(t2) / np.pi
    m, n = lattice(M)
    w = c * np.exp(-np.pi * t2 * (m * m + n * n) + 2j * np.pi * t1 * m * n)
    lam = w.sum().real
    # K
    eK = _poly_exp(-2 * np.pi * t2 * m + 2j * np.pi * t1 * n, 2j * np.pi * n, -np.pi * t2 * np.ones_like(m), N)
    Kc = np.einsum('p,pab->ab', w, eK)
    Ac = 2 * Kc
    Ac[0, 0] -= lam
    eB = _poly_exp(-2j * np.pi * (m - n) * tb, -2j * np.pi * (m - n), -2j * np.pi * tb * np.ones_like(m), N)
    Bc = np.einsum('p,pab->ab', w, eB)
    return lam, Ac, Bc

def sigma_taylor(tau, N=8, M=10):
    lam, Ac, Bc = taylor_AB(tau, N, M)
    P = np.zeros((N + 1, N + 1), complex)
    def mul(U, V):
        R = np.zeros_like(U)
        for a in range(N + 1):
            for b in range(N + 1 - a):
                R[a:, b:] += U[a, b] * V[:N + 1 - a, :N + 1 - b]
        return R
    S = mul(Ac, Ac) - mul(Bc, np.conj(Bc))
    return lam, Ac, Bc, S

def homog(S, d):
    """coefficients c_k of s^(d-k) xi^k in degree d part"""
    return np.array([S[d - k, k] for k in range(d + 1)])

def P_theta(S, d, th):
    co = homog(S, d)
    s = np.cos(th); x = np.sin(th)
    return sum(co[k].real * s ** (d - k) * x ** k for k in range(d + 1))

def quartic_min(tau, nth=2000, N=6):
    lam, Ac, Bc, S = sigma_taylor(tau, N)
    th = np.linspace(0, np.pi, nth, endpoint=False)
    p4 = P_theta(S, 4, th)
    i = np.argmin(p4)
    return p4.min(), th[i], S

if __name__ == "__main__":
    for tau in [complex(-0.4, 0.95), complex(0.4, 1.0), complex(-0.5, np.sqrt(3) / 2), complex(0.5, 1.04), complex(0.5, 1.052), 1j]:
        mn, th, S = quartic_min(tau)
        print(tau, "min P4 on unit circle", mn, "at angle", th, "low-order |coef|", [abs(S[a, b]) for a in range(4) for b in range(4 - a)][:3])
