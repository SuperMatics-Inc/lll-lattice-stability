"""Sample the LLL linearised symbol in double precision for general tau.
Run:  uv run --no-project --with numpy --with scipy python pocket_scan.py
Symbol: Sigma(s,xi)=A^2-|B|^2, A=2*K1hat^s(xi)-lam (real), B=K2hat^s(-xi).  Stable iff Sigma>=0 on the Bloch torus.
Kernels (Poisson-resummed Gaussian sums, verified against quadrature in the working notes):
  K1^s(d)=c*sum_u exp(i*pi*[(tau/2)(u-s)^2-(conj(tau)/2)(u-2d-s)^2]),  c=sqrt(tau2)/pi
  K2^s(j)=c*exp(-i*pi*conj(tau)*(j-2s)^2/2)*sum_k exp(2*pi*i*tau*(k-j/2)^2)
Long-wave comparison form (physical covolume-1 lattice, shear second variation of the Gaussian lattice sum):
  Q_tau(th)=2*pi*sum_w exp(-pi|w|^2) (khat.w)^2 [2*pi*(khat_perp.w)^2-1],  khat=(cos th,sin th).
Physical dimensionless wavevector: q=(s,(xi+s*tau1)/tau2); lattice frame of Q is the mirror image
(tau -> -conj(tau)) of the Bloch frame (mirror leaves stability unchanged).
"""
import numpy as np, time
from scipy.optimize import brentq

def kern(tau, s):
    t2 = tau.imag; c = np.sqrt(t2) / np.pi
    M = int(np.ceil(np.sqrt(85 / (np.pi * t2)))) + 3
    idx = np.arange(-M, M + 1)
    K1 = np.zeros(len(idx), complex); K2 = np.zeros(len(idx), complex)
    kk = np.arange(-2 * M - 6, 2 * M + 7)
    for i, d in enumerate(idx):
        u = np.round(d + s) + np.arange(-M - 2, M + 3)
        K1[i] = c * np.exp(1j * np.pi * (0.5 * tau * (u - s) ** 2 - 0.5 * np.conj(tau) * (u - 2 * d - s) ** 2)).sum()
        K2[i] = c * np.exp(-1j * np.pi * np.conj(tau) * (d - 2 * s) ** 2 / 2) * np.exp(2j * np.pi * tau * (kk - d / 2) ** 2).sum()
    return idx, K1, K2

def ft(idx, K, xi):
    xi = np.atleast_1d(xi)
    return (K[None, :] * np.exp(-2j * np.pi * np.outer(xi, idx))).sum(1)

def sigma(tau, s, xi):
    idx, K1, K2 = kern(tau, s)
    lam = ft(idx, kern(tau, 0.0)[1], 0.0)[0].real
    A = 2 * ft(idx, K1, xi) - lam
    B = ft(idx, K2, -np.atleast_1d(xi))
    return (A ** 2 - np.abs(B) ** 2).real

def Q(tau, th, N=7):
    t1, t2 = tau.real, tau.imag
    m, n = np.meshgrid(np.arange(-N, N + 1), np.arange(-N, N + 1))
    x = ((m + n * t1) / np.sqrt(t2)).ravel(); y = (n * t2 / np.sqrt(t2)).ravel()
    c, s = np.cos(th), np.sin(th)
    ko = c * x + s * y; et = -s * x + c * y
    w = np.exp(-np.pi * (x * x + y * y))
    return 2 * np.pi * np.sum(w * ko ** 2 * (2 * np.pi * et ** 2 - 1))

def minQ(tau, n=361):
    return min(Q(tau, t) for t in np.linspace(0, np.pi, n))

def S4(tau, th, kk=0.03):
    qx, qy = kk * np.cos(th), kk * np.sin(th)
    s = qx; xi = qy * tau.imag - s * tau.real
    if s < 0: s, xi = -s, -xi
    return sigma(tau, s, xi)[0] / (qx * qx + qy * qy) ** 2

def polar_min(tau, kk=0.03, n=91):
    return min(S4(tau, th, kk) for th in np.linspace(0, np.pi, n)[:-1])

def finite_min(tau, ns=32, nx=64, qcut=0.1):
    mn = 1e9
    for i in range(ns):
        s = i / ns; s_ = s - 1 if s > 0.5 else s
        xi = np.arange(nx) / nx; xs = np.where(xi > 0.5, xi - 1, xi)
        v = sigma(tau, s, xi)
        qx = s_; qy = (xs + s_ * tau.real) / tau.imag
        ok = (qx * qx + qy * qy) >= qcut ** 2
        if ok.any(): mn = min(mn, v[ok].min())
    return mn

if __name__ == '__main__':
    t0 = time.time()
    hexa = complex(-0.5, np.sqrt(3) / 2)
    v = np.array([S4(hexa, th, 0.03) for th in np.linspace(0, np.pi, 37)[:-1]])
    print('hex Sigma/|q|^4 over angle (|q|=.03): min %.6f max %.6f  (isotropy check)' % (v.min(), v.max()))
    print('hex Sigma/|k|^4 with k=(2pi/gamma)q, gamma^2=pi/tau2: %.6f' % (v.mean() / (2 * np.pi) ** 4 * (np.pi / hexa.imag) ** 2))
    print('--- pocket map over fundamental domain (S stable both; L long-wave unstable only; F finite-k only; B both;')
    print('    lowercase = sign(Q) disagrees with direct small-k sign of Sigma) ---')
    agree = total = 0
    for t1 in np.arange(0.0, 0.5001, 0.05):
        line = ''
        for t2 in np.arange(np.sqrt(1 - t1 * t1) + 1e-9, 1.61, 0.04):
            tau = complex(t1, t2)
            q = minQ(tau, 181); p = polar_min(tau); f = finite_min(tau)
            lw = (q > 0 and p > 0); fin = f >= -1e-12
            c = 'S' if (lw and fin) else ('L' if (not lw and fin) else ('F' if (lw and not fin) else 'B'))
            total += 1
            if (q > 0) == (p > 0): agree += 1
            else: c = c.lower()
            line += c
        print('tau1=%.2f tau2 from %.3f step .04: %s  (%.0fs)' % (t1, np.sqrt(1 - t1 * t1), line, time.time() - t0), flush=True)
    print('sign(minQ)==sign(min small-k Sigma/|q|^4) at %d of %d grid points' % (agree, total))
    print('--- long-wave boundary tau2*(tau1) (zero of min_theta Q) and full-torus 64x128 min Sigma (|q|>=.05) on both sides ---')
    for t1, lo, hi in [(0.5, 1.0, 1.2), (0.45, 0.95, 1.2), (0.4, 0.95, 1.2), (0.35, 0.94, 1.1), (0.3, 0.954, 1.0)]:
        r = brentq(lambda t2: minQ(complex(t1, t2), 721), lo, hi, xtol=1e-7)
        print('tau1=%.2f tau2*=%.5f |  inside(-0.004): %.3e  outside(+0.004): %.3e' % (
            t1, r, finite_min(complex(t1, r - 0.004), 64, 128, 0.05), finite_min(complex(t1, r + 0.004), 64, 128, 0.05)), flush=True)
    print('--- ratio Sigma/|q|^4 : Q(mirror) vs angle (constancy => proportional), selected tau ---')
    for tau in [-0.45 + 0.9j, -0.4 + 1.0j, 1j]:
        r = [S4(tau, th, 0.02) / Q(complex(-tau.real, tau.imag), th) for th in np.linspace(0.1, np.pi - 0.1, 12)]
        print(tau, 'ratio min %.4f max %.4f' % (min(r), max(r)))
    print('--- single strip s=0 : min_xi Sigma on 511 pts ---')
    xi = np.linspace(0, 1, 513)[1:-1]
    for tau in [hexa, -0.45 + 0.9j, -0.4 + 0.95j, 0.5 + 1.04j, 0.5 + 2.0j, 1j]:
        print(np.round(tau, 4), '%.3e' % sigma(tau, 0.0, xi).min())
    print('done %.0fs' % (time.time() - t0))
