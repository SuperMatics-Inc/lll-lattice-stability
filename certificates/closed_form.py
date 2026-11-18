"""Closed Gaussian lattice-sum forms of the LLL linearised symbol for general tau.
NUMERICAL cross-check (float64) of the closed forms against pocket_scan.kern/sigma.

 A(s,xi) = 2*Khat1 - lam,   Khat1 = c * sum_{P,d in Z} exp(-pi t2[(P+s)^2+d^2]) exp(-2 pi i d (xi + t1 (P+s)))   [(P,d)-form]
 B(s,xi) = c * e^{4 pi i xi s} * sum_{x in s+Z, y in -s+Z} exp(-pi t2 (x^2+y^2) + 2 pi i t1 x y + 2 pi i xi (y-x))  [(x,y)-form]
 c = sqrt(t2)/pi,   Sigma = A^2 - |B|^2.
Run:  uv run --no-project --with numpy --with scipy python closed_form.py
"""
import numpy as np
from pocket_scan import kern, ft, sigma

def A_cf(tau, s, xi, M=14):
    t1, t2 = tau.real, tau.imag; c = np.sqrt(t2) / np.pi
    P, d = np.meshgrid(np.arange(-M, M + 1), np.arange(-M, M + 1), indexing='ij')
    x = P + s
    def K1hat(s_, xi_):
        x_ = P + s_
        return (c * np.exp(-np.pi * t2 * (x_ ** 2 + d ** 2) - 2j * np.pi * d * (xi_ + t1 * x_))).sum()
    return 2 * K1hat(s, xi).real - K1hat(0.0, 0.0).real, K1hat(s, xi).imag

def B_cf(tau, s, xi, M=14, sgn=+1):
    t1, t2 = tau.real, tau.imag; c = np.sqrt(t2) / np.pi
    m, n = np.meshgrid(np.arange(-M, M + 1), np.arange(-M, M + 1), indexing='ij')
    x = m + s; y = n - s
    return c * np.exp(4j * np.pi * xi * s) * np.exp(-np.pi * t2 * (x ** 2 + y ** 2) + sgn * 2j * np.pi * t1 * x * y + 2j * np.pi * xi * (y - x)).sum()

if __name__ == "__main__":
    rng = np.random.default_rng(1)
    for tau in [np.exp(2j * np.pi / 3) + 1.0 * 0, -0.4 + 0.95j, 0.4 + 1.0j, 1j]:
        if abs(tau - np.exp(2j * np.pi / 3)) < 1e-12: tau = -0.5 + np.sqrt(3) / 2 * 1j
        worst = {}
        for trial in range(12):
            s = rng.uniform(-.5, .5); xi = rng.uniform(-.5, .5)
            idx, K1, K2 = kern(tau, s)
            lam = ft(idx, kern(tau, 0.0)[1], 0.0)[0].real
            A_ref = (2 * ft(idx, K1, xi) - lam)[0]
            B_ref = ft(idx, K2, -np.atleast_1d(xi))[0]
            a, aim = A_cf(tau, s, xi)
            for sgn in (+1, -1):
                for sx in (+1, -1):
                    b = abs(B_cf(tau, s, sx * xi, sgn=sgn))
                    key = ('B', sgn, sx)
                    worst[key] = max(worst.get(key, 0), abs(b - abs(B_ref)))
            worst['A'] = max(worst.get('A', 0), abs(a - A_ref.real))
            worst['Aim'] = max(worst.get('Aim', 0), abs(aim))
        print(f"tau={tau}: ", {k: f'{v:.2e}' for k, v in worst.items()})
