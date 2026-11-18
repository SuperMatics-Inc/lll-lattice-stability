"""Check identities for the LLL Bloch symbol numerically, using pocket_scan in double and multiprecision arithmetic.
 1. A(k)=A(-k) real, |B(k)|=|B(-k)|, B(k)conj(B(-k)) = |B|^2 real  (pairing (s,xi)<->(-s,-xi))
 2. Poisson/Gaussian-lattice-sum form of K1hat^s(xi)
 3. SL2(Z) covariance of Sigma
 4. Rectangular lattices tau=i*b: min_theta Q<0 and strip Sigma<0 (GST thmunstable consistency)
 5. Hex strip symbol sqrt(Sigma(s=0, xi=1/2)) vs certified mu(1/2) = 0.0957890177120495 (cert_mu.log)
"""
import numpy as np
from pocket_scan import kern, ft, minQ, sigma

rng = np.random.default_rng(3)
print("== 1. symmetries of the symbol")
worst = [0, 0, 0]
for tau in [-0.5 + np.sqrt(3) / 2 * 1j, 0.35 + 0.99j, -0.2 + 1.05j, 0.1 + 1.3j]:
    for _ in range(3):
        s = rng.uniform(-0.45, 0.45); xi = rng.uniform(-0.45, 0.45)
        idx, K1, K2 = kern(tau, s); idm, K1m, K2m = kern(tau, -s)
        lam = ft(idx, kern(tau, 0.0)[1], 0.0)[0].real
        A = 2 * ft(idx, K1, xi)[0] - lam
        Am = 2 * ft(idm, K1m, -xi)[0] - lam
        B = ft(idx, K2, -xi)[0]; Bm = ft(idm, K2m, xi)[0]
        prod = B * np.conj(Bm)
        worst[0] = max(worst[0], abs(A - Am), abs(A.imag))
        worst[1] = max(worst[1], abs(abs(B) - abs(Bm)))
        worst[2] = max(worst[2], abs(prod.imag), abs(prod.real - abs(B) ** 2))
print("max |A(k)-A(-k)|,|Im A| = %.2e ; max ||B(k)|-|B(-k)|| = %.2e ; max |Im prod|,|Re prod-|B|^2| = %.2e" % tuple(worst))
assert max(worst) < 1e-14

print("== 2. Poisson form of K1hat")
def gl(tau, s, xi, N=12):
    m, d = np.meshgrid(np.arange(-N, N + 1), np.arange(-N, N + 1)); m = m.ravel(); d = d.ravel()
    w = np.exp(-np.pi * np.abs(m - d * tau) ** 2 / tau.imag)
    return (w * np.exp(-2j * np.pi * (d * xi + m * s))).sum() / np.pi
mx = 0
for tau in [complex(-0.5, np.sqrt(3) / 2), 0.3 + 1.1j, -0.4 + 0.95j]:
    for s, xi in [(0.0, 0.13), (0.0, 0.37), (0.21, 0.13), (0.4, 0.3)]:
        idx, K1, K2 = kern(tau, s)
        mx = max(mx, abs(ft(idx, K1, xi)[0] - gl(tau, s, xi)))
print("max |K1hat - lattice sum| over 12 points = %.2e" % mx)
assert mx < 1e-12

print("== 3. SL2 covariance of Sigma")
def sigma_k(tau, s, xi):
    return float(np.real(sigma(tau, s, np.array([xi]))[0]))
mx1 = mx2 = 0
for tau in [0.35 + 0.99j, 0.1 + 1.3j, -0.2 + 1.05j]:
    for _ in range(3):
        s = rng.uniform(-0.45, 0.45); xi = rng.uniform(-0.45, 0.45)
        S0 = sigma_k(tau, s, xi)
        mx1 = max(mx1, abs(sigma_k(tau + 1, s, xi - s) - S0))
        mx2 = max(mx2, abs(sigma_k(-1 / tau, xi, -s) - S0))
print("max |Sigma_{tau+1}(s,xi-s)-Sigma_tau| = %.2e ; max |Sigma_{-1/tau}(xi,-s)-Sigma_tau| = %.2e" % (mx1, mx2))
assert max(mx1, mx2) < 1e-12

print("== 4. rectangular lattices tau=i b")
bs = sorted(set(list(np.linspace(0.5, 3, 26)) + [np.sqrt(3) / 2, np.sqrt(3)]))
mq = max(minQ(complex(0, b), 361) for b in bs)
print("max over b in [0.5,3] (28 values) of min_theta Q = %.4f (<0 => unstable)" % mq)
assert mq < 0
xi = np.linspace(0, 1, 513)[1:-1]
for b in (0.6, 0.8, 1.0, 1.3, 1.7, 2.0):
    m = sigma(complex(0, b), 0.0, xi).min()
    print("strip s=0 min Sigma tau=%gi: %.2e" % (b, m))
    assert m < 0
print("== 5. hex strip symbol vs certified mu(1/2)")
hexa = complex(-0.5, np.sqrt(3) / 2)
mu_half = np.sqrt(sigma(hexa, 0.0, np.array([0.5]))[0])
ref = 0.0957890177120495
print("sqrt(Sigma_hex(0,1/2)) = %.16f ; certified mu(1/2) = %.16f ; diff = %.2e" % (mu_half, ref, abs(mu_half - ref)))
assert abs(mu_half - ref) < 1e-13
print("ALL OK")
