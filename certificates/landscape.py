"""Sample Sigma(s,xi)=A^2-|B|^2 on the torus in double precision using the closed Gaussian sums
(A: 2*sum G(m+s,n) e^{2 pi i xi n} - lam ; B: sum G(m+s,n-s) e^{2 pi i xi (n-m)}, G(x,y)=c exp(-pi t2 (x^2+y^2)+2 pi i t1 x y)).
Used to design the Arb cell certificate (cert_sigma.py): reports min Sigma, min Sigma/rho^4 (rho=max(|s|,|xi|)) on rings.
Run: uv run --no-project --with numpy python landscape.py
"""
import numpy as np, sys

def sigma_grid(tau, S, X, M=7):
    t1, t2 = tau.real, tau.imag; c = np.sqrt(t2) / np.pi
    m = np.arange(-M, M + 1)[:, None, None, None]
    n = np.arange(-M, M + 1)[None, :, None, None]
    s = S[None, None, :, :]; xi = X[None, None, :, :]
    GA = np.exp(-np.pi * t2 * ((m + s) ** 2 + n ** 2) + 2j * np.pi * t1 * (m + s) * n + 2j * np.pi * xi * n)
    GB = np.exp(-np.pi * t2 * ((m + s) ** 2 + (n - s) ** 2) + 2j * np.pi * t1 * (m + s) * (n - s) + 2j * np.pi * xi * (n - m))
    lam = (c * np.exp(-np.pi * t2 * (m[:, :, 0, 0] ** 2 + n[:, :, 0, 0] ** 2) + 2j * np.pi * t1 * m[:, :, 0, 0] * n[:, :, 0, 0])).sum().real
    A = 2 * (c * GA).sum((0, 1)).real - lam
    B = (c * GB).sum((0, 1))
    return lam, A, B, A ** 2 - np.abs(B) ** 2

if __name__ == "__main__":
    for tau in [-0.4 + 0.95j, 0.4 + 1.0j, np.exp(2j * np.pi / 3), 0.5 + 1.04j]:
        tau = complex(tau)
        N = 400
        g = (np.arange(N) + 0.5) / N - 0.5
        S, X = np.meshgrid(g, g, indexing='ij')
        lam, A, B, Sg = sigma_grid(tau, S, X)
        rho = np.maximum(np.abs(S), np.abs(X))
        print(f"\ntau={tau} lam={lam:.10f} A range [{A.min():.4f},{A.max():.4f}] |B| range [{abs(B).min():.4f},{abs(B).max():.4f}]")
        for lo, hi in [(0.0, 0.02), (0.02, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.5)]:
            mk = (rho >= lo) & (rho < hi)
            r = Sg[mk] / rho[mk] ** 4
            i = np.argmin(Sg[mk] / rho[mk] ** 4)
            print(f"  rho in [{lo},{hi}): min Sigma={Sg[mk].min():.4e}  min Sigma/rho^4={r.min():.4f} max={r.max():.3f}  neg count={(Sg[mk]<0).sum()}")
        sys.stdout.flush()
