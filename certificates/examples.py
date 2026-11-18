"""Archived table: full-torus (64x128, |q|>=0.05) min Sigma, small-k quartic ratio min Sigma/|q|^4 at |q|=0.03, and min_theta Q
for selected tau. NUMERICAL (float64)."""
import numpy as np
from pocket_scan import finite_min, polar_min, minQ

hexa = complex(-0.5, np.sqrt(3) / 2)
rows = [hexa, -0.45 + 0.9j, -0.4 + 0.95j, 0.5 + 1.04j, 0.4 + 1.0j,
        0.5 + 1.052j, 0.5 + 1.06j, 0.4 + 1.035j, 0.3 + 1.0j, 1j, 0.5 + 2j]
print("%-16s %14s %16s %12s" % ("tau", "min Sigma(64x128,|q|>=.05)", "min Sigma/|q|^4@.03", "min_th Q"))
for t in rows:
    fm = finite_min(t, 64, 128, 0.05)
    pm = polar_min(t, 0.03, 91)
    q = minQ(t, 361)
    print("%-16s %14.3e %16.4f %12.4f" % (np.round(t, 4), fm, pm, q), flush=True)
