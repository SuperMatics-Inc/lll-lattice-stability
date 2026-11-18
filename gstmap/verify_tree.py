"""Independent floating-point S sampling helpers; not a rigorous second cover."""
import json
from fractions import Fraction as F


def load(path):
    d = {}
    for ln in open(path):
        ln = ln.strip()
        if ln:
            e = json.loads(ln)
            d[(e['c1'], e['c2'], e['w'])] = e
    return list(d.values())


def job_s_replay(e):
    import indep
    c1, c2, w = float(F(e['c1'])), float(F(e['c2'])), float(F(e['w']))
    vals = indep.box_replay(c1, c2, w, 'S', n_rand=3, nrad=40, nth=48, ngrid=81)
    return ('Srep', e['c1'], e['c2'], e['w'], min(vals), all(v > 0 for v in vals))
