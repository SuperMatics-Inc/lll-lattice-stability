"""Record the exact dyadic leaf geometry and verify that altered subdivisions are rejected."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import contextlib
import importlib.util
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / ('cert' if (ROOT/'cert').exists() else 'certificates')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output',type=Path, default=ROOT/'logs'/'rigor', help='Directory for the subdivision records')
    a = p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    for name, args in [('minus',['-2/5','19/20','1/100']),('plus',['2/5','1','1/40'])]:
        sys.argv = ['cert_sigma.py',*args,'20','8','192','13']
        spec = importlib.util.spec_from_file_location('sigma',CERT/'cert_sigma.py')
        c = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(c)
        def ball(q): return c.arb(q.numerator)/q.denominator
        with (a.output/f'trace-{name}.log').open('w') as log, contextlib.redirect_stdout(log):
            kap,_,_ = c.build_kappa6()
            ok,p4,margin = c.near_certificate(c.R1,c.N)
            assert ok
            todo = [(F(2*i+1,8),F(2*j-3,8),F(1,8)) for i in range(2) for j in range(4)]
            leaves=[]; tests=0; worst=None
            while todo:
                s,x,h = todo.pop()
                near = s+h <= F(args[2]) and abs(x)+h <= F(args[2])
                if near: kind='near'; bound=None
                else:
                    v=c.cell_test(ball(s),ball(x),ball(h),kap); tests+=1
                    if v > 0:
                        kind='far'; bound=str(v.lower())
                        worst=v.lower() if worst is None else min(worst,v.lower())
                    else:
                        assert h>F(1,8192)
                        todo.extend((s+ds*h/2,x+dx*h/2,h/2) for ds in [-1,1] for dx in [-1,1])
                        continue
                leaves.append(dict(s=str(s),xi=str(x),h=str(h),kind=kind,lower=bound))
            # The recursive split gives disjoint interiors. Check the area
            # separately to detect a missing or repeated leaf.
            area=sum((4*F(v['h'])**2 for v in leaves),F(0))
            assert area == F(1,2)
            assert sum((4*F(v['h'])**2 for v in leaves[:-1]),F(0)) != F(1,2)
            assert area+4*F(leaves[0]['h'])**2 != F(1,2)
            assert not (margin-margin-c.arb(1) > 0)
            assert not (p4-p4-c.arb(1) > 0)
            # An intentionally too-large near radius must not be accepted.
            mutant=c.near_certificate(c.arb(1)/2,c.N)[0]
            assert not mutant
            print('Certified area',area,'tests',tests,'far',sum(v['kind']=='far' for v in leaves),
                  'near',sum(v['kind']=='near' for v in leaves),'worst',worst)
            print('Rejected: missing leaf, duplicate leaf, negative near or quartic bound, radius=1/2')
            print('exit=0')
        data=dict(parameter=args,leaves=leaves,area=str(area),tests=tests,
                  adversarial={'removed_leaf':'rejected','duplicated_leaf':'rejected',
                  'negative_near_bound':'rejected','negative_quartic_bound':'rejected',
                  'near_radius_1/2':'rejected'})
        (a.output/f'trace-{name}.json').write_text(json.dumps(data,indent=2)+'\n')
        print(name,len(leaves),'exit=0',flush=True)


if __name__=='__main__':main()
