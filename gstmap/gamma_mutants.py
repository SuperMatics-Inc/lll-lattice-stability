"""Reject false root/patch-join claims; controls use exact rational perturbations."""
import json,sys,os
from fractions import Fraction as F
import gamma_implicit as G
import verify_gamma_chain as V

if __name__=='__main__':
    if os.environ.get('SM_OUT'):
        sys.stdout=open(os.path.join(os.environ['SM_OUT'],'gamma_mutants.log'),'w',buffering=1)
    x=F(2,5); y,r=G.seed(x); hy=F(1,2**22); hr=2*hy
    cases=[]
    for name,yy,rr,expected in [
        ('valid double root',y,r,True),
        ('wrong y centre',y+F(1,1000),r,False),
        ('wrong null direction',y,r+F(1,1000),False)]:
        ans=V.root_in_box(x,yy,hy,rr,hr)
        cases.append((name,ans['ok'],expected))
    source=sys.argv[1] if len(sys.argv)>1 else 'gamma_implicit.jsonl'
    cells=[e for ln in open(source) if (e:=json.loads(ln)).get('ok') is True]
    cells.sort(key=lambda e:F(e['x']))
    cases.append(('valid adjacent join',V.glue(cells[:2])['ok'],True))
    cases.append(('missing intermediate patch',V.glue([cells[0],cells[2]])['ok'],False))
    shifted=dict(cells[1]); shifted['x']=str(F(shifted['x'])+F(1,2**20))
    cases.append(('shifted common endpoint',V.glue([cells[0],shifted])['ok'],False))
    for name,actual,expected in cases:
        print(json.dumps(dict(case=name,accepted=actual,expected=expected)))
    fails=sum(a!=e for _,a,e in cases)
    print(json.dumps(dict(summary=True,n=len(cases),fails=fails)))
    sys.exit(bool(fails))
