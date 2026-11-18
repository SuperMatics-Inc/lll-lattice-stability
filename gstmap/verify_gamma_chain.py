"""Glue saved parametric Krawczyk patches by an existence proof in each overlap.
Also enclose the unique unit-arc crossing by exact rational x endpoints.
Usage: python verify_gamma_chain.py gamma_implicit.jsonl [workers=4].
"""
import json,sys,os,time,math
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from mpmath import iv
import gamma_implicit as G
import cancel_iv as C


def box(e):
    return [F(e[k]) for k in ('x','hx','y','hy','r','hr')]


def root_in_box(x,y,hy,r,hr):
    c0,_,_=G.quartic(x,y,F(0),F(0),False)
    _,cx,cy=G.quartic(x,y,F(0),hy,True)
    cb=[a+C.I(hy)*C.rad(abs(d)) for a,d in zip(c0,cy)]
    rr=C.hull(r-hr,r+hr); r0=C.I(r)
    J=[[G.poly(cy,rr),G.poly(G.diff(cb),rr)],
       [G.poly(G.diff(cy),rr),G.poly(G.diff(G.diff(cb)),rr)]]
    JJ=np.array([[float(v.mid) for v in row] for row in J])
    R=[[F(float(v)).limit_denominator(2**30) for v in row] for row in np.linalg.inv(JJ)]
    f=[G.poly(c0,r0),G.poly(G.diff(c0),r0)]
    E=[[iv.mpf(i==j)-sum(C.I(R[i][k])*J[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    radii=[C.I(hy),C.I(hr)]
    K=[-sum(C.I(R[i][k])*f[k] for k in range(2))+sum(E[i][j]*C.rad(radii[j]) for j in range(2)) for i in range(2)]
    cn=[sum(abs(E[i][j]).b*radii[j]/radii[i] for j in range(2)) for i in range(2)]
    ok=all(abs(K[i]).b<radii[i] and cn[i].b<1 for i in range(2))
    yy=C.I(y)+K[0]
    rr=C.I(r)+K[1]
    norm=C.I(x)**2+yy**2-1
    lower=lambda v:math.nextafter(float(v.a),-math.inf)
    upper=lambda v:math.nextafter(float(v.b),math.inf)
    return dict(ok=bool(ok),x=str(x),hy=str(hy),hr=str(hr),
                contraction=max(float(v.b) for v in cn),norm_lo=lower(norm),norm_hi=upper(norm),
                y_lo=lower(yy),y_hi=upper(yy),r_lo=lower(rr),r_hi=upper(rr))


def glue(pair):
    a,b=map(box,pair); x=a[0]+a[1]
    if x!=b[0]-b[1]: return dict(ok=False,error='x gap/overlap')
    yl=max(a[2]-a[3],b[2]-b[3]); yh=min(a[2]+a[3],b[2]+b[3])
    rl=max(a[4]-a[5],b[4]-b[5]); rh=min(a[4]+a[5],b[4]+b[5])
    if not (yl<yh and rl<rh): return dict(ok=False,error='empty intersection')
    return root_in_box(x,(yl+yh)/2,(yh-yl)/2,(rl+rh)/2,(rh-rl)/2)


def endpoint(x,cells,bits=22):
    y,r=G.seed(x); hy=F(1,2**bits); hr=2*hy
    cs=[box(e) for e in cells if F(e['x'])-F(e['hx'])<=x<=F(e['x'])+F(e['hx'])]
    contained=any(abs(y-a[2])+hy<a[3] and abs(r-a[4])+hr<a[5] for a in cs)
    ans=root_in_box(x,y,hy,r,hr); ans['contained']=contained
    ans['ok']=ans['ok'] and contained
    return ans


if __name__=='__main__':
    src=sys.argv[1]; nw=int(sys.argv[2]) if len(sys.argv)>2 else 4
    if os.environ.get('SM_OUT'):
        sys.stdout=open(os.path.join(os.environ['SM_OUT'],'gamma_chain.jsonl'),'w',buffering=1)
    cells=[e for ln in open(src) if (e:=json.loads(ln)).get('ok') is True]
    cells.sort(key=lambda e:F(e['x']))
    assert F(cells[0]['x'])-F(cells[0]['hx'])==F(37,128)
    assert F(cells[-1]['x'])+F(cells[-1]['hx'])==F(1,2)
    assert all(e['norm_deriv_lower']>0 for e in cells if F(e['x'])-F(e['hx'])<F(3,10))
    assert all(-1<F(e['r'])-F(e['hr'])<F(e['r'])+F(e['hr'])<0 for e in cells)
    t=time.time()
    with ProcessPoolExecutor(nw) as ex:
        results=list(ex.map(glue,zip(cells,cells[1:])))
    for r in results: print(json.dumps(r),flush=True)
    end=[endpoint(F(297,1000),cells),endpoint(F(2971,10000),cells)]
    for r in end: print(json.dumps(dict(endpoint=True,**r)),flush=True)
    point=[endpoint(x,cells,30) for x in [F(3,10),F(2,5),F(1,2)]]
    for r in point: print(json.dumps(dict(point=True,**r)),flush=True)
    fails=sum(not r['ok'] for r in results+end+point)
    fails+=not (end[0]['norm_hi']<0<end[1]['norm_lo'])
    print(json.dumps(dict(summary=True,patches=len(cells),glues=len(results),fails=fails,
                         endpoint_bracket=['297/1000','2971/10000'],seconds=time.time()-t)),flush=True)
    sys.exit(bool(fails))
