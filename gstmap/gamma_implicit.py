"""Enclose the quartic transition curve Gamma4 with a parametric Krawczyk argument.
The unknowns are (y,r), f(x,y,r)=P4_(x+iy)(1,r), F=(f,f_r).
We interpolate rational matrices and use mpmath.iv for the bounds. The starting guesses
are floating-point values. The calculation locates a long-wave transition.
Usage: python gamma_implicit.py [workers=4] [first_cells=0]
"""
import json, sys, time, os
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import root
import gstsym
import cancel_iv as C
from mpmath import iv
S=F(1,64)
NODES=[F(-1),F(-1,2),F(0),F(1,2),F(1)]


def inverse(A):
    n=len(A); A=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for j in range(n):
        k=next(k for k in range(j,n) if A[k][j])
        A[j],A[k]=A[k],A[j]; a=A[j][j]; A[j]=[x/a for x in A[j]]
        for k in range(n):
            if k!=j:
                a=A[k][j]; A[k]=[x-a*y for x,y in zip(A[k],A[j])]
    return [row[n:] for row in A]


V= inverse([[r**k for k in range(5)] for r in NODES])


def interpolate(vals):
    return [sum(C.I(v)*a for v,a in zip(row,vals)) for row in V]


def poly(c,r):
    out=iv.mpf(0)
    for a in reversed(c): out=out*r+a
    return out


def diff(c): return [i*c[i] for i in range(1,len(c))]


def quartic(x,y,hx,hy,deriv):
    vals=[]; gx=[]; gy=[]
    for r in NODES:
        l,A,B,_,_,g=C.coefficients(x,y,hx,S,S*r,4,4,deriv=deriv,remainder=False,wy=hy)
        vals.append((2*l*(A[4]-B[4].real)-B[2].imag**2)/C.I(S)**4)
        if deriv:
            gx.append(g[0]/C.I(S)**4); gy.append(g[1]/C.I(S)**4)
    return interpolate(vals),interpolate(gx) if deriv else None,interpolate(gy) if deriv else None


def seed(x):
    def fun(z):
        c=gstsym.homog(gstsym.sigma_taylor(complex(float(x),z[0]),4)[3],4).real
        return [np.polynomial.polynomial.polyval(z[1],c), np.polynomial.polynomial.polyval(z[1],np.polynomial.polynomial.polyder(c))]
    ans=root(fun,[1.0491262-2.3*(float(x)-.5)**2,-.5+1.6*(.5-float(x))])
    if not ans.success or max(abs(np.array(fun(ans.x))))>1e-7: raise ValueError('seed failed')
    y,r=map(lambda z:F(round(float(z)*2**40),2**40),ans.x)
    return y,r


def cell(x,hx):
    t0=time.time(); y,r=seed(x)
    hy=4*hx; hr=8*hx
    c0,_,_=quartic(x,y,F(0),F(0),False)
    cb,cx,cy=quartic(x,y,hx,hy,True)
    cb=[c0[k]+C.I(hx)*C.rad(abs(cx[k]))+C.I(hy)*C.rad(abs(cy[k])) for k in range(5)]
    rr=C.hull(r-hr,r+hr); r0=C.I(r)
    J=[[poly(cy,rr),poly(diff(cb),rr)], [poly(diff(cy),rr),poly(diff(diff(cb)),rr)]]
    # Any exact rational approximate inverse is admissible.
    JJ=np.array([[float(poly(cy,r0).mid),float(poly(diff(c0),r0).mid)],
                 [float(poly(diff(cy),r0).mid),float(poly(diff(diff(c0)),r0).mid)]])
    R=[[F(float(v)).limit_denominator(2**30) for v in row] for row in np.linalg.inv(JJ)]
    F0=[poly(c0,r0)+C.I(hx)*C.rad(abs(poly(cx,r0))),
        poly(diff(c0),r0)+C.I(hx)*C.rad(abs(poly(diff(cx),r0)))]
    E=[[iv.mpf(i==j)-sum(C.I(R[i][k])*J[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    radii=[C.I(hy),C.I(hr)]
    K=[-sum(C.I(R[i][k])*F0[k] for k in range(2))+sum(E[i][j]*C.rad(radii[j]) for j in range(2)) for i in range(2)]
    cn=[sum(abs(E[i][j]).b*radii[j]/radii[i] for j in range(2)) for i in range(2)]
    contraction=max(float(v.b) for v in cn)
    inclusion=all(abs(K[i]).b<radii[i] for i in range(2))
    # A quartic with a double real root is (r-r*)^2 times this quadratic.
    qa=cb[4]; qb=cb[3]+2*rr*cb[4]; qc=cb[2]+2*rr*cb[3]+3*rr**2*cb[4]
    discr=4*qa*qc-qb**2
    l,A,B,*_=C.coefficients(x,y,hx,S,C.I(S)*rr,6,4,remainder=False,wy=hy)
    p6=(2*l*(A[6]-B[6].real)+2*A[2]*(A[4]-B[4].real)-2*B[2].imag*B[4].imag)/C.I(S)**6
    yp=-poly(cx,rr)/J[0][0]
    rp=-(poly(diff(cx),rr)+J[1][0]*yp)/J[1][1]
    dc=[cx[k]+cy[k]*yp for k in range(5)]
    db=dc[3]+2*rp*cb[4]+2*rr*dc[4]
    dq=dc[2]+2*rp*cb[3]+2*rr*dc[3]+6*rr*rp*cb[4]+3*rr**2*dc[4]
    discr_deriv=4*dc[4]*qc+4*qa*dq-2*qb*db
    norm_deriv=2*C.hull(x-hx,x+hx)+2*C.hull(y-hy,y+hy)*yp
    norm_at_x=C.I(x)**2+C.hull(y-hy,y+hy)**2-1
    # Below the unit arc this analytic branch is not the global minimum branch.
    # Near the arc, D' > 0 and (|tau|^2)' > 0 transfer the exact D=0 identity
    # at the unique unit-arc crossing to D>0 on the F+ side (PROOFS section 10).
    near_arc=x-hx<F(3,10)
    quotient_ok=(discr_deriv.a>0 and norm_deriv.a>0) if near_arc else discr.a>0
    ok=inclusion and all(v.b<1 for v in cn) and J[0][0].b<0 and J[1][1].a>0 and qa.a>0 and quotient_ok and p6.a>0
    return dict(ok=bool(ok),x=str(x),hx=str(hx),y=str(y),hy=str(hy),r=str(r),hr=str(hr),
                contraction=contraction,inclusion=bool(inclusion),fy_upper=float(J[0][0].b),frr_lower=float(J[1][1].a),
                leading_lower=float(qa.a),quotient_discr_lower=float(discr.a),p6_lower=float(p6.a),
                discr_deriv_lower=float(discr_deriv.a),norm_deriv_lower=float(norm_deriv.a),
                norm_at_x_lo=float(norm_at_x.a),norm_at_x_hi=float(norm_at_x.b),seconds=time.time()-t0)


def run_cell(args):
    x,hx,depth=args
    try: out=cell(x,hx)
    except Exception as ex: out=dict(ok=False,x=str(x),hx=str(hx),error=repr(ex))
    out['depth']=depth
    return out


if __name__=='__main__':
    if os.environ.get('SM_OUT'):
        sys.stdout=open(os.path.join(os.environ['SM_OUT'],'gamma_implicit.jsonl'),'w',buffering=1)
    workers=int(sys.argv[1]) if len(sys.argv)>1 else 4
    limit=int(sys.argv[2]) if len(sys.argv)>2 else 0
    # Exact closed parameter interval [37/128,1/2], disjoint interiors.
    todo=[(F(37,128)+F(2*i+1,2048),F(1,2048),0) for i in range(216)]
    if limit: todo=todo[:limit]
    t0=time.time(); leaves=[]; failed=[]
    with ProcessPoolExecutor(workers) as ex:
        while todo:
            out=list(ex.map(run_cell,todo)); todo=[]
            for d in out:
                if d['ok']: leaves.append(d); print(json.dumps(d),flush=True)
                elif d['depth']<7:
                    x,h=F(d['x']),F(d['hx'])/2
                    todo.extend((x+i*h,h,d['depth']+1) for i in (-1,1))
                else: failed.append(d); print(json.dumps(d),flush=True)
    summary=dict(summary=True,n=len(leaves),fails=len(failed),seconds=time.time()-t0,
                 fy_upper=max([d['fy_upper'] for d in leaves],default=0),
                 frr_lower=min([d['frr_lower'] for d in leaves],default=0),
                 p6_lower=min([d['p6_lower'] for d in leaves],default=0),
                 max_hy=max([float(F(d['hy'])) for d in leaves],default=0))
    print(json.dumps(summary),flush=True)
    sys.exit(bool(failed))
