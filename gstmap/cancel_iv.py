"""Independent mpmath.iv scalar-q Taylor certificate, with exact low-order identities.
No Arb, gsteng, boxcert, unst or gstsym imports. See PROOFS.md section 9.
Usage: python cancel_iv.py tree.jsonl [workers=4] [limit=0]
Verifies supplied rational q0 over every unstable leaf by interval series.
"""
import json, sys, time, math
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor
from mpmath import iv
iv.dps = 35


def I(x):
    if hasattr(x, '_mpi_'):
        return x
    x = F(x)
    return iv.mpf(x.numerator) / x.denominator


def hull(a, b):
    return iv.mpf([I(a).a, I(b).b])


def rad(e):
    return iv.mpf([-e.b, e.b])


def mag(z):
    return abs(z).b


def series(L, Q, N):
    out = [iv.mpc(1), L]
    for j in range(2, N + 1):
        out.append((L * out[-1] + 2 * Q * out[-2]) / j)
    return out


def rem(l, q, N):
    # Coefficientwise majorant exp(l*t + q*t^2), |t|<=1.
    # exp(x)-sum_{0..k} x^j/j! <= exp(x)*x^(k+1)/(k+1)! for x>=0.
    kmax = N // 2
    out = iv.mpf(0)
    for k in range(kmax + 1):
        j = N - 2 * k + 1
        out += q**k / math.factorial(k) * l**j / math.factorial(j)
    return iv.exp(l) * out + iv.exp(l + q) * q**(kmax + 1) / math.factorial(kmax + 1)


def tail(x, y, s, xi, M, half=False):
    # Simultaneously bounds omitted K/B values on complex |t|<=1 and hence,
    # by Cauchy, each omitted Taylor coefficient. Norms replaced by l1 bounds.
    a = iv.pi * y.a / (2 if half else 1)
    beta = (2 * iv.pi * (y.b * abs(s) + abs(x * s + xi))).b
    q = (2 * iv.pi * (y.b + abs(x).b) * s*s).b
    r0 = iv.exp(-3*a + beta)
    rM = iv.exp(-(2*M + 3)*a + beta)
    if not (r0.b < 1 and rM.b < 1):
        raise ValueError('tail geometric ratios do not contract')
    full = 1 + 2*iv.exp(-a + beta)/(1-r0)
    omitted = 2*iv.exp(-a*(M+1)**2 + beta*(M+1))/(1-rM)
    return (2 * iv.sqrt(y).b / iv.pi * iv.exp(q) * omitted * full).b


def coefficients(c1, c2, w, s, xi, N=10, M=4, deriv=False, remainder=True, wy=None):
    x = hull(c1-w, c1+w)
    wy = w if wy is None else wy
    y = hull(c2-wy, c2+wy)
    s, xi = I(s), I(xi)
    c = iv.sqrt(y) / iv.pi
    K = [iv.mpc(0) for _ in range(N+1)]
    B = [iv.mpc(0) for _ in range(N+1)]
    epsK = iv.mpf(0); epsB = iv.mpf(0)
    dK = [[iv.mpc(0) for _ in range(5)] for _ in range(2)]
    dB = [[iv.mpc(0) for _ in range(5)] for _ in range(2)]
    for m in range(-M, M+1):
        for n in range(-M, M+1):
            weight = c * iv.exp(iv.mpc(-iv.pi*y*(m*m+n*n), 2*iv.pi*x*m*n))
            LK = iv.mpc(-2*iv.pi*y*m*s, 2*iv.pi*n*(x*s+xi))
            QK = -iv.pi*y*s*s
            LB = iv.mpc(2*iv.pi*y*(n-m)*s, 2*iv.pi*(n-m)*(x*s+xi))
            QB = iv.mpc(-2*iv.pi*y*s*s, -2*iv.pi*x*s*s)
            kser = series(LK, QK, N); bser = series(LB, QB, N)
            for j in range(0, N+1, 2):
                K[j] += weight*kser[j]
                B[j] += weight*bser[j]
            if deriv:
                dw = (iv.mpc(0,2*iv.pi*m*n), 1/(2*y)-iv.pi*(m*m+n*n))
                dlK = (iv.mpc(0,2*iv.pi*n*s), -2*iv.pi*m*s)
                dqK = (iv.mpc(0), -iv.pi*s*s)
                dlB = (iv.mpc(0,2*iv.pi*(n-m)*s), 2*iv.pi*(n-m)*s)
                dqB = (iv.mpc(0,-2*iv.pi*s*s), -2*iv.pi*s*s)
                for axis in range(2):
                    for j in (0,2,4):
                        dk = dw[axis]*kser[j]
                        db = dw[axis]*bser[j]
                        if j:
                            dk += dlK[axis]*kser[j-1]+dqK[axis]*kser[j-2]
                            db += dlB[axis]*bser[j-1]+dqB[axis]*bser[j-2]
                        dK[axis][j] += weight*dk
                        dB[axis][j] += weight*db
            # Use exact-parameter upper magnitude bounds, not the interval weight
            # magnitude (which would include arbitrary complex phase wrapping).
            if remainder:
                wm = (c.b * iv.exp(-iv.pi*y.a*(m*m+n*n))).b
                epsK += wm * rem(mag(LK), mag(QK), N)
                epsB += wm * rem(mag(LB), mag(QB), N)
    T = tail(x, y, s, xi, M)
    # Infinite tail: value plus the first N+1 Taylor coefficients removed.
    epsK = (epsK + (N+2)*T).b
    epsB = (epsB + (N+2)*T).b
    for j in range(0, N+1, 2):
        K[j] += iv.mpc(rad(T), rad(T))
        B[j] += iv.mpc(rad(T), rad(T))
    lam = K[0].real
    A = [2*k.real for k in K]
    A[0] = lam; B[0] = iv.mpc(lam, 0)
    grad = None
    if deriv:
        Td = (tail(x,y,s,xi,M,half=True) *
              (4/(y.a*iv.e)+4*iv.pi*s*s+1/(2*y.a))).b
        for axis in range(2):
            for j in (0,2,4):
                dK[axis][j] += iv.mpc(rad(Td),rad(Td))
                dB[axis][j] += iv.mpc(rad(Td),rad(Td))
        grad = [2*dK[a][0].real*(A[4]-B[4].real)
                +2*lam*(2*dK[a][4].real-dB[a][4].real)
                -2*B[2].imag*dB[a][2].imag for a in range(2)]
    return lam, A, B, 2*epsK, epsB, grad


def upper(c1, c2, w, q0, N=10, M=4):
    s,xi=F(q0[0]),F(q0[1])
    lam, A, B, epsA, epsB, grad = coefficients(c1,c2,w,s,xi,N,M,deriv=True)
    # Centre only the quartic; all remaining coefficients retain rho^d scaling.
    lc,Ac,Bc,*_=coefficients(c1,c2,F(0),s,xi,4,M,remainder=False)
    ans = 2*lc*(Ac[4]-Bc[4].real)-Bc[2].imag**2
    ans += I(w)*(rad(abs(grad[0]))+rad(abs(grad[1])))
    for d in range(6, 2*N+1, 2):
        for j in range(max(0,d-N), min(N,d)+1, 2):
            k = d-j
            ans += A[j]*A[k] - B[j].real*B[k].real - B[j].imag*B[k].imag
    absA = sum(abs(a).b for a in A)
    absB = sum(mag(b) for b in B)
    err = 2*absA*epsA + epsA*epsA + 2*absB*epsB + epsB*epsB
    return (ans + err).b


def verify(e, maxdepth=7):
    t0=time.time(); n=0; worst=-float('inf')
    c1,c2,w = (F(e[k]) for k in ('c1','c2','w'))
    stack=[(c1,c2,w,0)]
    try:
        while stack:
            a,b,h,d=stack.pop()
            up=upper(a,b,h,e['q0']); n+=1
            if up < 0:
                worst=max(worst,float(up)); continue
            if d>=maxdepth:
                return dict(ok=False,c1=e['c1'],c2=e['c2'],w=e['w'],upper=float(up),n=n,seconds=time.time()-t0)
            hh=h/2
            stack.extend((a+i*hh,b+j*hh,hh,d+1) for i in (-1,1) for j in (-1,1))
        return dict(ok=True,c1=e['c1'],c2=e['c2'],w=e['w'],upper=worst,n=n,seconds=time.time()-t0)
    except Exception as ex:
        return dict(ok=False,c1=e['c1'],c2=e['c2'],w=e['w'],error=repr(ex),n=n,seconds=time.time()-t0)


if __name__=='__main__':
    ev={}
    for ln in open(sys.argv[1]):
        e=json.loads(ln); ev[(e['c1'],e['c2'],e['w'])]=e
    leaves=sorted([e for e in ev.values() if e['ev']=='U'],key=lambda e:-e['m4c'])
    workers=int(sys.argv[2]) if len(sys.argv)>2 else 4
    limit=int(sys.argv[3]) if len(sys.argv)>3 else 0
    if limit: leaves=leaves[:limit]
    t0=time.time(); fails=nodes=0; worst=-float('inf')
    with ProcessPoolExecutor(workers) as ex:
        for r in ex.map(verify,leaves):
            print(json.dumps(r),flush=True)
            fails+=not r['ok']; nodes+=r['n']
            if r['ok']: worst=max(worst,r['upper'])
    print(json.dumps(dict(summary=True,n=len(leaves),fails=fails,nodes=nodes,worst_upper=worst,seconds=time.time()-t0)),flush=True)
    sys.exit(bool(fails))
