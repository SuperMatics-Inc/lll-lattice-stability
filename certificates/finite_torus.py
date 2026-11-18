"""Check the LLL spectrum numerically using a separate physical-space Galerkin calculation.
No Gaussian Bloch symbol or cert_sigma imports. Primitive magnetic translations
pair alpha,-alpha. Their pointwise phases cancel in the anomalous integral.
Memory is bounded by batching; M=256 means 65,536 flux quanta.
"""
import argparse
import json
import numpy as np


def wave(z,tau,K):
    a=np.sqrt(np.pi/tau.imag); h=np.pi/a
    x=z.real; y=z.imag
    centre=np.floor(-y/h).astype(int)
    out=np.zeros_like(z,dtype=complex)
    for k in range(-K,K+1):
        n=centre+k
        out+=np.exp(-(y+n*h)**2+1j*np.pi*tau.real*n*n+2j*np.pi*n*x/a)
    return (2/(np.pi*a*a))**0.25*np.exp(1j*x*y)*out


def check(tau,M,Q,K):
    a=np.sqrt(np.pi/tau.imag)
    v=(np.arange(Q)+0.5)/Q
    u,w=np.meshgrid(v,v,indexing='ij')
    z=(a*(u+tau*w)).reshape(1,-1)
    u0=wave(z,tau,K)
    weight=np.pi/(Q*Q)
    norm=float(weight*np.sum(abs(u0)**2))
    lam=float(weight*np.sum(abs(u0)**4)/norm)
    m=np.arange(-M//2,M//2)
    p,q=np.meshgrid(m,m,indexing='ij')
    alpha=(a*(p+tau*q)/M).ravel()
    max_growth=0.; min_frequency=float('inf'); max_norm_error=0.
    for start in range(0,len(alpha),128):
        av=alpha[start:start+128,None]
        up=wave(z+av,tau,K); um=wave(z-av,tau,K)
        dp=weight*np.sum(abs(up)**2*abs(u0)**2,axis=1)
        dm=weight*np.sum(abs(um)**2*abs(u0)**2,axis=1)
        e=weight*np.sum(np.conj(up)*np.conj(um)*u0*u0,axis=1)
        H=np.zeros((len(av),2,2),complex)
        H[:,0,0]=2*dp-lam; H[:,1,1]=-(2*dm-lam)
        H[:,0,1]=e; H[:,1,0]=-np.conj(e)
        eig=np.linalg.eigvals(H)
        nonzero=abs(av[:,0])>1e-14
        max_growth=max(max_growth,float(np.max(abs(eig[nonzero].imag))))
        min_frequency=min(min_frequency,float(np.min(abs(eig[nonzero]))))
        max_norm_error=max(max_norm_error,float(np.max(abs(weight*np.sum(abs(up)**2,axis=1)-1))))
    return dict(tau=[tau.real,tau.imag],M=M,flux=M*M,Q=Q,K=K,
                lambda_value=lam,max_growth=max_growth,min_nonzero_frequency=min_frequency,
                max_norm_error=max_norm_error,label='NUMERICAL')


def main():
    p=argparse.ArgumentParser();p.add_argument('--M',type=int,default=256)
    a=p.parse_args()
    rows=[]
    for tau in [1j,-0.5+np.sqrt(3)/2*1j,-0.4+0.95j,0.4+1j,0.5+1.052j]:
        row=check(tau,a.M,16,8);rows.append(row)
        print(json.dumps(row),flush=True)
    # Known unstable controls must remain detectable at the enlarged torus.
    assert rows[0]['max_growth']>0.07
    assert rows[4]['max_growth']>0.00009
    for row in rows[1:4]: assert row['max_growth']<1e-9
    for row in rows: assert row['max_norm_error']<1e-10
    print('exit=0')


if __name__=='__main__':main()
