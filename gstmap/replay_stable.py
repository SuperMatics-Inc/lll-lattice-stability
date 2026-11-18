"""Check all stable leaves in the supplied tree.
The arb mode checks the near and far interval bounds. The sample mode uses a separate floating-point calculation.
"""
import json,sys,time
from fractions import Fraction as F
from concurrent.futures import ProcessPoolExecutor


def job(arg):
    e,mode=arg
    if mode=='arb':
        from boxcert import try_stable
        ok,info=try_stable(F(e['c1']),F(e['c2']),F(e['w']),N=e['info']['N'],rmax=F(e['info']['rmax']))
        return dict(key=[e[k] for k in ('c1','c2','w')],ok=bool(ok),info=info)
    from verify_tree import job_s_replay
    r=job_s_replay(e)
    return dict(key=[e[k] for k in ('c1','c2','w')],ok=bool(r[5]),sample_minimum=r[4])


def main():
    from verify_tree import load
    leaves=[e for e in load(sys.argv[1]) if e['ev']=='S']
    assert len(leaves)==4178
    nw=int(sys.argv[2]); assert 1<=nw<=6
    mode=sys.argv[3]; assert mode in ('arb','sample')
    start=time.monotonic(); failures=0; minimum=float('inf')
    with ProcessPoolExecutor(nw) as pool:
        for row in pool.map(job,[(e,mode) for e in leaves]):
            print(json.dumps(row),flush=True)
            failures+=int(not row['ok'])
            minimum=min(minimum,row.get('sample_minimum',float('inf')))
    print(json.dumps(dict(summary=True,mode=mode,leaves=len(leaves),failures=failures,seconds=time.monotonic()-start,sampled_minimum=minimum if mode=='sample' else None)),flush=True)
    return int(failures>0)

if __name__=='__main__': sys.exit(main())
