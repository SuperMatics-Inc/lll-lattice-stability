"""Check that the interval verifier accepts valid witnesses and rejects altered ones."""
import json,sys
from fractions import Fraction as F
import cancel_iv as C

if __name__=='__main__':
    q=['1/64','-1/128']
    cases=[('unstable control',F(53,50),F(1,1000),q,4,True),
           ('wrong-region stable box',F(26,25),F(1,1000),q,4,False),
           ('quartic-curve straddling box',F(10491,10000),F(1,1000),q,4,False),
           ('zero Bloch witness',F(53,50),F(1,1000),['0','0'],4,False),
           ('M=0 inadequate truncation',F(53,50),F(1,1000),q,0,False)]
    fails=0
    for name,y,w,q,M,want in cases:
        up=C.upper(F(1,2),y,w,q,M=M)
        accepted=bool(up<0); fails+=accepted!=want
        print(json.dumps(dict(case=name,upper=float(up),accepted=accepted,expected=want)),flush=True)
    print(json.dumps(dict(summary=True,n=len(cases),fails=fails)),flush=True)
    sys.exit(bool(fails))
