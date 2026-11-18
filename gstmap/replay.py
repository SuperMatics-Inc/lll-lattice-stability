"""Run the shape-map checks using the original Arb code for the stable boxes.
uv run --no-project --with python-flint --with mpmath --with numpy --with scipy python replay.py --output OUT --workers 6
--evidence omits the full Arb check of the stable boxes.
--ci-subset checks the partition, Gamma4 patches and joins, negative controls and cusp.
It also omits the full unstable cover and the stable-box sampling.
Build the separate Lean project with: cd lean && lake build.
"""
import argparse, hashlib, json, os, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True,help='Directory for results and logs')
    ap.add_argument('--workers',type=int,default=6,help='Number of worker processes (default: 6)')
    ap.add_argument('--evidence',action='store_true',help='Omit the full Arb check of the stable boxes')
    ap.add_argument('--ci-subset',action='store_true',help='Run the smaller set of checks used by GitHub Actions')
    a=ap.parse_args(); assert 1<=a.workers<=6
    out=a.output.resolve(); out.mkdir(parents=True,exist_ok=True)
    tree=ROOT/'data/tree_smallq.jsonl'
    if not tree.exists(): tree=ROOT/'tree_smallq.jsonl'
    assert digest(tree)=='647585a453c630987b8ee98a03f2af03065d6e6598636f11e4dd2f18e49012ee'
    env=dict(os.environ,SM_OUT=str(out)); nw=str(a.workers)
    commands=[('partition',['gen_data.py',str(tree),str(out/'Data.lean')]),
              ('independent-U',['cancel_iv.py',str(tree),nw]),
              ('gamma-implicit',['gamma_implicit.py',nw]),
              ('gamma-chain',['verify_gamma_chain.py',str(out/'gamma_implicit.jsonl'),nw]),
              ('U-mutants',['cancel_iv_mutants.py']),
              ('gamma-mutants',['gamma_mutants.py',str(out/'gamma_implicit.jsonl')]),
              ('sample-S',['replay_stable.py',str(tree),nw,'sample']),
              ('cusp',['cusp_arb.py','3'])]
    skipped=[]
    if a.ci_subset:
        commands=[(label,args) for label,args in commands if label not in ('independent-U','sample-S')]
        skipped.extend(['independent U cover','independent S sampling'])
    if a.evidence or a.ci_subset:
        skipped.append('primary Arb S replay')
    else:
        commands.append(('primary-Arb-S',['replay_stable.py',str(tree),nw,'arb']))
    receipts=[]
    for label,args in commands:
        start=time.monotonic(); command=[sys.executable,*args]
        log=out/(label+'.log')
        with log.open('w') as f:
            p=subprocess.run(command,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
        row=dict(label=label,command=command,exit_code=p.returncode,seconds=time.monotonic()-start,log=log.name,sha256=digest(log))
        receipts.append(row); print(json.dumps(row),flush=True)
    expected=ROOT/'lean/GSTMap/Data.lean'
    # The release checkout and generated exact data must both be LF bytes.
    data_equal=(out/'Data.lean').read_bytes()==expected.read_bytes()
    receipt=dict(input_sha256=digest(tree),checks=receipts,lean_data_sha256=digest(out/'Data.lean'),lean_data_matches=data_equal,skipped=skipped,output_hashes={p.name:digest(p) for p in out.iterdir() if p.is_file() and p.name!='replay-receipt.json'})
    receipt['lean_data_comparison']='Exact LF source bytes'
    receipt['exit_code']=int(not data_equal or any(r['exit_code'] for r in receipts))
    (out/'replay-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt['exit_code']

if __name__=='__main__': sys.exit(main())
