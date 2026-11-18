"""Check three certificates with mpmath.iv and save their output and exit codes.
This compares the arithmetic backends using the same analytic bounds.
Run from any directory.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / ('cert' if (ROOT/'cert').exists() else 'certificates')


def child(name):
    sys.path.insert(0, str(CERT))
    import interval_backend
    import runpy
    sys.modules['flint'] = interval_backend
    if name == 'mu':
        sys.argv = ['cert_mu.py']
        runpy.run_path(str(CERT/'cert_mu.py'), run_name='__main__')
    else:
        args = ['-2/5','19/20','1/100'] if name == 'minus' else ['2/5','1','1/40']
        sys.argv = ['cert_sigma.py', *args, '20','8','192','13']
        runpy.run_path(str(CERT/'cert_sigma.py'), run_name='__main__')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--child', choices=['mu','minus','plus'], help='Check only the selected certificate')
    p.add_argument('--output', type=Path, default=ROOT/'logs'/'rigor', help='Directory for calculation logs')
    a = p.parse_args()
    if a.child:
        child(a.child)
        return
    a.output.mkdir(parents=True, exist_ok=True)
    rows = []
    for name in ['mu','minus','plus']:
        cmd = [sys.executable, str(Path(__file__).resolve()), '--child', name]
        start = time.monotonic()
        r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf8', timeout=18000)
        text = 'command: '+' '.join(cmd)+'\n'+r.stdout+r.stderr+f'\nexit={r.returncode}\n'
        dest = a.output/f'iv-{name}.log'
        dest.write_text(text, encoding='utf8')
        row = dict(name=name,exit=r.returncode,seconds=time.monotonic()-start,
                   log=dest.name,sha256=hashlib.sha256(dest.read_bytes()).hexdigest())
        rows.append(row)
        print(json.dumps(row), flush=True)
    (a.output/'iv-replay.json').write_text(json.dumps(rows,indent=2)+'\n')
    if any(r['exit'] for r in rows): raise SystemExit(1)


if __name__ == '__main__': main()
