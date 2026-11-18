"""Replay the complete reproducibility supplement (Python checks; Lean is separate).
Install uv first; this command fetches public numerical packages into uv's cache.
Outputs are generated calculations, not experimental measurements.
"""
from pathlib import Path
import argparse
import json
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / ('cert' if (ROOT/'cert').is_dir() else 'certificates')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=ROOT/'replay-output')
    a = p.parse_args()
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    runs = [
        ('arb-and-numerics', [sys.executable, str(ROOT/'scripts/reproduce.py'), '--output', str(out/'original')]),
        ('independent-intervals', ['uv','run','--no-project','--with-requirements',str(ROOT/'requirements.txt'),'python',str(ROOT/'scripts/rigor_replay.py'),'--output',str(out/'independent')]),
        ('exact-leaf-trace', ['uv','run','--no-project','--with-requirements',str(ROOT/'requirements.txt'),'python',str(ROOT/'scripts/trace_certificates.py'),'--output',str(out/'trace')]),
        ('finite-torus-256', ['uv','run','--no-project','--with-requirements',str(ROOT/'requirements.txt'),'python',str(CERT/'finite_torus.py'),'--M','256'])]
    records = []
    for name, command in runs:
        start = time.monotonic()
        with (out/f'{name}.log').open('w', encoding='utf-8') as log:
            log.write('command: '+ ' '.join(command)+'\n'); log.flush()
            r = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, timeout=3600)
            log.write(f'\nexit={r.returncode}\n')
        records.append(dict(name=name, command=command, exit=r.returncode, seconds=time.monotonic()-start))
        print(name, 'exit=', r.returncode, flush=True)
    (out/'verification.json').write_text(json.dumps(records, indent=2)+'\n')
    if any(x['exit'] for x in records): raise SystemExit(1)


if __name__ == '__main__': main()
