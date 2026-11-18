"""Run the numerical programs and record commands, exit codes and checksums.

Run from any directory with Python and uv installed:
  python scripts/reproduce.py --output /path/to/logs
The calculations are mathematical/numerical, not experimental measurements.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import subprocess
import time

PAPER = Path(__file__).resolve().parents[1]
CERT = PAPER / ('cert' if (PAPER / 'cert').is_dir() else 'certificates')
PREFIX = ['uv', 'run', '--no-project', '--with-requirements', str(PAPER/'requirements.txt'), 'python']
RUNS = [
    ('cert_mu', ['cert_mu.py']),
    ('closed_form', ['closed_form.py']),
    ('quartic_decomp', ['quartic_decomp.py']),
    ('taylor_quartic', ['taylor_quartic.py']),
    ('symbol_checks', ['symbol_checks.py']),
    ('examples', ['examples.py']),
    ('cert_sigma_m04_095', ['cert_sigma.py', '-2/5', '19/20', '1/100', '20', '8', '192', '13']),
    ('cert_sigma_p04_100', ['cert_sigma.py', '2/5', '1', '1/40', '20', '8', '192', '13']),
    ('pocket_scan', ['pocket_scan.py']),
    ('landscape', ['landscape.py']),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='Directory for calculation logs')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    records = []
    versions = PREFIX + ['-c', 'import sys,numpy,scipy,mpmath,flint,sympy; '
                         'print("Python",sys.version); '
                         'print("numpy",numpy.__version__); '
                         'print("scipy",scipy.__version__); '
                         'print("mpmath",mpmath.__version__); '
                         'print("python-flint",flint.__version__); '
                         'print("sympy",sympy.__version__)']
    for name, command in [('versions', versions)] + [(name, PREFIX + argv) for name, argv in RUNS]:
        started = time.monotonic()
        result = subprocess.run(command, cwd=CERT, capture_output=True, text=True,
                                encoding='utf-8', errors='replace', timeout=1200)
        log = output / f'{name}.log'
        with log.open('w', encoding='utf-8', newline='\n') as stream:
            stream.write('command: ' + ' '.join(command) + '\n'
                         + result.stdout + result.stderr + f'\nexit={result.returncode}\n')
        records.append({'name': name, 'command': command, 'exit_code': result.returncode,
                        'elapsed_seconds': time.monotonic() - started,
                        'log': log.name, 'sha256': hashlib.sha256(log.read_bytes()).hexdigest()})
        print(f'{name}: exit={result.returncode}', flush=True)
    manifest = {'platform': platform.platform(),
                'sources': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in sorted(CERT.glob('*.py'))},
                'runs': records}
    (output / 'replay.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    if any(row['exit_code'] != 0 for row in records):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
