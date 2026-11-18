"""Verify every released file against SHA256SUMS (LF checkout bytes)."""
from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parents[1]
if __name__ == '__main__':
    count = 0
    for line in (ROOT / 'SHA256SUMS').read_text(encoding='utf-8').splitlines():
        expected, name = line.split('  ', 1)
        path = ROOT / name
        if not path.is_relative_to(ROOT) or '..' in Path(name).parts:
            raise ValueError('Unsafe manifest path')
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f'Hash mismatch: {name}')
        count += 1
    print(f'Verified {count} LF-checkout file hashes')
