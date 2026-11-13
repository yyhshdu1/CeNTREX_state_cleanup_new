"""Forward-velocity bookkeeping for the cached simulation results under tmp/.

The models run at the measured forward velocity VZ_MEASURED = 156 m/s (TOF 153.6-156.8 m/s). Results computed at
another vz live in a subfolder vz<vz>/ of the same cache directory, with unchanged file names; everything computed
before 2026-11-09 was at 184 m/s and is in vz184/. A subfolder (not a file-name suffix) keeps the 184 m/s files out
of every "{stem}_*.npz" glob, so they can never be mixed into a 156 m/s result.

Writers: use vz_dir(root, vz) as the output directory (root itself at 156 m/s).
Readers of results that exist only at 184 m/s: find()/glob() fall back to vz184/ and say so on stderr.
"""
import sys
from pathlib import Path

VZ_MEASURED = 156.0
OLD = "vz184"
_warned = set()


def vz_dir(root, vz):
    root = Path(root)
    return root if float(vz) == VZ_MEASURED else root / f"vz{float(vz):g}"


def _warn(root):
    if root not in _warned:
        _warned.add(root)
        print(f"[vz_paths] no 156 m/s result in {root}; using the 184 m/s one in {root / OLD}", file=sys.stderr)


def find(path, root=None):
    """path if it exists, else the same file under <root>/vz184/ (root defaults to path's directory)."""
    path = Path(path)
    if path.exists():
        return path
    root = Path(root) if root is not None else path.parent
    old = root / OLD / path.relative_to(root)
    if old.exists():
        _warn(root)
        return old
    return path


def glob(root, pattern):
    """sorted root.glob(pattern), or the matches under <root>/vz184/ if there are none."""
    root = Path(root)
    hits = sorted(root.glob(pattern))
    if hits:
        return hits
    hits = sorted((root / OLD).glob(pattern))
    if hits:
        _warn(root)
    return hits
