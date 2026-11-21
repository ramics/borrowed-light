"""Run manifests: if a result cannot be regenerated from its manifest alone,
it does not count (SPEC.md)."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

PACKAGES = ("numpy", "scipy", "pandas", "scikit-learn", "pyyaml", "pyreadr")


def _git(*args) -> str | None:
    try:
        return subprocess.check_output(["git", *args], cwd=Path(__file__).parent, text=True,
                                       stderr=subprocess.DEVNULL).strip()
    except Exception:
        return None


def file_digest(path: Path) -> str:
    """SHA-256 of a file, or of a directory's files (relative names + bytes,
    in sorted order)."""
    path = Path(path)
    if path.is_dir():
        h = hashlib.sha256()
        for f in sorted(p for p in path.rglob("*") if p.is_file()):
            h.update(str(f.relative_to(path)).encode())
            h.update(bytes.fromhex(file_digest(f)))
        return h.hexdigest()
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _redact(argv: list[str]) -> list[str]:
    """Drop the data location from the recorded command line."""
    out, skip = [], False
    for a in argv:
        if skip:
            out.append("<data>")
            skip = False
        elif a == "--data":
            out.append(a)
            skip = True
        elif a.startswith("--data="):
            out.append("--data=<data>")
        else:
            out.append(a)
    return out


def build(config: dict, data_path: Path | None, argv: list[str]) -> dict:
    versions = {}
    for p in PACKAGES:
        try:
            versions[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            pass
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "config": config,
        "code": {"commit": _git("rev-parse", "HEAD"),
                 "dirty": bool(_git("status", "--porcelain")),
                 "argv": _redact(list(argv))},
        # The data's fingerprint, never its location or contents.
        "data": {"sha256": file_digest(data_path) if data_path else None,
                 "filename": data_path.name if data_path else None},
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": versions,
    }


def write(manifest: dict, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
