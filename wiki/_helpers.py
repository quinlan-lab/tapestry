"""
Shared helpers for tapestry's wiki generator.

Permalink utilities (`head_sha`, `permalink`) — pin GitHub source links
to the current HEAD SHA at generate time.
"""

from __future__ import annotations

import subprocess

TAPESTRY_REPO = "quinlan-lab/tapestry"


def head_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip()
    except Exception:
        return "main"


def permalink(path: str, line: int, sha: str, repo: str = TAPESTRY_REPO) -> str:
    return f"https://github.com/{repo}/blob/{sha}/{path}#L{line}"
