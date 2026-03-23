from __future__ import annotations

import sys
from typing import Callable


def run_cli(main_fn: Callable[[], int | None]) -> int:
    try:
        result = main_fn()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0 if result is None else int(result)
