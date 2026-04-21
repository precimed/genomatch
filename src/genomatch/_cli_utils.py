from __future__ import annotations

import logging
import sys
from typing import Callable


LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def configure_cli_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format=LOG_FORMAT,
        datefmt=LOG_DATE_FORMAT,
        stream=sys.stderr,
        force=True,
    )


def run_cli(main_fn: Callable[[], int | None]) -> int:
    configure_cli_logging()
    try:
        result = main_fn()
    except Exception as exc:
        logging.error("Error: %s", exc)
        return 1
    return 0 if result is None else int(result)
