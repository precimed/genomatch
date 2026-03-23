#!/usr/bin/env python3
from ._cli_utils import run_cli
from .apply_vmap_bfile import main


def cli_main() -> int:
    return run_cli(main)


if __name__ == "__main__":
    raise SystemExit(cli_main())
