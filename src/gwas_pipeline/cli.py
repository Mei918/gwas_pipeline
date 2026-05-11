from __future__ import annotations

import argparse
import sys
from typing import List, Optional

from gwas_pipeline.doctor import run_doctor
from gwas_pipeline.steps import STEP_MODULES, get_step_main


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gwas-pipeline",
        description="Run the packaged GWAS tutorial workflow steps.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("doctor", help="Check GWAS server dependency availability.")

    for step_name in STEP_MODULES:
        step_parser = subparsers.add_parser(
            step_name,
            help=f"Run GWAS tutorial {step_name}.",
            description=f"Delegate to the packaged GWAS tutorial {step_name} module.",
        )
        step_parser.add_argument(
            "step_args",
            nargs=argparse.REMAINDER,
            help="Arguments forwarded to the underlying step module.",
        )

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    raw_argv = list(sys.argv[1:] if argv is None else argv)

    if raw_argv and raw_argv[0] in STEP_MODULES:
        return get_step_main(raw_argv[0])(raw_argv[1:])
    if raw_argv and raw_argv[0] == "doctor":
        return run_doctor(raw_argv[1:])

    parser = build_parser()
    parser.parse_args(raw_argv)
    parser.print_help(sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
