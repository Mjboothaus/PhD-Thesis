"""Command-line interface: python -m charged_fluids.bulk.cli [SYMBOL | --list | --all]."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from charged_fluids.bulk.generate import (
    c_short_tail_max,
    generate_all_bulk_cr,
    generate_bulk_cr,
    list_bulk_symbols,
    wall_usable_c_short,
)


def _print_result(symbol: str, result, cr_path: Path | None, meta_path: Path | None) -> None:
    status = "converged" if result.converged else "NOT converged"
    closure = result.meta.get("closure", result.meta.get("fluid", "?"))
    print(
        f"{symbol}: {status} after {result.n_iter} iterations "
        f"(residual={result.residual:.3e}, source={closure})"
    )
    tail = c_short_tail_max(result)
    print(
        f"  |c_short| tail max = {tail:.6e}  n_pair={result.c_short.shape[0]}  n_r={result.r.size}"
    )
    note = result.meta.get("note")
    if note:
        print(f"  note: {note}")
    if cr_path is not None:
        print(f"  wrote {cr_path}")
        print(f"  meta  {meta_path}")


def _exit_code(symbol: str, result) -> int:
    if result.converged:
        return 0
    if wall_usable_c_short(result):
        print("  note: residual not tight, but c_short tail OK for wall-fluid use")
        return 0
    return 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate bulk direct correlation c(r) artefacts for fluid symbols."
    )
    parser.add_argument(
        "symbol",
        nargs="?",
        default=None,
        help="Fluid symbol (e.g. lj1, kcl, 1_1). Omit with --list / --all.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List symbols that support bulk generation and exit.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Generate bulk artefacts for every supported fluid symbol.",
    )
    parser.add_argument(
        "--n-points",
        type=int,
        default=None,
        help="Radial grid size (default: fluid [symbol.bulk] TOML or 1024)",
    )
    parser.add_argument(
        "--delta-r",
        type=float,
        default=None,
        help="Radial spacing in A (default: fluid TOML or 0.05)",
    )
    parser.add_argument(
        "--closure",
        choices=("hnc", "py"),
        default=None,
        help="OZ closure (default: fluid TOML)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for .cr.dat / .meta.json (default: data/bulk)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Solve but do not write data/bulk artefacts",
    )
    args = parser.parse_args(argv)

    if args.list:
        symbols = list_bulk_symbols()
        for s in symbols:
            print(s)
        return 0 if symbols else 1

    if args.all:
        results = generate_all_bulk_cr(
            bulk_dir=args.output_dir,
            n_points=args.n_points,
            delta_r=args.delta_r,
            closure=args.closure,
            write=not args.dry_run,
        )
        exit_code = 0
        for symbol, result, cr_path, meta_path in results:
            _print_result(symbol, result, cr_path, meta_path)
            code = _exit_code(symbol, result)
            if code != 0:
                exit_code = code
        return exit_code

    if not args.symbol:
        parser.error("symbol is required unless --list or --all is given")

    result, cr_path, meta_path = generate_bulk_cr(
        args.symbol,
        bulk_dir=args.output_dir,
        n_points=args.n_points,
        delta_r=args.delta_r,
        closure=args.closure,
        write=not args.dry_run,
    )
    _print_result(args.symbol, result, cr_path, meta_path)
    return _exit_code(args.symbol, result)


if __name__ == "__main__":
    sys.exit(main())
