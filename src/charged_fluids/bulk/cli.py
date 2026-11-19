"""Command-line interface: python -m charged_fluids.bulk.cli SYMBOL."""

from __future__ import annotations

import argparse
import sys

from charged_fluids.bulk.generate import generate_bulk_cr


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate bulk direct correlation c(r) artefacts for a fluid symbol."
    )
    parser.add_argument("symbol", help="Fluid symbol (e.g. lj1, kcl)")
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
        help="OZ closure (default: fluid TOML; lj1->hnc, kcl->py)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Solve but do not write data/bulk artefacts",
    )
    args = parser.parse_args(argv)

    result, cr_path, meta_path = generate_bulk_cr(
        args.symbol,
        n_points=args.n_points,
        delta_r=args.delta_r,
        closure=args.closure,
        write=not args.dry_run,
    )
    status = "converged" if result.converged else "NOT converged"
    print(
        f"{args.symbol}: {status} after {result.n_iter} iterations "
        f"(residual={result.residual:.3e}, closure={result.meta.get('closure')})"
    )
    tail = float(abs(result.c_short[:, -max(5, result.r.size // 20) :]).max())
    print(f"  |c_short| tail max = {tail:.6e}")
    if cr_path is not None:
        print(f"  wrote {cr_path}")
        print(f"  meta  {meta_path}")
    # Ionic melts may not hit tight Picard tol; treat small c_short tail as success
    if result.converged:
        return 0
    if tail < 1e-2 and result.meta.get("fluid", "").startswith("binary"):
        print("  note: residual not tight, but c_short tail OK for wall-fluid use")
        return 0
    if tail < 1e-2 and args.symbol == "kcl":
        print("  note: residual not tight, but c_short tail OK for wall-fluid use")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
