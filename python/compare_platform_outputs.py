"""
Numerically compares two GWASH/LDSC sim-output CSVs produced by the same driver script
and arguments (same seed) on different machines/platforms, to check that a run moved to a
different machine (e.g. Windows -> Unix) reproduces the same results up to floating-point
noise from differing BLAS/LAPACK backends.

Usage: python compare_platform_outputs.py <csv_a> <csv_b>
Exits 0 if the two files match (same header/row count, max relative diff < 1e-6), else 1.
"""
import csv
import sys

TOL_REL = 1e-6


def load(path):
    with open(path, newline="") as f:
        return list(csv.reader(f))


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <csv_a> <csv_b>", file=sys.stderr)
        sys.exit(2)

    path_a, path_b = sys.argv[1], sys.argv[2]
    a, b = load(path_a), load(path_b)

    print(f"a: {path_a} ({len(a)} rows)")
    print(f"b: {path_b} ({len(b)} rows)")

    header_match = a[0] == b[0]
    print(f"header match: {header_match}")

    if len(a) != len(b):
        print(f"ROW COUNT MISMATCH: {len(a)} vs {len(b)}")
        sys.exit(1)

    max_abs_diff = 0.0
    max_rel_diff = 0.0
    max_rel_loc = None
    nrows = len(a)
    ncols = len(a[0])

    for r in range(1, nrows):
        for c in range(ncols):
            va, vb = float(a[r][c]), float(b[r][c])
            d = abs(va - vb)
            rel = d / (abs(va) + 1e-12)
            if d > max_abs_diff:
                max_abs_diff = d
            if rel > max_rel_diff:
                max_rel_diff = rel
                max_rel_loc = (r, a[0][c])

    print(f"max abs diff: {max_abs_diff:.3e}")
    print(f"max rel diff: {max_rel_diff:.3e}" + (f"  (row {max_rel_loc[0]}, col {max_rel_loc[1]})" if max_rel_loc else ""))

    if header_match and max_rel_diff < TOL_REL:
        print(f"PASS: outputs match within tolerance (rel diff < {TOL_REL:.0e})")
        sys.exit(0)
    else:
        print(f"FAIL: outputs diverge beyond tolerance (rel diff >= {TOL_REL:.0e})")
        sys.exit(1)


if __name__ == "__main__":
    main()
