"""Verify saved trace integrity, Ising constraints and reported raw moments.

python experiments/verify_results.py experiments/results/demo
Requires NumPy and Python 3.11+; verification does not establish equilibration.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from pathlib import Path
import sys
import tomllib

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument(
        "--allow-source-drift", action="store_true",
        help="inspect historical data with different source hashes; retain all data checks",
    )
    args = parser.parse_args()
    root = args.directory.resolve()
    with (root / "metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    with (root / "summary.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == metadata["runs"]
    assert len({r["seed"] for r in rows}) == len(rows), "reused seed"
    package = Path(__file__).resolve().parents[1]
    source_mismatches = []
    for name, digest in metadata["source_sha256"].items():
        source = package / name if "/" in name else package / "src" / name
        if not source.is_file() or hashlib.sha256(source.read_bytes()).hexdigest() != digest:
            source_mismatches.append(name)
    if source_mismatches:
        message = "Current source differs from recorded source: " + ", ".join(source_mismatches)
        if not args.allow_source_drift:
            parser.error(message + ". Use --allow-source-drift only to inspect historical data.")
        print("Warning: " + message + ". Checking historical data without claiming current-source reproduction.", file=sys.stderr)
    count = 0
    for row in rows:
        path = root / "traces" / f"{row['run_id']}.csv"
        if path.exists():
            data = np.genfromtxt(path, delimiter=",", names=True)
        else:
            with gzip.open(path.with_suffix(".csv.gz"), "rt") as stream:
                data = np.genfromtxt(stream, delimiter=",", names=True)
        n = int(row["samples"])
        N = int(row["L"]) ** 2
        E, M = data["total_energy"], data["total_magnetization"]
        assert len(data) == n
        assert np.array_equal(data["measurement"], np.arange(1, n + 1))
        assert np.all(np.isfinite(E)) and np.all(np.isfinite(M))
        assert np.all(np.abs(E) <= 2 * N) and np.all(np.abs(M) <= N)
        assert np.all((E + 2 * N) % 4 == 0), "invalid Ising energy grid"
        assert np.all((M + N) % 2 == 0), "invalid magnetization parity"
        assert np.all((data["flipped_sites"] >= 0) & (data["flipped_sites"] <= N))
        beta = 1 / float(row["T"])
        m2 = np.mean(M * M)
        expected = {
            "energy_density": E.mean() / N,
            "abs_magnetization_density": np.abs(M).mean() / N,
            "heat_capacity": beta**2 * E.var() / N,
            "susceptibility": beta * M.var() / N,
            "susceptibility_abs": beta * (m2 - np.abs(M).mean()**2) / N,
            "binder": 1 - np.mean(M**4) / (3 * m2**2),
        }
        for key, expected_value in expected.items():
            assert np.isclose(float(row[key]), expected_value, rtol=1e-10, atol=1e-11), (row["run_id"], key)
            assert float(row[key + "_stderr"]) >= 0
            assert float(row[key + "_stderr_large_block"]) >= 0
        assert float(row["production_seconds"]) > 0
        count += n
    print(f"Verified {len(rows)} independent-seed traces, {count:,} measurements, six observables per run and physical bounds/parities.")
    if not source_mismatches:
        print("All recorded source SHA256 hashes match the current files.")


if __name__ == "__main__":
    main()
