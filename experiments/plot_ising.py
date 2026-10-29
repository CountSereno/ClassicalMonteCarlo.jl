"""Plot recorded runs; optionally gzip the raw CSV traces losslessly.

python experiments/plot_ising.py experiments/results/demo --compress
Requires NumPy and Matplotlib. No simulation is performed by this script.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import platform
from pathlib import Path
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

TC = 2 / np.log(1 + np.sqrt(2))
COLORS = ["#2563eb", "#e87924", "#0d9488", "#9333ea", "#dc2626", "#64748b"]


def load_rows(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    text = {"run_id", "experiment", "algorithm", "initial", "window_truncated", "block_warning"}
    return [{k: (v if k in text else float(v)) for k, v in row.items()} for row in rows]


def trace(root, row):
    path = root / "traces" / f"{row['run_id']}.csv"
    if path.exists():
        return np.genfromtxt(path, delimiter=",", names=True)
    with gzip.open(path.with_suffix(".csv.gz"), "rt") as stream:
        return np.genfromtxt(stream, delimiter=",", names=True)


def acf(x, maxlag=200):
    x = np.asarray(x) - np.mean(x)
    n = len(x)
    f = np.fft.rfft(x, n=2 * n)
    cov = np.fft.irfft(f * f.conj())[: maxlag + 1]
    return cov / cov[0]


def combine(rows, key):
    """Descriptive mean and conservative diagnostic SE, not a coverage guarantee."""
    values = np.array([r[key] for r in rows])
    within = np.sqrt(sum(max(r[key + "_stderr"], r[key + "_stderr_large_block"]) ** 2 for r in rows)) / len(rows)
    between = values.std(ddof=1) / np.sqrt(len(rows))
    return values.mean(), max(within, between)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--compress", action="store_true")
    args = parser.parse_args()
    root = args.directory
    rows = load_rows(root / "summary.csv")
    figdir = root / "figures"
    figdir.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "figure.facecolor": "white", "savefig.facecolor": "white", "axes.titleweight": "bold"})
    bench = [r for r in rows if r["experiment"] == "benchmark"]
    sizes = sorted({int(r["L"]) for r in rows})
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8), constrained_layout=True)
    for algorithm, color, marker in [("metropolis", COLORS[0], "o"), ("wolff", COLORS[1], "s")]:
        selected = [r for r in bench if r["algorithm"] == algorithm]
        for ax, key in zip(axs, ["tau_E", "tau_E_spin_sweeps", "ESS_E_per_second"]):
            for r in selected:
                ax.scatter(r["L"], r[key], color=color, marker=marker, alpha=0.5, s=25)
            means = [np.mean([r[key] for r in selected if r["L"] == L]) for L in sizes]
            ax.plot(sizes, means, color=color, marker=marker, label=algorithm.title())
            ax.set(xlabel="Linear size L", yscale="log", xticks=sizes)
            ax.grid(alpha=0.18)
    axs[0].set(title="Energy autocorrelation", ylabel=r"$\tau_{int}$ in recorded updates")
    axs[1].set(title="Approximate spin-work units", ylabel=r"$\tau_{int}$ in spin sweeps")
    axs[2].set(title="Effective samples / elapsed time", ylabel="Estimated energy ESS / second")
    axs[0].legend(frameon=False)
    fig.suptitle(r"Critical Ising: $T_c=2/\log(1+\sqrt{2})$; independent random/up starts", fontsize=13)
    fig.savefig(figdir / "algorithm_comparison.png", dpi=170)
    plt.close(fig)

    fig, axs = plt.subplots(2, 2, figsize=(12, 6.5), constrained_layout=True)
    Lmax = max(sizes)
    for algorithm, color in [("metropolis", COLORS[0]), ("wolff", COLORS[1])]:
        selected = [r for r in bench if r["algorithm"] == algorithm and r["L"] == Lmax]
        for r in selected:
            t = trace(root, r)
            style = "-" if r["initial"] == "random" else "--"
            label = f"{algorithm}, {r['initial']} start"
            axs[0, 0].plot(t["measurement"][::20], t["total_energy"][::20]/Lmax**2, style, lw=0.6, alpha=0.65, color=color, label=label)
            axs[0, 1].plot(acf(t["total_energy"]), style, color=color, label=label)
            axs[1, 0].plot(t["measurement"][::20], t["total_magnetization"][::20]/Lmax**2, style, lw=0.65, alpha=0.6, color=color)
            axs[1, 1].plot(acf(np.abs(t["total_magnetization"])), style, color=color)
    axs[0, 0].set(title=f"Post-burn energy trace, L={Lmax}", xlabel="Recorded update", ylabel="E/N")
    axs[0, 1].set(title="Energy autocorrelation", xlabel="Lag in recorded updates", ylabel=r"$\rho_E(t)$")
    axs[1, 0].set(title="Signed magnetization trace", xlabel="Recorded update", ylabel="M/N")
    axs[1, 1].set(title="Absolute-magnetization autocorrelation", xlabel="Lag in recorded updates", ylabel=r"$\rho_{|M|}(t)$")
    axs[0, 1].legend(frameon=False, fontsize=8)
    for ax in (axs[0, 1], axs[1, 1]):
        ax.axhline(0, color="#999", lw=0.7)
    fig.savefig(figdir / "chain_diagnostics.png", dpi=170)
    plt.close(fig)

    scan = [r for r in rows if r["experiment"] == "scan"]
    fig, axs = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    columns = ["binder", "susceptibility", "heat_capacity", "abs_magnetization_density"]
    aggregate = []
    for L, color in zip(sizes, COLORS):
        temps = sorted({r["T"] for r in scan if r["L"] == L})
        stats = {}
        for key in columns:
            values = [combine([r for r in scan if r["L"] == L and r["T"] == T], key) for T in temps]
            stats[key] = np.array(values)
        for j, T in enumerate(temps):
            item = {"L": L, "T": T}
            for key in columns:
                item[key], item[key + "_diagnostic_se"] = map(float, stats[key][j])
            aggregate.append(item)
        for ax, key in zip(axs.flat, columns):
            ax.errorbar(temps, stats[key][:, 0], yerr=stats[key][:, 1], fmt="o-", lw=1.1, ms=4, capsize=2, color=color, label=f"L={L}")
            ax.axvline(TC, color="#aaa", ls=":", lw=0.8)
            ax.set_xlabel("Temperature T")
        scaled_x = (np.array(temps) - TC) / TC * L
        for ax, key, factor in [(axs[1, 1], "abs_magnetization_density", L**(1/8)),
                                (axs[1, 2], "susceptibility", L**(-7/4))]:
            ax.errorbar(scaled_x, stats[key][:, 0]*factor, yerr=stats[key][:, 1]*factor, fmt="o-", lw=1.1, ms=4, capsize=2, color=color)
            ax.set_xlabel(r"$[(T-T_c)/T_c] L^{1/\nu}$, known $\nu=1$")
    for ax, title, ylabel in zip(axs.flat, ["Binder crossing diagnostic", "Signed-M susceptibility", "Specific heat", "Absolute magnetization", "Known-exponent magnetization collapse", "Known-exponent susceptibility collapse"],
                                [r"$U_4$", r"$\chi=\beta(\langle M^2\rangle-\langle M\rangle^2)/N$", r"$C_V/N$", r"$\langle |M|\rangle/N$", r"$m L^{1/8}$", r"$\chi L^{-7/4}$"]):
        ax.set(title=title, ylabel=ylabel)
        ax.grid(alpha=0.14)
    axs[0, 0].legend(frameon=False)
    fig.suptitle("Finite-size diagnostics — reference exponents inserted, no exponent fitting", fontsize=13)
    fig.savefig(figdir / "finite_size_diagnostics.png", dpi=170)
    plt.close(fig)
    with (root / "aggregate.json").open("w") as stream:
        json.dump({"error_rule": "max(quadrature of within-chain SE/nchains, between-chain SE); within-chain SE=max(block,4*block). Diagnostic only; two starts do not establish equilibration.",
                   "values": aggregate}, stream, indent=2)

    warnings = [{k: r[k] for k in ["run_id", "experiment", "L", "T", "algorithm", "initial", "tau_E", "tau_absM", "block_warning", "window_truncated"]}
                for r in rows if r["window_truncated"] == "true" or r["block_warning"] == "true"]
    report = {"runs": len(rows), "measurements": int(sum(r["samples"] for r in rows)),
              "production_seconds": sum(r["production_seconds"] for r in rows),
              "warnings": warnings, "python_version": platform.python_version(),
              "numpy_version": np.__version__, "matplotlib_version": matplotlib.__version__}
    with (root / "diagnostics.json").open("w") as stream:
        json.dump(report, stream, indent=2)
    if args.compress:
        for path in sorted((root / "traces").glob("*.csv")):
            with path.open("rb") as source, path.with_suffix(".csv.gz").open("wb") as raw:
                with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as target:
                    shutil.copyfileobj(source, target)
            path.unlink()
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
