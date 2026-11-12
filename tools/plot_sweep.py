"""Plot one or more `oscnet --sweep` outputs: mean order parameter and frequency spread vs the swept value."""

import argparse
import csv

import matplotlib.pyplot as plt


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="+", help="label=path.csv or path.csv")
    ap.add_argument("--out", default="docs/sweep.png")
    args = ap.parse_args()

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    for spec in args.csv:
        label, path = spec.rsplit("=", 1) if "=" in spec else (spec, spec)
        with open(path) as f:
            rows = list(csv.DictReader(f))
        key = list(rows[0])[0]
        v = [float(r[key]) for r in rows]
        a1.plot(v, [float(r["R_mean"]) for r in rows], "o-", label=label)
        a2.plot(v, [float(r["frequency_std"]) for r in rows], "o-", label=label)
        a1.set_xlabel(key)
        a2.set_xlabel(key)
    a1.set_ylabel("mean order parameter R (second half)")
    a2.set_ylabel("std of mean frequencies")
    a2.set_yscale("log")
    for a in (a1, a2):
        a.grid(alpha=0.3)
        a.legend()
    fig.tight_layout()
    fig.savefig(args.out, dpi=110)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
