"""Plot an `oscnet --out run.csv` file: order parameter R(t) and x_i(t) of all oscillators as an image."""

import argparse
import csv

import matplotlib.pyplot as plt


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--out", default="docs/run.png")
    ap.add_argument("--title", default="")
    args = ap.parse_args()

    with open(args.csv) as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = [[float(v) for v in r] for r in reader]
    xcols = [i for i, h in enumerate(header) if h.startswith("x")]
    t = [r[0] for r in rows]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True, gridspec_kw={"height_ratios": [1, 2]})
    a1.plot(t, [r[1] for r in rows])
    a1.set_ylabel("R")
    a1.set_ylim(0, 1.02)
    a1.grid(alpha=0.3)
    if args.title:
        a1.set_title(args.title)
    img = [[r[c] for r in rows] for c in xcols]
    im = a2.imshow(
        img, aspect="auto", origin="lower", extent=(t[0], t[-1], 0, len(xcols)), cmap="RdBu_r", vmin=-1, vmax=1
    )
    a2.set_xlabel("t")
    a2.set_ylabel("oscillator")
    fig.colorbar(im, ax=a2, label="x_i", pad=0.01)
    fig.tight_layout()
    fig.savefig(args.out, dpi=110)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
