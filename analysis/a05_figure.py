"""Build the summary figure.

Run it with:

    python analysis/a05_figure.py
"""

import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ldlib import (
    SAMPLE,
    LimbDarkeningLibrary,
    load_facular_contrast,
    load_sample_uncertainties,
    load_table3,
)

FIG = pathlib.Path(__file__).resolve().parents[1] / "results" / "figures"
STEPS = {"FeH": 0.10, "Teff": 100.0, "logg": 0.10}


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    rows = load_table3()
    kepler = [r for r in rows if r["passband"] == "Kepler"]
    refld = next(r for r in kepler if r["is_refld"])

    fig, axes = plt.subplots(1, 4, figsize=(18, 4.6))

    # 1. Every library against the one that was chosen.
    ax = axes[0]
    order = np.argsort([r["dh1"] for r in kepler])
    labels, vals, errs, is_ref = [], [], [], []
    for i in order:
        r = kepler[i]
        labels.append(r["library"].replace(" (", "\n("))
        vals.append(r["dh1"]); errs.append(r["dh1_err"]); is_ref.append(r["is_refld"])
    colours = ["#C0392B" if f else "#7F8C8D" for f in is_ref]
    ax.barh(range(len(vals)), vals, xerr=errs, color=colours, height=0.7,
            error_kw={"lw": 1, "ecolor": "0.3"})
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.axvline(0, color="k", lw=0.8)
    ax.axvline(refld["dh1"], color="#C0392B", lw=0.8, ls="--")
    ax.set_xlabel("$\\Delta h_1$  (observed $-$ model)")
    ax.set_title("Kepler: the same measurement against\nnine non-magnetic libraries",
                 fontsize=10)
    ax.grid(alpha=0.25, axis="x")

    # 2. Signal against the systematics that could mimic it.
    ax = axes[1]
    d1 = np.array([r["dh1"] for r in kepler])
    lib_set1 = LimbDarkeningLibrary("set1", "Kepler")
    lib_set2 = LimbDarkeningLibrary("set2", "Kepler")
    h1a, _ = lib_set1.h(SAMPLE["FeH"], SAMPLE["Teff"], SAMPLE["logg"])
    h1b, _ = lib_set2.h(SAMPLE["FeH"], SAMPLE["Teff"], SAMPLE["logg"])
    unc = load_sample_uncertainties()
    s1_feh, _ = lib_set1.derivative(SAMPLE["FeH"], SAMPLE["Teff"], SAMPLE["logg"],
                                    "FeH", STEPS["FeH"])
    s1_teff, _ = lib_set1.derivative(SAMPLE["FeH"], SAMPLE["Teff"], SAMPLE["logg"],
                                     "Teff", STEPS["Teff"])
    names = ["signal\n(attributed to\nmagnetism)", "abundance\nscale\n(set 1 vs 2)",
             "[Fe/H] error\n(1 sigma)", "$T_{\\rm eff}$ error\n(1 sigma)",
             "library\nspread"]
    sizes = [abs(refld["dh1"]), abs(h1a - h1b),
             abs(s1_feh * unc["[Fe/H]"]["uncertainty"]),
             abs(s1_teff * unc["Teff"]["uncertainty"]),
             d1.max() - d1.min()]
    cols = ["#C0392B", "#E67E22", "#95A5A6", "#95A5A6", "#2C3E50"]
    bars = ax.bar(names, sizes, color=cols)
    for rect, v in zip(bars, sizes, strict=True):
        ax.text(rect.get_x() + rect.get_width() / 2, v * 1.08, f"{v:.4f}",
                ha="center", fontsize=7.5)
    ax.set_yscale("log")
    ax.set_ylabel("size in $h_1$")
    ax.set_title("What could mimic the signal,\nand what could not", fontsize=10)
    ax.tick_params(axis="x", labelsize=7)
    ax.grid(alpha=0.25, axis="y")

    # 3. The facular chord test.
    ax = axes[2]
    h1, h2 = lib_set1.h(SAMPLE["FeH"], SAMPLE["Teff"], SAMPLE["logg"])
    i2 = h1 - h2
    need = (refld["dh1"] - refld["dh2"]) / i2
    contrasts = np.linspace(0.01, 0.20, 200)
    ax.plot(contrasts * 100, need / contrasts, color="#2C3E50", lw=2)
    for entry in load_facular_contrast():
        c = entry["contrast_peak"]
        ax.plot(c * 100, need / c, "o", ms=7)
        ax.annotate(f"{entry['source'].split(' et')[0]}\n{entry['feature']}",
                    (c * 100, need / c), textcoords="offset points",
                    xytext=(8, 6), fontsize=6.5)
    ax.axhspan(0.3, 1.5, color="#C0392B", alpha=0.07)
    ax.text(11, 0.75, "not credible for these\nquiet, bright transit hosts",
            fontsize=7, color="#C0392B")
    ax.set_xlabel("peak facular contrast at $\\mu=1/3$ (%)")
    ax.set_ylabel("filling factor needed on the chord")
    ax.set_ylim(0, 1.5)
    ax.set_title("Bright features on the chord cannot\ndo it at any measured contrast",
                 fontsize=10)
    ax.grid(alpha=0.25)

    # 4. Is h2 as good as h1?
    ax = axes[3]
    kep = next(r for r in kepler if r["is_refld"])
    groups = ["$h_1$", "$h_2$"]
    offsets = [abs(kep["dh1"]), abs(kep["dh2"])]
    reanalysis = [0.000, 0.010]          # Maxted 2023 section 4.2
    predicted = [0.007, 0.005]           # Norris 2017 MURaM 100 G, via Maxted 4.3.4
    x = np.arange(2)
    w = 0.26
    ax.bar(x - w, offsets, w, label="measured offset", color="#C0392B")
    ax.bar(x, predicted, w, label="MURaM 100 G prediction", color="#2C6E49")
    ax.bar(x + w, reanalysis, w, label="shift from reanalysis alone", color="#E67E22")
    for xi, vals in zip(x, zip(offsets, predicted, reanalysis, strict=True), strict=True):
        for dx, v in zip((-w, 0, w), vals, strict=True):
            ax.text(xi + dx, v + 0.0004, f"{v:.3f}", ha="center", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=12)
    ax.set_ylabel("absolute size")
    ax.set_ylim(0, 0.016)
    ax.set_title("$h_1$ carries the result.\n$h_2$ is the size of its own systematic",
                 fontsize=10)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.25, axis="y")

    fig.tight_layout()
    path = FIG / "systematics_budget.png"
    fig.savefig(path, dpi=140)
    print(f"Figure written to {path}")


if __name__ == "__main__":
    main()
