"""Every number and the figure in the RNAAS note, from the data in data/.

The note asks whether two non-magnetic explanations can produce the limb-darkening
offset that Maxted (2023) measured and Kostogryz et al. (2024) attribute to surface
magnetic fields. The explanations are a systematic error in the adopted stellar
parameter scale, and bright magnetic features covering the transit chord.

Maxted (2023) section 4.2 finds that h2' may carry an analysis systematic of about
0.01, while h1' is robust. So each test is run twice: once with the published h2'
offset, and once with h1' alone or with the h2' offset moved by that 0.01 towards
zero. A conclusion is only stated in the note if it survives the second version.

Run it from the repository root with:

    python rnaas/build.py
"""

import json
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))

from ldlib import SAMPLE, LimbDarkeningLibrary, load_table3

HERE = pathlib.Path(__file__).resolve().parent

# Maxted (2023) section 4.2: h2' may be affected by systematic errors of about 0.01.
H2_SYSTEMATIC = 0.010

# Tayar et al. (2022): interferometric angular diameters and bolometric fluxes set
# a systematic floor of about 2.4 per cent on the effective temperature scale.
TEFF_FLOOR_FRACTION = 0.024

STEPS = {"FeH": 0.10, "Teff": 100.0, "logg": 0.10}

# Peak facular contrast near mu = 0.3 at HMI 6173 A. Pietrow et al. (2026) measure
# about 4 per cent for faculae and state it is a lower limit. They quote peak values
# approaching 10 per cent for strong-field features from Yeo et al. (2013).
CONTRAST_FACULAE = 0.04
CONTRAST_STRONG = 0.10

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dc"


def parameter_shifts(band):
    """The parameter scale error that would reproduce the offset, three ways."""
    lib = LimbDarkeningLibrary("set1", band)
    row = next(r for r in load_table3() if r["passband"] == band and r["is_refld"])
    d1, e1, d2, e2 = row["dh1"], row["dh1_err"], row["dh2"], row["dh2_err"]
    out = {}
    for name in ("FeH", "Teff", "logg"):
        s1, s2 = lib.derivative(SAMPLE["FeH"], SAMPLE["Teff"], SAMPLE["logg"],
                                name, STEPS[name])
        denom = s1**2 / e1**2 + s2**2 / e2**2
        joint = (d1 * s1 / e1**2 + d2 * s2 / e2**2) / denom
        d2_shifted = d2 + H2_SYSTEMATIC * (1 if d2 < 0 else -1)
        joint_shifted = (d1 * s1 / e1**2 + d2_shifted * s2 / e2**2) / denom
        out[name] = {
            "h1_only": d1 / s1, "h1_only_err": e1 / abs(s1),
            "joint": joint, "joint_err": denom**-0.5,
            "joint_h2_shifted": joint_shifted,
        }
    return out


def facular_requirement(band):
    """Relative brightening the offset needs at mu = 2/3 and 1/3, and the chord
    filling factor that bright features at a given peak contrast would need."""
    lib = LimbDarkeningLibrary("set1", band)
    row = next(r for r in load_table3() if r["passband"] == band and r["is_refld"])
    h1, h2 = lib.h(SAMPLE["FeH"], SAMPLE["Teff"], SAMPLE["logg"])
    out = {}
    for label, d2 in (("published", row["dh2"]),
                      ("h2_shifted", row["dh2"] + H2_SYSTEMATIC * (1 if row["dh2"] < 0 else -1))):
        excess_23 = row["dh1"] / h1                 # I(2/3) = h1
        excess_13 = (row["dh1"] - d2) / (h1 - h2)   # I(1/3) = h1 - h2
        out[label] = {
            "dh2": d2,
            "excess_two_thirds": excess_23,
            "excess_one_third": excess_13,
            # mu = 1/3 sits at the measured contrast peak, so f = excess / peak.
            "f_faculae": excess_13 / CONTRAST_FACULAE,
            "f_strong": excess_13 / CONTRAST_STRONG,
            "f_at_30_percent": excess_13 / 0.30,
            # h1' alone: contrast at mu = 2/3 cannot exceed the peak, so this is
            # a lower limit on f that does not use h2' at all.
            "f_min_h1_faculae": excess_23 / CONTRAST_FACULAE,
            "f_min_h1_strong": excess_23 / CONTRAST_STRONG,
        }
    return out


def figure(numbers, path):
    # AAS asks for common fonts and embedded TrueType rather than Type 3.
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    })
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.5),
                             gridspec_kw={"width_ratios": [1, 1, 1.5]})
    colours = {"Kepler": BLUE, "TESS": ORANGE}

    for ax, name, unit, label in ((axes[0], "FeH", "dex", "[Fe/H] shift needed (dex)"),
                                  (axes[1], "Teff", "K", "$T_{\\rm eff}$ shift needed (K)")):
        for i, band in enumerate(("Kepler", "TESS")):
            p = numbers["parameter_shifts"][band][name]
            ax.errorbar(i, p["h1_only"], yerr=p["h1_only_err"], fmt="o", ms=5,
                        color=colours[band], elinewidth=1.2, capsize=0)
            ax.plot(i + 0.22, p["joint"], marker="D", ms=4, ls="none",
                    mfc="white", mec=colours[band], mew=1.2)
        ax.axhline(0, color=MUTED, lw=0.6)
        ax.set_xticks([0.1, 1.1], ["Kepler", "TESS"])
        ax.set_xlim(-0.5, 1.6)
        ax.set_ylabel(label)
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    floor = numbers["teff_floor_K"]
    axes[1].axhspan(-floor, floor, color=GRID, lw=0, zorder=0)
    axes[1].text(1.55, -floor * 0.92, "2.4% scale\nfloor", ha="right", va="bottom",
                 fontsize=6.5, color=MUTED)
    axes[0].text(-0.45, -0.02, "filled: $h_1'$ alone\nopen: $h_1'$ and $h_2'$",
                 fontsize=6.5, color=MUTED, va="top")

    ax = axes[2]
    contrast = [c / 100 for c in range(2, 31)]
    for band in ("Kepler", "TESS"):
        for label, style in (("published", "-"), ("h2_shifted", "--")):
            e = numbers["facular"][band][label]["excess_one_third"]
            ax.plot([c * 100 for c in contrast], [e / c for c in contrast],
                    color=colours[band], lw=1.6, ls=style)
        # Name each band on its solid line, so identity does not rest on colour.
        e = numbers["facular"][band]["published"]["excess_one_third"]
        ax.text(29.5, e / 0.295 * 1.12, band, color=INK, fontsize=7,
                ha="right", va="bottom")
    for c in (CONTRAST_FACULAE, CONTRAST_STRONG):
        ax.axvline(c * 100, color=MUTED, lw=0.6, ls=":")
    ax.text(CONTRAST_FACULAE * 100 + 0.4, 1.02, "faculae", fontsize=6.5, color=MUTED)
    ax.text(CONTRAST_STRONG * 100 + 0.4, 1.02, "strong field", fontsize=6.5, color=MUTED)
    ax.set_yscale("log")
    ax.set_ylim(0.004, 1.2)
    ax.set_xlim(2, 30)
    ax.set_xlabel("peak facular contrast (%)")
    ax.set_ylabel("chord filling factor needed")
    ax.grid(axis="y", color=GRID, lw=0.5, which="both")
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.plot([], [], color=MUTED, lw=1.2, ls="-", label="published $h_2'$")
    ax.plot([], [], color=MUTED, lw=1.2, ls="--", label="$h_2'$ moved by 0.01")
    ax.legend(fontsize=6.5, frameon=False, loc="lower left")

    for ax, tag in zip(axes, "abc"):
        ax.text(-0.02, 1.04, f"({tag})", transform=ax.transAxes, fontsize=8,
                color=INK, ha="right")
    fig.tight_layout(w_pad=1.5)
    fig.savefig(path)
    fig.savefig(path.with_suffix(".png"), dpi=200)


def main():
    numbers = {
        "sample": SAMPLE,
        "h2_systematic": H2_SYSTEMATIC,
        "teff_floor_K": TEFF_FLOOR_FRACTION * SAMPLE["Teff"],
        "parameter_shifts": {b: parameter_shifts(b) for b in ("Kepler", "TESS")},
        "facular": {b: facular_requirement(b) for b in ("Kepler", "TESS")},
    }

    print(f"Teff systematic floor at the sample mean: {numbers['teff_floor_K']:.0f} K")
    for band in ("Kepler", "TESS"):
        print(f"\n{band}")
        for name, p in numbers["parameter_shifts"][band].items():
            print(f"  {name:<5} h1' alone {p['h1_only']:+9.3f} +/- {p['h1_only_err']:.3f}"
                  f"   joint {p['joint']:+9.3f} +/- {p['joint_err']:.3f}"
                  f"   joint, h2' moved {p['joint_h2_shifted']:+9.3f}")
        for label, f in numbers["facular"][band].items():
            print(f"  faculae, {label:<10} dh2' {f['dh2']:+.3f}"
                  f"  excess {f['excess_two_thirds'] * 100:.2f}% at 2/3,"
                  f" {f['excess_one_third'] * 100:.2f}% at 1/3"
                  f"  f = {f['f_faculae']:.2f} (4%), {f['f_strong']:.2f} (10%),"
                  f" {f['f_at_30_percent']:.3f} (30%)"
                  f"  h1' alone f >= {f['f_min_h1_faculae']:.2f} (4%),"
                  f" {f['f_min_h1_strong']:.3f} (10%)")

    (HERE / "numbers.json").write_text(json.dumps(numbers, indent=2) + "\n")
    figure(numbers, HERE / "figure.pdf")
    print("\nWrote rnaas/numbers.json, rnaas/figure.pdf and rnaas/figure.png")


if __name__ == "__main__":
    main()
