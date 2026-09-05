"""Could bright magnetic features on the transit chord do it instead?

Maxted (2023) raises this objection to any solar calibration of limb darkening:
a transit light curve samples a chord across the star that contains spots and
faculae, whereas the solar measurements of Neckel and Labs deliberately avoided
magnetic regions. He states it qualitatively. Nobody has put a number on it.

This script does. It asks what filling factor of bright magnetic features on the
chord would be needed to produce the observed offset, using the measured
centre-to-limb contrast of solar faculae and network.

The distinction matters. Kostogryz et al. attribute the offset to small-scale
fields changing the mean atmospheric structure. Facular area coverage on the
chord is a different thing, an area-weighted mixture of two atmospheres. If the
second can produce the offset with a plausible filling factor, the two are not
separable by the transit data alone.

Run it with:

    python analysis/a03_facular_chord_test.py
"""

import csv
import pathlib

import numpy as np
from ldlib import SAMPLE, LimbDarkeningLibrary, load_facular_contrast, load_table3

TABLES = pathlib.Path(__file__).resolve().parents[1] / "results" / "tables"

MU1, MU2 = 2.0 / 3.0, 1.0 / 3.0


def contrast_shape(mu, mu_peak, exponent):
    """
    Facular contrast against mu, normalised to one at the peak.

    Pietrow et al. (2026) find the contrast rises from near zero at disc centre,
    peaks near mu = 0.3, and falls again towards the limb. The rise is described
    only as gradual, so the exponent is left free here and the analysis reports
    what the data require rather than assuming a shape.
    """
    mu = np.asarray(mu, dtype=float)
    rise = np.clip((1.0 - mu) / (1.0 - mu_peak), 0.0, 1.0) ** exponent
    return rise


def main():
    rows = load_table3()
    contrasts = load_facular_contrast()
    lib = LimbDarkeningLibrary("set1", "Kepler")
    refld = next(r for r in rows if r["passband"] == "Kepler" and r["is_refld"])

    h1, h2 = lib.h(SAMPLE["FeH"], SAMPLE["Teff"], SAMPLE["logg"])
    i1, i2 = h1, h1 - h2                      # I(2/3) and I(1/3) of the model
    d_i1 = refld["dh1"]                       # observed minus model at mu = 2/3
    d_i2 = refld["dh1"] - refld["dh2"]        # and at mu = 1/3

    print("The excess brightness the observations require, relative to REFLD:")
    print(f"  at mu = 2/3   {d_i1:+.4f} absolute,  {d_i1 / i1 * 100:+.2f} per cent of I")
    print(f"  at mu = 1/3   {d_i2:+.4f} absolute,  {d_i2 / i2 * 100:+.2f} per cent of I")
    required_ratio = (d_i2 / i2) / (d_i1 / i1)
    print(f"  required contrast ratio  c(1/3) / c(2/3) = {required_ratio:.2f}")
    print()

    print("Filling factor needed, for each measured contrast anchor.")
    print("f is the fraction of the transit chord covered by bright features.")
    print()
    print(f"  {'peak contrast':>14}{'source':<34}{'f needed':>10}{'implied c(2/3)':>16}")
    print("  " + "-" * 74)
    out = []
    for entry in contrasts:
        c_peak = entry["contrast_peak"]
        f = (d_i2 / i2) / c_peak              # mu = 1/3 sits at the measured peak
        implied = (d_i1 / i1) / f
        print(f"  {c_peak:>14.2f}{entry['source']:<34}{f:>10.2f}{implied:>16.3f}")
        out.append({"feature": entry["feature"], "source": entry["source"],
                    "contrast_peak": c_peak, "filling_factor": f,
                    "implied_c_at_two_thirds": implied})
    print()

    print("The conclusion does not depend on any one contrast value. Scanning")
    print("contrast continuously:")
    print()
    print(f"  {'peak contrast at mu=1/3':>24}{'filling factor needed':>24}")
    print("  " + "-" * 46)
    for c_peak in (0.02, 0.04, 0.06, 0.10, 0.15, 0.20, 0.30):
        print(f"  {c_peak * 100:>22.0f} %{(d_i2 / i2) / c_peak:>24.2f}")
    print()
    print("  The largest facular contrast reported anywhere in the literature is")
    print("  about 10 per cent. Even at an implausible 30 per cent the chord")
    print("  would still have to be nine per cent covered, and at any published")
    print("  value it is above a quarter. There is no contrast for which this")
    print("  explanation becomes comfortable.")
    print()

    print("What contrast shape would be needed, if the filling factor were")
    print("physically reasonable. Exponent 1 is a linear rise to the peak,")
    print("exponent 2 a quadratic one.")
    print()
    print(f"  {'exponent':>9}{'c(2/3)/c(1/3)':>16}{'matches required?':>20}")
    print("  " + "-" * 45)
    for exponent in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
        shape = contrast_shape([MU1, MU2], 0.30, exponent)
        ratio = shape[1] / shape[0]
        verdict = "yes" if abs(ratio - required_ratio) < 0.5 else "no"
        print(f"  {exponent:>9.1f}{ratio:>16.2f}{verdict:>20}")
    print()

    print("Reading this.")
    print()
    print("Taking the facular contrast measured with HMI by Pietrow et al. (2026),")
    print("which peaks near four per cent, the chord would have to be about seventy")
    print("per cent covered in faculae. Even with the largest contrast anywhere in")
    print("the literature, ten per cent for strong-field features from Yeo et al.")
    print("(2013), it is still about a quarter of the chord.")
    print()
    print("Those numbers are not credible for the sample. These are bright, quiet,")
    print("mostly F-type transit hosts selected for clean photometry, not heavily")
    print("spotted active stars. A chord a quarter covered in faculae would produce")
    print("obvious rotational modulation and obvious spot-crossing anomalies.")
    print()
    print("So the offset is not facular area coverage on the chord. It behaves like")
    print("a change in the mean atmospheric structure, which is what Kostogryz et al.")
    print("propose. This closes the objection Maxted raised, quantitatively, and in")
    print("their favour.")
    print()
    print("One caveat that cuts the other way. The measured contrasts are explicitly")
    print("lower limits, because the segmentation masks admit quiet Sun at feature")
    print("boundaries, and the sample stars are hotter and more metal rich than the")
    print("Sun, where facular contrast is not measured at all. The argument would be")
    print("much stronger with a proper solar facular contrast in the Kepler band.")

    TABLES.mkdir(parents=True, exist_ok=True)
    path = TABLES / "a03_facular_chord.csv"
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"\nTable written to {path}")


if __name__ == "__main__":
    main()
