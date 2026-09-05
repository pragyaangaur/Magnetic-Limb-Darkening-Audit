"""Does the h2 offset carry the weight placed on it?

Kostogryz et al. (2024) present the joint fit to both steepness parameters as
the strength of their result: the magnetised models "simultaneously explain the
offsets in both limb darkening coefficients". Half of that joint constraint is
h2.

Maxted (2023), the paper the observed offsets come from, contains two warnings
about h2 that the 2024 paper does not repeat. This script puts numbers on both.

Run it with:

    python analysis/a06_h2_reliability.py
"""

import csv
import pathlib

from ldlib import DATA, load_table3

TABLES = pathlib.Path(__file__).resolve().parents[1] / "results" / "tables"


def load_h2_evidence():
    with (DATA / "maxted2023_h2_systematic.csv").open() as handle:
        rows = csv.DictReader(line for line in handle if not line.startswith("#"))
        return {r["quantity"]: {
            "value": float(r["value"]),
            "uncertainty": float(r["uncertainty"]) if r["uncertainty"] else None,
            "source": r["source"], "note": r["note"]} for r in rows}


def main():
    evidence = load_h2_evidence()
    rows = load_table3()

    print("Warning one: the same stars analysed twice.")
    print()
    print("Maxted (2023) section 4.2 compares sixteen stars that also appear in")
    print("Maxted (2018). The two analyses differ only in data processing and")
    print("fitting choices, not in the data. The shifts he finds are:")
    print()
    r1, r2 = evidence["repeat_analysis_dh1"], evidence["repeat_analysis_dh2"]
    print(f"  h1 shifted by {r1['value']:+.3f} +/- {r1['uncertainty']:.3f}")
    print(f"  h2 shifted by {r2['value']:+.3f} +/- {r2['uncertainty']:.3f}")
    print()
    print("  His conclusion: h1 is robust, h2 may be affected by systematic")
    print("  errors of about 0.01 depending on the details of the analysis.")
    print()

    print("Now set that against the offsets being interpreted.")
    print()
    print(f"  {'passband':<10}{'offset h1':>12}{'offset h2':>12}"
          f"{'h2 / systematic':>18}")
    print("  " + "-" * 52)
    out = []
    for band in ("Kepler", "TESS"):
        refld = next(r for r in rows if r["passband"] == band and r["is_refld"])
        ratio = abs(refld["dh2"]) / r2["value"]
        print(f"  {band:<10}{refld['dh1']:>+12.3f}{refld['dh2']:>+12.3f}{ratio:>17.1f}x")
        out.append({"passband": band, "dh1": refld["dh1"], "dh2": refld["dh2"],
                    "h2_analysis_systematic": r2["value"],
                    "ratio_h2_to_systematic": ratio})
    print()
    print("  The h2 offset is the same size as the shift that data-processing")
    print("  choices alone produced on the same stars. It is a mean over 24")
    print("  objects, so random scatter averages down, but an analysis choice")
    print("  applies to every star at once and does not.")
    print()

    print("Warning two: the magnetic prediction fits h1 and not h2.")
    print()
    m1, m2 = evidence["muram_100G_dh1"], evidence["muram_100G_dh2"]
    kepler = next(r for r in rows if r["passband"] == "Kepler" and r["is_refld"])
    print("  Maxted's own section 4.3.4 reads the MURaM 100 G calculation of")
    print("  Norris et al. (2017) at 611 nm, near the Kepler mean wavelength of")
    print("  630 nm, and quotes the predicted magnetic effect:")
    print()
    print(f"  {'':<12}{'predicted':>12}{'observed':>12}{'ratio':>10}")
    print("  " + "-" * 46)
    print(f"  {'h1':<12}{m1['value']:>+12.3f}{kepler['dh1']:>+12.3f}"
          f"{kepler['dh1'] / m1['value']:>10.2f}")
    print(f"  {'h2':<12}{m2['value']:>+12.3f}{kepler['dh2']:>+12.3f}"
          f"{kepler['dh2'] / m2['value']:>10.2f}")
    print()
    print("  h1 agrees to fifteen per cent. h2 is out by a factor of 2.4, in")
    print("  the sense that the observed offset is larger than the magnetic")
    print("  prediction. Maxted says as much: the h2 comparison is not")
    print("  straightforward, but the h1 conclusion is robust.")
    print()

    print("What follows.")
    print()
    print("The h1 result is solid on every count. It is stable against")
    print("reanalysis, it matches an independent MURaM prediction to fifteen")
    print("per cent, and no non-magnetic explanation I could construct survives.")
    print()
    print("The h2 result is not in the same condition. Its offset is the size")
    print("of a known analysis systematic, and the magnetic prediction that")
    print("matches h1 undershoots it by a factor of two. Presenting the joint")
    print("fit to both parameters as the strength of the case rests on the")
    print("weaker half.")
    print()
    print("This is not a claim that the magnetic interpretation is wrong. It is")
    print("a claim that its evidence is one parameter rather than two, and that")
    print("the warning is in the paper the numbers were taken from.")

    TABLES.mkdir(parents=True, exist_ok=True)
    path = TABLES / "a06_h2_reliability.csv"
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"\nTable written to {path}")


if __name__ == "__main__":
    main()
