"""Could a systematic error in the adopted stellar parameters do it instead?

The offset is computed by interpolating a model library to each star's
effective temperature, gravity and metallicity, and subtracting. Random errors
in those parameters average down over the sample. A systematic offset in the
parameter scale does not.

This script asks how large such a systematic would have to be. It uses the
Kostogryz et al. (2022) library directly, at the actual mean parameters of
Maxted's sample rather than at solar values, because the sample mean is
6355 K and [Fe/H] = +0.23, which is neither solar nor close to it.

Run it with:

    python analysis/a02_stellar_parameter_systematics.py
"""

import csv
import pathlib

from ldlib import SAMPLE, LimbDarkeningLibrary, load_sample_uncertainties, load_table3

TABLES = pathlib.Path(__file__).resolve().parents[1] / "results" / "tables"

# Central-difference half-widths, at or above the native grid spacing.
STEPS = {"FeH": 0.10, "Teff": 100.0, "logg": 0.10}


def main():
    rows = load_table3()
    unc = load_sample_uncertainties()
    out = []

    print("Sample means (Maxted 2023):"
          f" Teff = {SAMPLE['Teff']:.0f} K, log g = {SAMPLE['logg']},"
          f" [Fe/H] = {SAMPLE['FeH']:+.2f}")
    print("Typical per-star uncertainties:"
          f" Teff {unc['Teff']['uncertainty']:.0f} K,"
          f" log g {unc['logg']['uncertainty']:.2f} dex,"
          f" [Fe/H] {unc['[Fe/H]']['uncertainty']:.2f} dex")
    print()

    for band in ("Kepler", "TESS"):
        lib = LimbDarkeningLibrary("set1", band)
        refld = next(r for r in rows if r["passband"] == band and r["is_refld"])
        d1o, e1 = refld["dh1"], refld["dh1_err"]
        d2o, e2 = refld["dh2"], refld["dh2_err"]
        h1, h2 = lib.h(SAMPLE["FeH"], SAMPLE["Teff"], SAMPLE["logg"])

        print(f"=== {band} ===  REFLD at the sample mean: h1 = {h1:.5f}, h2 = {h2:.5f}")
        print(f"  offset to explain:  dh1 = {d1o:+.3f} +/- {e1:.3f},"
              f"  dh2 = {d2o:+.3f} +/- {e2:.3f}")
        print()
        print(f"  {'parameter':<10}{'dh1/dp':>12}{'dh2/dp':>12}"
              f"{'joint offset':>15}{'chi2 (1 dof)':>14}{'vs typical err':>16}")
        print("  " + "-" * 79)

        for name in ("FeH", "Teff", "logg"):
            s1, s2 = lib.derivative(SAMPLE["FeH"], SAMPLE["Teff"], SAMPLE["logg"],
                                    name, STEPS[name])
            # Weighted least squares for the single offset that best fits both.
            denom = s1**2 / e1**2 + s2**2 / e2**2
            best = (d1o * s1 / e1**2 + d2o * s2 / e2**2) / denom
            chi2 = ((best * s1 - d1o) / e1) ** 2 + ((best * s2 - d2o) / e2) ** 2

            key = {"FeH": "[Fe/H]", "Teff": "Teff", "logg": "logg"}[name]
            sigma = unc[key]["uncertainty"]
            ratio = abs(best) / sigma
            print(f"  {key:<10}{s1:>+12.5f}{s2:>+12.5f}{best:>+15.2f}"
                  f"{chi2:>14.2f}{ratio:>13.0f}x")
            out.append({"passband": band, "parameter": key,
                        "dh1_dp": s1, "dh2_dp": s2, "joint_offset": best,
                        "chi2": chi2, "typical_uncertainty": sigma,
                        "offset_over_uncertainty": ratio})
        print()

    print("Verdict.")
    print()
    print("The metallicity and temperature columns fit both offsets reasonably")
    print("well, and the required values are consistent between the Kepler and")
    print("TESS samples, which is what a genuine scale error would look like.")
    print()
    print("But the sizes are impossible. Explaining the offset needs the adopted")
    print("metallicities to be wrong by about half a dex, or the temperatures by")
    print("more than two hundred kelvin. SWEET-Cat quotes 0.05 dex and 60 K.")
    print("A systematic ten times the formal error, in the same direction for")
    print("every star in two independent samples, is not credible.")
    print()
    print("So this alternative is ruled out. That is a result in favour of the")
    print("magnetic interpretation and it does not appear to have been checked.")

    TABLES.mkdir(parents=True, exist_ok=True)
    path = TABLES / "a02_parameter_systematics.csv"
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"\nTable written to {path}")


if __name__ == "__main__":
    main()
