"""Cross-checks against numbers published in the source papers.

Every check here compares something this code computes against a number printed
in a paper. If they all pass, the table readers, the interpolation, the h1 and
h2 definitions and the unit convention are all correct.

Run it with:

    python analysis/verify.py

It exits non-zero if any check fails.
"""

import sys

import numpy as np
from ldlib import (
    SAMPLE,
    LimbDarkeningLibrary,
    clv_intensity,
    load_solar_clv,
    load_table3,
)

CHECKS = []


def check(name, got, want, tol, source):
    ok = abs(got - want) <= tol
    CHECKS.append((name, got, want, tol, source, ok))
    return ok


def main():
    T, G = SAMPLE["Teff"], SAMPLE["logg"]

    # 1. Maxted (2023) section 4.3.1 quotes the metallicity correction he
    #    computed from Set 1 of Kostogryz et al. (2022): at Teff = 6355 K and
    #    log g = 4.39, using [Fe/H] = 0 instead of 0.23 makes h1 too high by
    #    0.0034 and h2 too low by 0.0041 in the Kepler band, and 0.0027 and
    #    0.0035 in TESS. He states he uses the tabulated intensities directly.
    for band, (want1, want2) in {"Kepler": (0.0034, -0.0041),
                                 "TESS": (0.0027, -0.0035)}.items():
        lib = LimbDarkeningLibrary("set1", band, method="clv")
        h_zero = lib.h(0.00, T, G)
        h_true = lib.h(0.23, T, G)
        check(f"Maxted metallicity correction, {band} h1",
              h_zero[0] - h_true[0], want1, 0.0002, "Maxted 2023 section 4.3.1")
        check(f"Maxted metallicity correction, {band} h2",
              h_zero[1] - h_true[1], want2, 0.0002, "Maxted 2023 section 4.3.1")

    # 2. Maxted (2023) Table 3 gives Set 1 and Set 2 offsets in the Kepler band
    #    as -0.006 and -0.009. The difference, 0.003, is a property of the two
    #    libraries alone and can be recomputed from the CDS tables.
    rows = load_table3()
    for band, tol in (("Kepler", 0.0015), ("TESS", 0.0015)):
        band_rows = [r for r in rows if r["passband"] == band]
        s1 = next(r for r in band_rows if "Set 1" in r["library"])
        s2 = next(r for r in band_rows if "Set 2" in r["library"])
        want = s2["dh1"] - s1["dh1"]
        lib1 = LimbDarkeningLibrary("set1", band, method="clv")
        lib2 = LimbDarkeningLibrary("set2", band, method="clv")
        got = lib1.h(SAMPLE["FeH"], T, G)[0] - lib2.h(SAMPLE["FeH"], T, G)[0]
        check(f"Set1 minus Set2 in h1, {band}", got, want, tol,
              "Maxted 2023 Table 3, two MPS-ATLAS rows")

    # 3. The unit convention. Kostogryz et al. (2024) quote the Kepler offsets
    #    as 0.6 per cent and -1.2 per cent. Maxted's abstract states the same
    #    offset as "dh1 ~ 0.006". If the percentages were relative they would
    #    be 0.005 and -0.0023 in absolute units, and the h2 value would not
    #    match Table 3 at all.
    kepler = [r for r in rows if r["passband"] == "Kepler"]
    refld = next(r for r in kepler if r["is_refld"])
    check("REFLD Kepler h1 offset equals the 0.6 per cent figure",
          refld["dh1"], 0.006, 1e-9, "Kostogryz 2024 Methods, Maxted 2023 Table 3")
    check("REFLD Kepler h2 offset equals the -1.2 per cent figure",
          refld["dh2"], -0.012, 1e-9, "Kostogryz 2024 Methods, Maxted 2023 Table 3")

    # 4. The solar profiles. Hestroffer and Magnan (1998) state that the
    #    Pierce and Slaughter and Neckel and Labs polynomials agree well for
    #    r/R below 0.9, which is mu above 0.436. Check that the two agree to
    #    better than 0.3 per cent there and disagree by more than 5 per cent at
    #    the limb, which is what their Figure 1 shows.
    profiles = load_solar_clv()
    ps = profiles["PS77"]["coefficients"]
    nl = profiles["NL94"]["coefficients"]
    mu = np.linspace(0.436, 1.0, 200)
    worst = np.max(np.abs(clv_intensity(mu, nl) - clv_intensity(mu, ps))
                   / clv_intensity(mu, ps))
    check("solar profiles agree inside r/R = 0.9", worst, 0.0, 0.003,
          "Hestroffer and Magnan 1998, their Figure 1")
    limb = abs(clv_intensity(0.0, nl) - clv_intensity(0.0, ps)) / clv_intensity(0.0, ps)
    check("solar profiles disagree at the limb", limb, 0.069, 0.01,
          "Hestroffer and Magnan 1998 Table 1, computed from the coefficients")

    # 5. Both solar polynomials are normalised so the coefficients sum to one.
    for name, coeffs in (("PS77", ps), ("NL94", nl)):
        check(f"{name} coefficients sum to one", float(sum(coeffs)), 1.0, 2e-5,
              "Hestroffer and Magnan 1998, their equation for P5")

    # 6. Sanity: h1 must lie between the intensity at mu = 1 and at the limb,
    #    and h2 must be positive for any limb-darkened star.
    lib = LimbDarkeningLibrary("set1", "Kepler", method="clv")
    h1, h2 = lib.h(SAMPLE["FeH"], T, G)
    check("h1 is between 0 and 1", float(0 < h1 < 1), 1.0, 0.0, "definition")
    check("h2 is positive", float(h2 > 0), 1.0, 0.0, "definition")

    width = max(len(c[0]) for c in CHECKS)
    print(f"{'check':<{width}}{'computed':>12}{'published':>12}{'diff':>11}  ")
    print("-" * (width + 38))
    for name, got, want, tol, source, ok in CHECKS:
        print(f"{name:<{width}}{got:>12.5f}{want:>12.5f}{got - want:>+11.5f}"
              f"  {'pass' if ok else 'FAIL'}")
    print()
    for name, *_, source, ok in CHECKS:
        if not ok:
            print(f"  failed: {name}   source: {source}")

    failed = sum(1 for c in CHECKS if not c[-1])
    print(f"{len(CHECKS) - failed} of {len(CHECKS)} checks pass.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
