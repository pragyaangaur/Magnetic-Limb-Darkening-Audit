"""Two checks on the reference side of the magnetic limb-darkening result.

Analysis a01 shows the magnetic offset is a residual against one library, and
a04 shows the Sun is sharp enough to rank the libraries and then stops. This
script takes the next step and finds that one half of it cannot be taken yet
and the other half produces a number worth asking about.

PART ONE. The solar residual against MPS-ATLAS, which is the reference the 2024
result is measured against. The answer is that this comparison cannot be made
from the data in hand, and the reason is quantified rather than asserted: the
measured solar profiles are at 579.88 nm, the library is tabulated in broad
passbands, and moving between passbands moves model h1 by six times the signal.
The passband conversion is therefore not a detail. It is the whole measurement.

PART TWO. A check that has no wavelength problem, because both sides come from
the same library in the same passband. Kostogryz et al. (2024) restrict their
MURaM simulations to a solar atmosphere, on the stated grounds that the sample
has near-solar parameters. The sample mean is 583 K hotter than the Sun and
0.23 dex more metal rich. This measures what that gap is worth in the currency
of the result.

Run it with:

    python analysis/a07_solar_residual.py
"""

import csv
import pathlib

import numpy as np
from ldlib import (
    PASSBANDS,
    SAMPLE,
    LimbDarkeningLibrary,
    clv_intensity,
    load_solar_clv,
    load_table3,
)

MU1, MU2 = 2.0 / 3.0, 1.0 / 3.0

# IAU 2015 nominal solar values. [Fe/H] = 0 by definition of the scale.
SUN = {"Teff": 5772.0, "logg": 4.438, "FeH": 0.0}


def solar_h(coefficients):
    """h1 and h2 straight from a measured polynomial, with no law fitted."""
    i1 = float(clv_intensity(MU1, coefficients))
    i2 = float(clv_intensity(MU2, coefficients))
    return i1, i1 - i2


def main():
    rows = load_table3()
    kep = next(r for r in rows if r["passband"] == "Kepler" and r["is_refld"])
    tes = next(r for r in rows if r["passband"] == "TESS" and r["is_refld"])

    # ---------------------------------------------------------------- part one
    profiles = load_solar_clv()
    measured = {name: solar_h(p["coefficients"]) for name, p in profiles.items()}
    h1_meas = float(np.mean([v[0] for v in measured.values()]))

    print("=" * 72)
    print("PART ONE. The solar residual cannot be computed at one wavelength.")
    print("=" * 72)
    print()
    print("Measured solar profile at 579.88 nm, no law fitted:")
    print(f"  {'':<10}{'h1':>10}{'h2':>10}")
    for name in ("PS77", "NL94"):
        print(f"  {name:<10}{measured[name][0]:>10.5f}{measured[name][1]:>10.5f}")
    print(f"  {'spread':<10}"
          f"{abs(measured['NL94'][0] - measured['PS77'][0]):>10.5f}"
          f"{abs(measured['NL94'][1] - measured['PS77'][1]):>10.5f}")
    print()

    print("MPS-ATLAS Set 1 at solar parameters, by passband:")
    print(f"  {'passband':<10}{'h1':>10}{'h2':>10}{'h1 minus measured':>22}")
    model = {}
    for which in ("set1", "set2"):
        for band in PASSBANDS:
            lib = LimbDarkeningLibrary(which, band, method="clv")
            model[(which, band)] = lib.h(SUN["FeH"], SUN["Teff"], SUN["logg"])
    for band in PASSBANDS:
        h1, h2 = model[("set1", band)]
        print(f"  {band:<10}{h1:>10.5f}{h2:>10.5f}{h1_meas - h1:>+22.5f}")

    band_h1 = [model[("set1", b)][0] for b in PASSBANDS]
    band_spread1 = max(band_h1) - min(band_h1)
    print()
    print(f"  Model h1 moves by {band_spread1:.4f} across the four passbands, which is")
    print(f"  {band_spread1 / abs(kep['dh1']):.1f} times the {abs(kep['dh1']):.3f} offset under discussion. Limb darkening")
    print("  is weaker in the red, so h1 rises with effective wavelength, and TESS")
    print("  sits highest as it should.")
    print()
    print("  The measured value is at 579.88 nm and every model value is a broad")
    print("  band. The apparent residuals in the last column are therefore of")
    print("  unknown sign and size until the measurement is carried into the band.")
    print("  NO SOLAR RESIDUAL IS CLAIMED HERE. What is established is that the")
    print("  passband conversion is worth several times the signal, so it is not")
    print("  a refinement. It is the measurement. That is the object Kostogryz")
    print("  et al. (2022) section 5.2 calls NL_Kepler, built for their own")
    print("  library and applied to no other.")
    print()

    d1 = np.array([r["dh1"] for r in rows if r["passband"] == "Kepler"])
    print(f"  For scale: nine libraries span {d1.max() - d1.min():.4f} in h1, the two solar")
    print(f"  measurements agree to {abs(measured['NL94'][0] - measured['PS77'][0]):.4f}, "
          f"and the passband choice is worth {band_spread1:.4f}.")
    print("  The Sun is sharp enough. The bridge to the passband is what is missing.")
    print()

    # ---------------------------------------------------------------- part two
    print("=" * 72)
    print("PART TWO. What the near-solar assumption is worth.")
    print("=" * 72)
    print()
    print("Kostogryz et al. (2024) run MURaM for a solar atmosphere only, on the")
    print("stated grounds that the sample has near-solar parameters. Both sides")
    print("below come from the same library in the same passband, so there is no")
    print("wavelength mismatch and no interpolation between different sources.")
    print()
    print(f"  Sun    Teff = {SUN['Teff']:.0f} K, logg = {SUN['logg']:.2f}, [Fe/H] = {SUN['FeH']:+.2f}")
    print(f"  Sample Teff = {SAMPLE['Teff']:.0f} K, logg = {SAMPLE['logg']:.2f}, [Fe/H] = {SAMPLE['FeH']:+.2f}")
    print(f"  Gap    {SAMPLE['Teff'] - SUN['Teff']:+.0f} K, "
          f"{SAMPLE['logg'] - SUN['logg']:+.2f} dex, {SAMPLE['FeH'] - SUN['FeH']:+.2f} dex")
    print()
    print(f"  {'passband':<10}{'dh1':>10}{'dh2':>10}{'vs offset h1':>16}{'vs offset h2':>16}")
    print("  " + "-" * 62)
    results = {}
    for band, obs in (("Kepler", kep), ("TESS", tes)):
        lib = LimbDarkeningLibrary("set1", band, method="clv")
        hs = lib.h(SAMPLE["FeH"], SAMPLE["Teff"], SAMPLE["logg"])
        hq = lib.h(SUN["FeH"], SUN["Teff"], SUN["logg"])
        d_h1, d_h2 = hs[0] - hq[0], hs[1] - hq[1]
        results[band] = (d_h1, d_h2)
        print(f"  {band:<10}{d_h1:>+10.5f}{d_h2:>+10.5f}"
              f"{abs(d_h1 / obs['dh1']):>15.1f}x{abs(d_h2 / obs['dh2']):>15.1f}x")
    print()

    dh1, dh2 = results["Kepler"]
    print("READING")
    print()
    print("  Going from solar parameters to the sample mean changes non-magnetic")
    print(f"  limb darkening by {dh1:+.4f} in h1 and {dh2:+.4f} in h2 in the Kepler band.")
    print(f"  That is {abs(dh1 / kep['dh1']):.1f} times the h1 offset and "
          f"{abs(dh2 / kep['dh2']):.1f} times the h2 offset being explained,")
    print("  and the Kepler and TESS values agree, so it is not a passband artefact.")
    print()
    print("  Note the direction. The parameter gap pushes h1 up and h2 down, which")
    print("  is the same direction as the signature attributed to magnetic fields.")
    print()
    print("  What this does NOT say. It does not explain the offset away. Maxted")
    print("  computes the model prediction at each star's own parameters, so the")
    print("  measured residual is already differential in this respect, and a02")
    print("  separately rules out a parameter scale error.")
    print()
    print("  What it does say. The atmospheric structure changes a lot across the")
    print("  gap between the Sun and this sample, by several times the effect")
    print("  being modelled. Kostogryz et al. compute the magnetic CORRECTION from")
    print("  a solar box and apply it to stars well outside solar parameters. The")
    print("  paper does not bound the error in that step, and their own 2026")
    print("  follow-up reports that the magnetic effect grows towards hotter and")
    print("  more metal-rich stars, which is the direction this sample lies in.")
    print()
    print("  That is a question with a short answer, and only the authors have it:")
    print("  how much does the MURaM 100 G prediction move between a solar box and")
    print("  one at 6355 K and [Fe/H] = 0.23?")


    _write_tables(model, measured, results)


def _write_tables(model, measured, results):
    """Regenerable tables for results/."""
    out = pathlib.Path(__file__).resolve().parents[1] / "results" / "tables"
    out.mkdir(parents=True, exist_ok=True)

    with (out / "a07_solar_passband_spread.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# MPS-ATLAS at solar parameters (Teff 5772, logg 4.438, [Fe/H] 0)."])
        w.writerow(["# Measured solar values are at 579.88 nm and are NOT in-band."])
        w.writerow(["# The last two columns are therefore not residuals. See a07 part one."])
        w.writerow(["set", "passband", "h1", "h2",
                    "measured_minus_model_h1", "measured_minus_model_h2"])
        h1m = sum(v[0] for v in measured.values()) / len(measured)
        h2m = sum(v[1] for v in measured.values()) / len(measured)
        for (which, band), (h1, h2) in model.items():
            w.writerow([which, band, f"{h1:.5f}", f"{h2:.5f}",
                        f"{h1m - h1:+.5f}", f"{h2m - h2:+.5f}"])

    with (out / "a07_near_solar_assumption.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# Change in non-magnetic MPS-ATLAS limb darkening between solar"])
        w.writerow(["# parameters and the Maxted (2023) sample mean, same library and"])
        w.writerow(["# same passband on both sides, so no wavelength mismatch."])
        w.writerow(["passband", "dh1_sample_minus_solar", "dh2_sample_minus_solar"])
        for band, (d1, d2) in results.items():
            w.writerow([band, f"{d1:+.5f}", f"{d2:+.5f}"])


if __name__ == "__main__":
    main()
