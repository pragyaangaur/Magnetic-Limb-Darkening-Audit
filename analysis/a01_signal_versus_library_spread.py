"""How big is the magnetic signal next to the spread between model libraries?

Kostogryz et al. (2024) take the offset between transit-measured limb darkening
and one non-magnetic model library, MPS-ATLAS Set 1, and attribute it to
small-scale surface magnetic fields. Maxted (2023) measured that offset against
nine different published libraries in the Kepler band and five in TESS. This
script puts the chosen offset next to the disagreement among the libraries.

Run it with:

    python analysis/a01_signal_versus_library_spread.py
"""

import csv
import pathlib

import numpy as np
from ldlib import load_table3

TABLES = pathlib.Path(__file__).resolve().parents[1] / "results" / "tables"


def main():
    rows = load_table3()
    out_rows = []

    for band in ("Kepler", "TESS"):
        band_rows = [r for r in rows if r["passband"] == band]
        refld = next(r for r in band_rows if r["is_refld"])
        d1 = np.array([r["dh1"] for r in band_rows])
        d2 = np.array([r["dh2"] for r in band_rows])

        print(f"=== {band} ===")
        print(f"  {'model library':<40}{'dh1':>9}{'dh2':>9}{'':>4}")
        print("  " + "-" * 62)
        for r in band_rows:
            mark = "  <- REFLD" if r["is_refld"] else ""
            print(f"  {r['library']:<40}{r['dh1']:>+9.3f}{r['dh2']:>+9.3f}{mark}")
        print("  " + "-" * 62)
        print(f"  {'range across libraries':<40}{d1.max() - d1.min():>9.3f}"
              f"{d2.max() - d2.min():>9.3f}")
        print(f"  {'standard deviation across libraries':<40}{d1.std(ddof=1):>9.3f}"
              f"{d2.std(ddof=1):>9.3f}")
        print()
        print("  The signal attributed to magnetic fields is the REFLD row:")
        print(f"      dh1 = {refld['dh1']:+.3f} +/- {refld['dh1_err']:.3f}"
              f"      dh2 = {refld['dh2']:+.3f} +/- {refld['dh2_err']:.3f}")
        print("  The full range across published non-magnetic libraries is")
        print(f"      {d1.max() - d1.min():.3f} in h1 and {d2.max() - d2.min():.3f} in h2,")
        print(f"  which is {(d1.max() - d1.min()) / abs(refld['dh1']):.1f} and "
              f"{(d2.max() - d2.min()) / abs(refld['dh2']):.1f} times the signal.")

        opposite = [r for r in band_rows if np.sign(r["dh1"]) != np.sign(refld["dh1"])]
        if opposite:
            print(f"  {len(opposite)} of {len(band_rows)} libraries give dh1 of the "
                  "opposite sign:")
            for r in opposite:
                print(f"      {r['library']} ({r['note'] or 'no note'})")
        print()

        out_rows.append({
            "passband": band,
            "refld_dh1": refld["dh1"], "refld_dh2": refld["dh2"],
            "range_dh1": d1.max() - d1.min(), "range_dh2": d2.max() - d2.min(),
            "std_dh1": d1.std(ddof=1), "std_dh2": d2.std(ddof=1),
            "n_libraries": len(band_rows),
            "ratio_dh1": (d1.max() - d1.min()) / abs(refld["dh1"]),
            "ratio_dh2": (d2.max() - d2.min()) / abs(refld["dh2"]),
        })

    print("What this does and does not say.")
    print()
    print("It does not say the measurement is wrong. Maxted's offsets are")
    print("measured consistently and the error on the mean is small.")
    print()
    print("It says that calling the residual against one particular library a")
    print("magnetic effect is a statement about that library. Another published")
    print("non-magnetic library gives a residual of the opposite sign in the")
    print("Kepler band. Nothing in the transit data chooses between them.")
    print()
    print("The TESS spread is much smaller, but only five libraries were compared")
    print("there and the two extreme Kepler libraries are not among them, so that")
    print("is not evidence that the TESS reference is better constrained.")

    TABLES.mkdir(parents=True, exist_ok=True)
    path = TABLES / "a01_library_spread.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out_rows[0]))
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"\nTable written to {path}")


if __name__ == "__main__":
    main()
