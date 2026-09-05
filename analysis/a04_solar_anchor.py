"""Is the Sun accurate enough to choose between the model libraries?

Analysis a01 shows the offset attributed to magnetic fields is smaller than the
disagreement among published non-magnetic libraries. Nothing in the transit data
breaks that. Something outside it has to.

The Sun is the only candidate. Its effective temperature, gravity and
composition are known independently, and its centre-to-limb variation is
measured rather than modelled. So the library that gets the Sun right, at
exactly solar parameters, is the one whose stellar residual can be trusted.

The question this script answers is whether the existing solar measurements are
precise enough for that test to discriminate.

Run it with:

    python analysis/a04_solar_anchor.py
"""

import numpy as np
from ldlib import clv_intensity, load_solar_clv, load_table3

MU1, MU2 = 2.0 / 3.0, 1.0 / 3.0


def solar_h(coefficients):
    i1 = float(clv_intensity(MU1, coefficients))
    i2 = float(clv_intensity(MU2, coefficients))
    return i1, i1 - i2


def main():
    profiles = load_solar_clv()
    rows = load_table3()
    ps = profiles["PS77"]["coefficients"]
    nl = profiles["NL94"]["coefficients"]

    h1_ps, h2_ps = solar_h(ps)
    h1_nl, h2_nl = solar_h(nl)
    spread1 = abs(h1_nl - h1_ps)
    spread2 = abs(h2_nl - h2_ps)

    print("The two measured solar centre-to-limb profiles, at 579.88 nm.")
    print("Steepness parameters computed directly, with no law fitted.")
    print()
    print(f"  {'':<12}{'h1':>10}{'h2':>10}")
    print(f"  {'PS77':<12}{h1_ps:>10.5f}{h2_ps:>10.5f}")
    print(f"  {'NL94':<12}{h1_nl:>10.5f}{h2_nl:>10.5f}")
    print(f"  {'difference':<12}{h1_nl - h1_ps:>+10.5f}{h2_nl - h2_ps:>+10.5f}")
    print()

    print("How that compares with what has to be resolved.")
    print()
    print(f"  {'quantity':<38}{'h1':>10}{'h2':>10}")
    print("  " + "-" * 58)
    print(f"  {'spread between the two solar profiles':<38}{spread1:>10.4f}{spread2:>10.4f}")
    for band in ("Kepler", "TESS"):
        band_rows = [r for r in rows if r["passband"] == band]
        refld = next(r for r in band_rows if r["is_refld"])
        d1 = np.array([r["dh1"] for r in band_rows])
        d2 = np.array([r["dh2"] for r in band_rows])
        print(f"  {'signal, ' + band:<38}{abs(refld['dh1']):>10.4f}{abs(refld['dh2']):>10.4f}")
        print(f"  {'library spread, ' + band:<38}{d1.max() - d1.min():>10.4f}"
              f"{d2.max() - d2.min():>10.4f}")
    print()

    kepler = [r for r in rows if r["passband"] == "Kepler"]
    refld = next(r for r in kepler if r["is_refld"])
    print("So, in the Kepler band:")
    print(f"  the solar profiles agree to {spread1:.4f} in h1, which is "
          f"{spread1 / abs(refld['dh1']) * 100:.0f} per cent of the signal")
    print(f"  and to {spread2:.4f} in h2, which is "
          f"{spread2 / abs(refld['dh2']) * 100:.0f} per cent of the signal")
    d1 = np.array([r["dh1"] for r in kepler])
    print(f"  while the libraries disagree by {d1.max() - d1.min():.4f} in h1, "
          f"{(d1.max() - d1.min()) / spread1:.0f} times the solar spread")
    print()
    print("The Sun is therefore sharp enough. It resolves the libraries with room")
    print("to spare, and it is a fifth of the signal itself. The measurement is not")
    print("the obstacle.")
    print()
    print("An important correction to an earlier version of this note. The Sun")
    print("HAS been compared to MURaM. Norris et al. (2017) Figure 4 puts the")
    print("Neckel and Labs polynomials against the MURaM band running from the")
    print("field-free case to 100 G, at 391, 611 and 1099 nm, and their text says")
    print("the measured solar gradient is slightly shallower than the field-free")
    print("simulation, which they attribute to the quiet Sun never being truly")
    print("field free. Kostogryz et al. (2024) Figure 2 makes the same comparison")
    print("for their SSD case and finds it reproduces the measurements. So the")
    print("solar test of the MURaM models exists and it comes out in their favour.")
    print()
    print("What has not been done is ranking the OTHER libraries against the Sun.")
    print("Maxted compares nine of them to the stars. Only MPS-ATLAS and MURaM")
    print("have been compared to the Sun. The Sun cannot say whether the magnetic")
    print("interpretation is right, since Norris already showed the quiet Sun")
    print("carries its own magnetic signature, but it can say which non-magnetic")
    print("library is the right zero point, and that is what sets the size of the")
    print("inferred field.")
    print()
    print("The obstacle is that solar profiles are tabulated by wavelength and the")
    print("libraries by passband. Kostogryz et al. (2022) section 5.2 built that")
    print("bridge, a Neckel and Labs profile carried into the Kepler band called")
    print("NL_Kepler, and used it only on their own library.")
    print()
    print("Caveats on this script. The solar numbers are at 579.88 nm and the")
    print("comparison has to be done in the passbands. The polynomial fits diverge")
    print("near the limb, which matters more for h2 than h1 because h2 reaches to")
    print("mu = 1/3. Neither caveat is large enough to change the conclusion that")
    print("the solar reference is sharper than the thing it would be resolving.")


if __name__ == "__main__":
    main()
