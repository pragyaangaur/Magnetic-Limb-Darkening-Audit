"""Turn a limb darkening error into a planet radius error.

Kostogryz et al. (2024, Nature Astronomy) showed that limb darkening from
non-magnetic stellar atmosphere models is systematically wrong, and that
small-scale surface magnetic fields explain the difference. Their follow up
(arXiv:2606.21912, June 2026) quantifies how the effect varies across the main
sequence and releases the intensity tables.

Neither paper converts that limb darkening error into an error on the retrieved
planet radius. This script does that conversion, so the size of the consequence
can be read off directly.

The experiment. A transit is generated with the limb darkening the star really
has. It is then fitted with the limb darkening a modeller would have used,
under the two strategies people actually follow:

  fixed  the coefficients are held at the wrong tabulated values, which is
         standard for shallow transits and low signal to noise data
  free   the coefficients float, which is standard for deep, well sampled
         transits

The output is the fractional error in the recovered radius ratio.

Run it with:

    python examples/07_limb_darkening_error_to_radius_bias.py
"""

import numpy as np
from scipy.optimize import least_squares

from solarbatman import occultation_flux
from solarbatman.geometry import linear_track

# Limb darkening the star really has. A solar type star in an optical band.
TRUE_U = (0.42, 0.24)

# Errors in the tabulated coefficients. The first two are the scale by which
# non-magnetic models are reported to miss. The third is the offset reported
# for PHOENIX coefficients against empirical values from TESS.
ERRORS = [
    ("small, du1=+0.03", (0.03, 0.0)),
    ("reported scale, du1=+0.05", (0.05, 0.0)),
    ("reported scale, du1=+0.10", (0.10, 0.0)),
    ("PHOENIX offset, du2=-0.20", (0.0, -0.20)),
    ("combined, du1=+0.08 du2=-0.10", (0.08, -0.10)),
]

PLANETS = [
    ("Earth around a G star", 0.009),
    ("mini Neptune", 0.030),
    ("Neptune", 0.055),
    ("hot Jupiter", 0.110),
]

IMPACTS = [0.0, 0.3, 0.6, 0.8]
N_POINTS = 1500


def light_curve(ratio, impact, coeffs):
    t = np.linspace(-1.5, 1.5, N_POINTS)
    z = linear_track(t, impact_parameter=impact)
    return z, occultation_flux(z, ratio, "quadratic", coeffs)


def fit_fixed_ld(z, observed, ratio_guess, coeffs):
    """Fit the radius ratio and a free baseline, with limb darkening held wrong."""

    def residual(p):
        return p[1] * occultation_flux(z, abs(p[0]), "quadratic", coeffs) - observed

    fit = least_squares(residual, [ratio_guess, 1.0], method="lm", max_nfev=20000)
    return abs(fit.x[0]), float(np.sqrt(np.mean(fit.fun**2)))


def fit_free_ld(z, observed, ratio_guess, coeffs):
    """Fit the radius ratio, both limb darkening coefficients and a baseline."""

    def residual(p):
        return (
            p[3] * occultation_flux(z, abs(p[0]), "quadratic", p[1:3]) - observed
        )

    start = [ratio_guess, coeffs[0], coeffs[1], 1.0]
    fit = least_squares(residual, start, method="lm", max_nfev=20000)
    return abs(fit.x[0]), float(np.sqrt(np.mean(fit.fun**2)))


def main():
    print(f"True limb darkening: u1 = {TRUE_U[0]:.2f}, u2 = {TRUE_U[1]:.2f}")
    print("The transit is generated with the truth and fitted with the error.")
    print("Numbers are the percentage error in the recovered radius ratio.")
    print()

    for label, (du1, du2) in ERRORS:
        wrong = (TRUE_U[0] + du1, TRUE_U[1] + du2)
        print(f"=== {label}   (fitted with u1={wrong[0]:.2f}, u2={wrong[1]:.2f}) ===")

        header = f"  {'planet':<24}{'b':>6}" + "".join(
            f"{h:>16}" for h in ("LD fixed", "LD free", "residual, fixed")
        )
        print(header)
        print("  " + "-" * (len(header) - 2))

        for planet_label, ratio in PLANETS:
            for b in IMPACTS:
                z, observed = light_curve(ratio, b, TRUE_U)
                fixed, resid_fixed = fit_fixed_ld(z, observed, ratio, wrong)
                free, _ = fit_free_ld(z, observed, ratio, wrong)

                e_fixed = (fixed - ratio) / ratio * 100
                e_free = (free - ratio) / ratio * 100
                name = planet_label if b == IMPACTS[0] else ""
                print(
                    f"  {name:<24}{b:>6.1f}{e_fixed:>15.3f} %{e_free:>15.3f} %"
                    f"{resid_fixed * 1e6:>13.1f} ppm"
                )
        print()

    print("Reading this.")
    print()
    print(
        "The radius error under fixed limb darkening does not depend much on "
        "planet size, because it is a fractional error on the depth rather than "
        "an absolute one. It does depend strongly on impact parameter, because "
        "that sets which part of the stellar disk the planet samples."
    )
    print()
    print(
        "Letting the coefficients float removes most of the bias for a deep "
        "transit. It cannot be done for a shallow one, which is exactly the "
        "regime where small planets live and where radii feed into occurrence "
        "rates and into the radius valley."
    )
    print()
    print(
        "The residual column is the warning signal. Compare it against the per "
        "point precision of the survey. TESS reaches a few hundred parts per "
        "million per two minute cadence for a bright star, and Kepler reached "
        "tens. A residual below that is invisible."
    )


if __name__ == "__main__":
    main()
