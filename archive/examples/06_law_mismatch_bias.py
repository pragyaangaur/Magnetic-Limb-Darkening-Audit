"""Do image derived and transit derived limb darkening agree?

Every transit fit describes limb darkening with a parametric law, and the
quadratic law is by far the most used. No real star follows any of these laws
exactly. So there are two different ways to get quadratic coefficients for the
same star, and they need not give the same answer.

The first way is to measure the intensity profile directly off a resolved image
and least squares fit the quadratic law to it. This is only possible for the
Sun.

The second way is to fit the quadratic law to a transit light curve. This is
what is done for every other star.

Both are least squares fits of the same wrong law to the same star, but they
minimise different things. The image fit weights every part of the disk by how
many pixels it covers. The transit fit weights it by how much the occulter
happens to cross, which depends on the impact parameter. If the law were exact,
both would return the same coefficients. It is not exact, so they do not.

This script measures the size of that disagreement, and the radius error that
follows from it. Truth is taken to be the four parameter nonlinear law, which
is the form used to describe model stellar atmospheres, and the fitted law is
quadratic.

Run it with:

    python examples/06_law_mismatch_bias.py
"""

import numpy as np
from scipy.optimize import least_squares

from solarbatman import occultation_flux
from solarbatman.geometry import linear_track
from solarbatman.limbdark import intensity

# Truth profiles. These are four parameter nonlinear coefficients of the kind
# tabulated for a solar type star in an optical band. They are stand-ins for a
# real tabulated set and are used only to have a truth that is not itself a
# quadratic.
TRUTHS = {
    "nonlinear A": ("nonlinear", (0.55, -0.10, 0.55, -0.25)),
    "nonlinear B": ("nonlinear", (0.70, -0.35, 0.85, -0.35)),
    "power2": ("power2", (0.52, 0.68)),
}

RATIO = 0.0311  # Venus sized, a typical hot Jupiter is a few times larger
IMPACTS = [0.0, 0.3, 0.5, 0.7, 0.85]
MU_MIN = 0.05


def image_quadratic(truth_law, truth_coeffs, area_weighted=True):
    """
    Quadratic coefficients from a resolved image of the disk.

    On a disk of unit radius the number of pixels between mu and mu + dmu goes
    as mu dmu, because the area element is r dr and r dr equals mu dmu. So an
    unbinned pixel by pixel fit carries that weight. A fit to a binned profile
    carries no weight unless one is applied, and the two answers differ.
    """
    mu = np.linspace(MU_MIN, 1.0, 4000)
    target = intensity(mu, truth_law, truth_coeffs)
    weight = np.sqrt(mu) if area_weighted else np.ones_like(mu)

    def residual(p):
        return (intensity(mu, "quadratic", p) - target) * weight

    return least_squares(residual, [0.4, 0.2], method="lm").x


def transit_quadratic(truth_law, truth_coeffs, ratio, impact):
    """Quadratic coefficients and radius ratio from fitting a transit light curve."""
    t = np.linspace(-1.4, 1.4, 900)
    z = linear_track(t, impact_parameter=impact)

    # The exact light curve of a star that follows the truth law.
    observed = occultation_flux(z, ratio, truth_law, truth_coeffs)

    def residual(p):
        model = occultation_flux(z, abs(p[0]), "quadratic", p[1:3])
        return model - observed

    fit = least_squares(residual, [ratio, 0.4, 0.2], method="lm", max_nfev=20000)
    return abs(fit.x[0]), fit.x[1:3], float(np.sqrt(np.mean(fit.fun**2)))


def main():
    print("Fitted law: quadratic. Occulter ratio:", RATIO)
    print("Truth is a law the quadratic cannot reproduce exactly, which is the")
    print("situation for every real star.")
    print()

    for name, (truth_law, truth_coeffs) in TRUTHS.items():
        u_image_w = image_quadratic(truth_law, truth_coeffs, area_weighted=True)
        u_image_u = image_quadratic(truth_law, truth_coeffs, area_weighted=False)

        print(f"=== truth: {name} ({truth_law}) ===")
        print(f"  image derived, area weighted   u1={u_image_w[0]:+.4f}  u2={u_image_w[1]:+.4f}")
        print(f"  image derived, unweighted      u1={u_image_u[0]:+.4f}  u2={u_image_u[1]:+.4f}")
        print(
            f"  the weighting choice alone moves u1 by "
            f"{u_image_u[0] - u_image_w[0]:+.4f}"
        )
        print()

        header = (
            f"  {'impact':>7}{'u1 transit':>13}{'u2 transit':>13}"
            f"{'u1 - image':>13}{'u2 - image':>13}{'radius error':>15}{'residual':>12}"
        )
        print(header)
        print("  " + "-" * (len(header) - 2))

        for b in IMPACTS:
            ratio_fit, u_fit, resid = transit_quadratic(
                truth_law, truth_coeffs, RATIO, b
            )
            d_ratio = (ratio_fit - RATIO) / RATIO * 100
            print(
                f"  {b:>7.2f}{u_fit[0]:>13.4f}{u_fit[1]:>13.4f}"
                f"{u_fit[0] - u_image_w[0]:>13.4f}{u_fit[1] - u_image_w[1]:>13.4f}"
                f"{d_ratio:>13.3f} %{resid:>12.2e}"
            )
        print()

    print("What to look at.")
    print()
    print(
        "The 'u1 - image' column is the disagreement between the two ways of "
        "measuring the same star. It is not noise. Both fits are noiseless and "
        "both converge. The disagreement exists because the quadratic law "
        "cannot describe the star, so the answer depends on which part of the "
        "disk the measurement weights."
    )
    print()
    print(
        "The 'radius error' column is what that costs. The occulter has one "
        "true size in every row, and the recovered size changes with impact "
        "parameter alone. Two planets of identical size around identical stars "
        "would be reported at different radii purely because they cross at "
        "different latitudes."
    )
    print()
    print(
        "The residual column says whether the light curve gives any warning. "
        "Compare it to the photometric precision of a real transit survey."
    )


if __name__ == "__main__":
    main()
