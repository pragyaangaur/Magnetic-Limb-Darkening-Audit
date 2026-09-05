"""How much does a sunspot bias a transit fit?

This is the question the Sun can answer and other stars cannot. Stellar
activity biases every transit measurement, and the correction is normally made
with assumptions about how much of the star is covered in spots. On the Sun the
spots are visible in the same images that produce the light curve, so the bias
can be attributed to a specific spot rather than estimated.

Here that is done on a rendered disk, where the truth is known exactly. A disk
is drawn with a spot on it, an occulter is walked across, and the light curve is
built by direct summation over the pixels. That curve is then fitted with the
clean model, which has no spot in it, and the error in the recovered parameters
is the bias.

The occulter is Venus sized, because that is the case where the limb darkening
signal is strongest and where a spot does the most damage relative to the depth.

Run it with:

    python examples/05_spot_crossing_bias.py
"""

import numpy as np

from solarbatman import add_spot, disk_image, occulted_disk_image
from solarbatman.fitting import fit_occultation
from solarbatman.geometry import linear_track

NPIX = 1400
RATIO = 0.0311  # Venus at transit
TRUE_U = (0.42, 0.24)
IMPACT = 0.0  # a central crossing, so the occulter passes over disk centre
N_SAMPLES = 90


def curve_over_image(image, radius, ratio, z_values, angle=0.0):
    """Occultation light curve by direct summation, so a spot is included."""
    total = image.sum()
    return np.array(
        [
            occulted_disk_image(image, radius, z=float(z), ratio=ratio, angle=angle).sum()
            / total
            for z in z_values
        ]
    )


def main():
    hours = np.linspace(-1.15, 1.15, N_SAMPLES)
    z = linear_track(hours, impact_parameter=IMPACT, crossing_time=1.0)

    clean, radius = disk_image(npix=NPIX, coeffs=TRUE_U, radius_fraction=0.98)
    baseline_curve = curve_over_image(clean, radius, RATIO, z)

    print(f"Grid                {NPIX} by {NPIX}")
    print(f"Occulter ratio      {RATIO} (Venus at transit)")
    print(f"Transit depth       {(1 - baseline_curve.min()) * 1e6:.0f} ppm")
    print(f"Samples             {N_SAMPLES}")
    print()

    # A reference fit to the unspotted curve. Any error here is the pixel grid
    # rather than the spot, so it is the floor everything else is measured
    # against.
    reference = fit_occultation(
        z, baseline_curve, ratio_guess=RATIO, coeffs_guess=(0.4, 0.2)
    )
    print("Unspotted disk, which measures the pixel grid rather than any physics:")
    report(reference, reference)
    print()

    print("Now with a spot. The occulter crosses the disk along y = 0, so a spot")
    print("at y = 0 is crossed directly and a spot at larger y is never covered.")
    print()

    header = (
        f"{'spot':<34}{'ratio error':>13}{'u1 error':>11}"
        f"{'u2 error':>11}{'residual':>12}"
    )
    print(header)
    print("-" * len(header))

    cases = [
        ("no spot", None),
        ("crossed, r=0.10, contrast 0.3", (0.0, 0.0, 0.10, 0.3)),
        ("crossed, r=0.05, contrast 0.3", (0.0, 0.0, 0.05, 0.3)),
        ("crossed, r=0.05, contrast 0.7", (0.0, 0.0, 0.05, 0.7)),
        ("off centre, r=0.10 at x=0.5", (0.5, 0.0, 0.10, 0.3)),
        ("never crossed, r=0.10 at y=0.5", (0.0, 0.5, 0.10, 0.3)),
        ("never crossed, r=0.20 at y=0.5", (0.0, 0.5, 0.20, 0.3)),
    ]

    for label, spot in cases:
        if spot is None:
            image = clean
        else:
            x, y, r, contrast = spot
            image = add_spot(clean, radius, x, y, r, contrast=contrast)

        curve = curve_over_image(image, radius, RATIO, z)
        fit = fit_occultation(z, curve, ratio_guess=RATIO, coeffs_guess=(0.4, 0.2))

        d_ratio = (fit["ratio"] - RATIO) / RATIO * 100
        d_u1 = fit["coeffs"][0] - TRUE_U[0]
        d_u2 = fit["coeffs"][1] - TRUE_U[1]
        print(
            f"{label:<34}{d_ratio:>11.2f} %{d_u1:>11.4f}{d_u2:>11.4f}"
            f"{fit['residual_rms']:>12.2e}"
        )

    print()
    print("Reading the table.")
    print()
    print(
        "A spot the occulter crosses pushes the recovered radius ratio down, "
        "because the occulter spends part of its time covering something "
        "already dark, so it appears to remove less light than it should. It "
        "wrecks the limb darkening completely. The coefficients come back off "
        "by more than one, which is far outside any physical range, because the "
        "fit is bending the intensity profile to absorb a bump that is not a "
        "profile at all."
    )
    print()
    print(
        "A spot the occulter never crosses pushes the ratio the other way. The "
        "spot dims the whole disk, so the same absolute loss of light becomes a "
        "larger fraction of a smaller total, and the transit looks deeper than "
        "it is. A free baseline does not help, because the curve is already "
        "normalised to the dimmed star. This is the known systematic that "
        "inflates measured exoplanet radii on spotted stars, and here it is "
        "with a number on it."
    )
    print()
    print(
        "The residual column separates the two cases and this is the useful "
        "part. A crossed spot raises the residual by two orders of magnitude, "
        "so it announces itself and can be cut. An uncrossed spot leaves the "
        "residual at the noise floor while still biasing the radius by over one "
        "percent. That is the dangerous case, because nothing in the light "
        "curve reveals it. On the Sun it is visible in the image."
    )


def report(fit, _reference):
    print(f"  ratio     {fit['ratio']:.6f}  (true {RATIO}, "
          f"error {(fit['ratio'] - RATIO) / RATIO * 100:+.3f} percent)")
    print(f"  u1, u2    {fit['coeffs'][0]:.4f}, {fit['coeffs'][1]:.4f}  "
          f"(true {TRUE_U[0]}, {TRUE_U[1]})")
    print(f"  residual  {fit['residual_rms']:.2e}")


if __name__ == "__main__":
    main()
