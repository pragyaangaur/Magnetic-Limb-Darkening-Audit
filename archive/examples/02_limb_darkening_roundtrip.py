"""Measure limb darkening from an image, then use it to model an occultation.

This is the pipeline that a real measurement follows, run end to end on a
synthetic disk so that the true answer is known. A disk is rendered with chosen
limb darkening coefficients, the radial profile is measured back off the image,
a law is fitted to that profile, and the fitted coefficients are then used to
model an occultation. Comparing that model against one built from the true
coefficients shows how much a limb darkening error costs in the light curve.

With a real HMI continuum image the only change is where the profile comes
from. Everything after that step is the same.

Run it with:

    python examples/02_limb_darkening_roundtrip.py
"""

import numpy as np

from solarbatman import disk_image, occultation_flux
from solarbatman.limbdark import fit_intensity_profile

TRUTH = (0.42, 0.24)
NPIX = 2000
VENUS = 0.0311


def measured_profile(image, radius_pixels, n_bins=200):
    """Radial profile straight off a rendered image, the same way a map is read."""
    npix = image.shape[0]
    axis = np.arange(npix, dtype=np.float64) - 0.5 * (npix - 1)
    x, y = np.meshgrid(axis, axis, indexing="xy")
    r = np.hypot(x, y) / radius_pixels

    on_disk = (r < 1.0) & (image > 0)
    mu_pixels = np.sqrt(1.0 - r[on_disk] ** 2)
    values = image[on_disk]

    edges = np.linspace(0.0, 1.0, n_bins + 1)
    index = np.clip(np.digitize(mu_pixels, edges) - 1, 0, n_bins - 1)
    centres = 0.5 * (edges[:-1] + edges[1:])

    profile = np.array(
        [np.median(values[index == b]) if np.any(index == b) else np.nan
         for b in range(n_bins)]
    )
    keep = np.isfinite(profile)
    return centres[keep], profile[keep]


def main():
    print(f"Rendering a {NPIX} by {NPIX} disk with u1={TRUTH[0]}, u2={TRUTH[1]}")
    image, radius = disk_image(npix=NPIX, coeffs=TRUTH, radius_fraction=0.95)
    mu, profile = measured_profile(image, radius)
    print(f"Measured {mu.size} usable bins in mu")
    print()

    print("Fitting each law to the same measured profile:")
    header = f"{'law':<14}{'coefficients':<46}{'residual rms':>14}"
    print(header)
    print("-" * len(header))

    fits = {}
    for law in ["linear", "quadratic", "squareroot", "power2", "nonlinear"]:
        result = fit_intensity_profile(mu, profile, law=law, mu_min=0.05)
        fits[law] = result
        coeffs = ", ".join(f"{c:+.4f}" for c in result["coeffs"])
        print(f"{law:<14}{coeffs:<46}{result['residual_rms']:>14.2e}")

    print()
    recovered = fits["quadratic"]["coeffs"]
    print(f"True quadratic coefficients:      {TRUTH[0]:.4f}, {TRUTH[1]:.4f}")
    print(f"Recovered from the image:         {recovered[0]:.4f}, {recovered[1]:.4f}")
    print(f"Difference:                       {recovered[0] - TRUTH[0]:+.4f}, "
          f"{recovered[1] - TRUTH[1]:+.4f}")
    print()
    print(
        "The quadratic law is the right law here by construction, so its residual "
        "is set by the pixel grid alone. On a real image the residuals tell you "
        "which law actually describes the Sun, and the answer is usually not the "
        "quadratic one."
    )

    # What does a limb darkening error cost in the light curve?
    print()
    print("Cost of using the wrong coefficients, for a Venus sized occulter:")
    z = np.linspace(0.0, 1.05, 400)
    true_curve = occultation_flux(z, VENUS, "quadratic", TRUTH)

    for label, coeffs in [
        ("recovered from the image", recovered),
        ("off by 0.05 in u1", (TRUTH[0] + 0.05, TRUTH[1])),
        ("off by 0.10 in u1", (TRUTH[0] + 0.10, TRUTH[1])),
        ("no limb darkening at all", (0.0, 0.0)),
    ]:
        curve = occultation_flux(z, VENUS, "quadratic", coeffs)
        worst_ppm = np.abs(curve - true_curve).max() * 1e6
        print(f"  {label:<28} {worst_ppm:>8.1f} ppm at worst")

    print()
    print(
        "A Venus transit is about 1180 ppm deep, so an error of 0.10 in u1 moves "
        "the curve by a few percent of its own depth. That is the size of the "
        "systematic that solar measurements can remove from exoplanet work."
    )


if __name__ == "__main__":
    main()
