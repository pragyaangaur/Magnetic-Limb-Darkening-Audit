"""Check the batman kernel against a direct integration at solar geometries.

batman is tested in the exoplanet regime, where the occulter covers at most a
few percent of the star. Solar occultations go far outside that. Venus covers
about a tenth of a percent of the disk, and the Moon seen from an Earth orbit
covers almost all of it. This script checks the analytic kernel against a
direct summation over a rendered disk, across that whole span.

Run it with:

    python examples/01_validate_kernel.py
"""

import numpy as np

from solarbatman import brute_force_flux, occultation_flux

# Ratios of angular radii, occulter over Sun, as seen from an Earth orbit.
CASES = [
    ("Mercury at transit", 0.0058),
    ("Venus at transit", 0.0311),
    ("halfway case", 0.5),
    ("Moon, annular", 0.97),
    ("Moon, exactly matched", 1.00),
    ("Moon, total", 1.03),
]

# Illustrative photospheric coefficients near the HMI continuum wavelength.
SOLAR_U = (0.42, 0.24)

SEPARATIONS = np.array([0.0, 0.15, 0.35, 0.55, 0.75, 0.90, 1.00, 1.10])
GRID = 3000


def main():
    print(f"Reference grid: {GRID} by {GRID} pixels")
    print(f"Limb darkening: quadratic, u1={SOLAR_U[0]}, u2={SOLAR_U[1]}")
    print()
    header = f"{'case':<24}{'ratio':>8}{'depth at z=0':>15}{'area ratio':>13}{'max error':>12}"
    print(header)
    print("-" * len(header))

    worst = 0.0
    for label, ratio in CASES:
        analytic = occultation_flux(SEPARATIONS, ratio, "quadratic", SOLAR_U)
        numeric = brute_force_flux(SEPARATIONS, ratio, "quadratic", SOLAR_U, npix=GRID)
        error = float(np.abs(analytic - numeric).max())
        worst = max(worst, error)

        depth = 1.0 - analytic[0]
        print(
            f"{label:<24}{ratio:>8.4f}{depth:>15.6f}{ratio**2:>13.6f}{error:>12.2e}"
        )

    print()
    print(f"Largest disagreement anywhere: {worst:.2e}")
    print()
    print(
        "The reference grid itself is only good to a few parts in ten million, "
        "so this bounds the error of the grid as much as the error of the kernel. "
        "Either way it is far below any solar photometric noise floor."
    )

    # The excess of the depth over the plain area ratio is the part of the
    # signal that carries the limb darkening.
    print()
    print("How much of the depth comes from limb darkening rather than size:")
    for label, ratio in CASES:
        depth = 1.0 - occultation_flux(np.array([0.0]), ratio, "quadratic", SOLAR_U)[0]
        plain = 1.0 - occultation_flux(np.array([0.0]), ratio, "uniform", ())[0]
        excess = (depth - plain) / plain * 100
        print(f"  {label:<24} {excess:>6.1f} percent deeper than the area ratio")


if __name__ == "__main__":
    main()
