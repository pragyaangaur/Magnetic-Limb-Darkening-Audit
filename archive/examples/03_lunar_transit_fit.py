"""Forward model a lunar transit of the Sun, then fit it back.

SDO passes through the Moon's shadow several times a year and keeps observing
while the Moon crosses the solar disk. Seen from a geosynchronous orbit the
Moon and the Sun have almost the same angular size, so the crossing removes
almost all of the light. This script builds that light curve, adds noise at a
level a real disk integrated light curve might reach, fits it back and writes
the summary figure used in the README.

The geometry here is a straight track rather than a real ephemeris, because the
point is to test the fitter. Swap in
``solarbatman.geometry.separation_from_ephemeris`` to run it on a real event.

Run it with:

    python examples/03_lunar_transit_fit.py
"""

import pathlib

import numpy as np

from solarbatman import brute_force_flux, occultation_flux
from solarbatman.fitting import fit_occultation
from solarbatman.geometry import linear_track
from solarbatman.limbdark import intensity

TRUE_RATIO = 0.985  # angular radius ratio, a near total lunar transit
TRUE_U = (0.42, 0.24)
IMPACT = 0.25
NOISE = 2e-3  # a hopeful but not absurd disk integrated photometric noise floor
OUTPUT = pathlib.Path(__file__).resolve().parents[1] / "docs" / "img"


def build_curve(seed=20260903):
    hours = np.linspace(-2.2, 2.2, 660)  # 12 second cadence would give far more
    z = linear_track(hours, t0=0.0, impact_parameter=IMPACT, crossing_time=1.0)
    clean = occultation_flux(z, TRUE_RATIO, "quadratic", TRUE_U)
    noisy = clean + np.random.default_rng(seed).normal(0, NOISE, clean.size)
    return hours, z, clean, noisy


def main():
    hours, z, clean, noisy = build_curve()

    print(f"True ratio            {TRUE_RATIO}")
    print(f"True limb darkening   u1={TRUE_U[0]}, u2={TRUE_U[1]}")
    print(f"Impact parameter      {IMPACT}")
    print(f"Noise per sample      {NOISE:.1e}")
    print(f"Samples               {hours.size}")
    print(f"Deepest point         {1 - clean.min():.4f} of the light removed")
    print()

    sigma = np.full(noisy.size, NOISE)

    print("Fit A, everything free. This is the exoplanet situation.")
    a = fit_occultation(z, noisy, ratio_guess=0.9, coeffs_guess=(0.3, 0.3), uncertainty=sigma)
    print(f"  ratio       {a['ratio']:.5f}  (true {TRUE_RATIO}, error "
          f"{a['ratio'] - TRUE_RATIO:+.5f})")
    print(f"  u1, u2      {a['coeffs'][0]:.4f}, {a['coeffs'][1]:.4f}")
    print(f"  residual    {a['residual_rms']:.2e}")
    print()

    print("Fit B, ratio fixed by the ephemeris. This is the solar situation.")
    b = fit_occultation(
        z, noisy, ratio_guess=TRUE_RATIO, fit_ratio=False,
        coeffs_guess=(0.3, 0.3), uncertainty=sigma,
    )
    print(f"  u1          {b['coeffs'][0]:.4f}  (true {TRUE_U[0]}, error "
          f"{b['coeffs'][0] - TRUE_U[0]:+.4f})")
    print(f"  u2          {b['coeffs'][1]:.4f}  (true {TRUE_U[1]}, error "
          f"{b['coeffs'][1] - TRUE_U[1]:+.4f})")
    print(f"  residual    {b['residual_rms']:.2e}")
    if b["errors"]:
        print(f"  formal error on u1  {b['errors']['c0']:.4f}")
    print()
    print(
        "Fixing the ratio removes one free parameter that is strongly correlated "
        "with the limb darkening, so the limb darkening comes back tighter. That "
        "is the whole reason to do this on the Sun rather than on a star."
    )

    make_figure(hours, z, clean, noisy, b)


def make_figure(hours, z, clean, noisy, fit):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nmatplotlib is not installed, so no figure was written.")
        return

    OUTPUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.2))
    fig.suptitle(
        "batman occultation kernels applied to solar geometries", fontsize=13
    )

    # Panel 1: light curves across the whole range of solar occulters.
    ax = axes[0, 0]
    grid = np.linspace(0.0, 1.15, 500)
    for label, ratio in [("Mercury, 0.0058", 0.0058), ("Venus, 0.0311", 0.0311),
                         ("Moon, 0.985", 0.985)]:
        depth = 1.0 - occultation_flux(grid, ratio, "quadratic", TRUE_U)
        ax.semilogy(grid, np.clip(depth, 1e-7, None), label=label)
    ax.set_xlabel("projected separation (solar radii)")
    ax.set_ylabel("fraction of light removed")
    ax.set_title("Occultation depth against separation")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)

    # Panel 2: the kernel against a direct integration.
    ax = axes[0, 1]
    probe = np.linspace(0.0, 1.1, 12)
    for ratio in (0.0311, 0.5, 0.97):
        analytic = occultation_flux(probe, ratio, "quadratic", TRUE_U)
        numeric = brute_force_flux(probe, ratio, "quadratic", TRUE_U, npix=2000)
        ax.plot(probe, (analytic - numeric) * 1e6, "o-", ms=3,
                label=f"ratio {ratio}")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("projected separation (solar radii)")
    ax.set_ylabel("difference (ppm), grid limited")
    ax.set_title("Kernel minus a 2000 pixel direct integration")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)

    # Panel 3: the limb darkening laws that were fitted.
    ax = axes[1, 0]
    mu = np.linspace(0.0, 1.0, 300)
    ax.plot(mu, intensity(mu, "quadratic", TRUE_U), label="true profile", lw=2)
    ax.plot(mu, intensity(mu, "quadratic", fit["coeffs"]), "--",
            label="recovered from the fit")
    ax.plot(mu, intensity(mu, "uniform", ()), ":", color="grey",
            label="no limb darkening")
    ax.set_xlabel("mu, cosine of the angle from the normal")
    ax.set_ylabel("relative intensity")
    ax.set_title("Limb darkening recovered from the light curve")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)

    # Panel 4: the lunar transit fit and its residuals.
    ax = axes[1, 1]
    ax.plot(hours, noisy, ".", ms=2, color="0.6", label="simulated data")
    ax.plot(hours, fit["model"], lw=1.6, color="C3", label="fitted model")
    ax.plot(hours, (noisy - fit["model"]) * 10 + 0.05, ".", ms=2, color="C0",
            label="residuals, 10x, offset")
    ax.axhline(0.05, color="k", lw=0.5)
    ax.set_xlabel("hours from mid transit")
    ax.set_ylabel("relative flux")
    ax.set_title("Lunar transit, fitted with the ratio fixed")
    ax.legend(fontsize=8, loc="center right")
    ax.grid(alpha=0.25)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = OUTPUT / "overview.png"
    fig.savefig(path, dpi=130)
    print(f"\nFigure written to {path}")


if __name__ == "__main__":
    main()
