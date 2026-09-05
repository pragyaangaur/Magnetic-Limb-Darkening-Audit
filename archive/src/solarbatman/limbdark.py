"""Limb darkening laws, and fitting them to a measured solar profile.

The laws here are the same ones batman uses, written as intensity against
``mu``, the cosine of the angle between the line of sight and the local surface
normal. On a disk of unit radius, ``mu = sqrt(1 - r**2)`` where ``r`` is the
distance from disk centre. The disk centre has ``mu = 1`` and the limb has
``mu = 0``.

Most of the laws give ``I(mu=1) = 1`` by construction. The exponential law is
the exception, because its second term does not vanish at disk centre. batman
normalises internally by the integral of the intensity over the disk, so the
overall scale never reaches the occultation flux and only the shape matters.

One warning about the exponential law. Its second term goes as one over
``1 - exp(mu)``, which has a pole at ``mu = 0``. The intensity therefore
diverges at the limb rather than falling to a finite value. In exoplanet work
the occulter is small and crosses the limb quickly, so this rarely shows up. A
solar occulter is large and spends a long time on the limb, and there the
divergence produces a light curve that is not even monotonic. Use the
quadratic, square root, power-2 or nonlinear laws for solar work.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from .kernel import LIMB_DARKENING_LAWS, n_coefficients

__all__ = [
    "intensity",
    "fit_intensity_profile",
    "SOLAR_REFERENCE",
]

# Rough quadratic coefficients for the solar photosphere, for orientation only.
# They come from the shape of published solar limb darkening curves and are
# meant as starting guesses. Measure your own from an HMI continuum image
# before quoting a number anywhere.
SOLAR_REFERENCE = {
    # wavelength in nanometres: (u1, u2)
    400.0: (0.72, 0.14),
    500.0: (0.55, 0.21),
    617.3: (0.42, 0.24),  # HMI continuum
    700.0: (0.36, 0.25),
    1600.0: (0.16, 0.22),
}


def intensity(mu, law="quadratic", coeffs=()):
    """
    Intensity of a limb darkened disk at a given ``mu``.

    Parameters
    ----------
    mu : array_like
        Cosine of the angle from the surface normal. One at disk centre, zero
        at the limb.
    law : str, optional
        Limb darkening law, matching the names in
        `solarbatman.kernel.LIMB_DARKENING_LAWS`.
    coeffs : array_like, optional
        Coefficients for the law.

    Returns
    -------
    numpy.ndarray
        Intensity normalised to one at disk centre.

    Examples
    --------
    >>> import numpy as np
    >>> float(np.round(intensity(0.0, "quadratic", (0.42, 0.24)), 4))
    0.34
    """
    mu = np.asarray(mu, dtype=np.float64)
    expected = n_coefficients(law)
    c = np.atleast_1d(np.asarray(coeffs, dtype=np.float64)).ravel()
    if c.size != expected:
        raise ValueError(f"law {law!r} takes {expected} coefficients, got {c.size}")

    if law == "uniform":
        return np.ones_like(mu)
    if law == "linear":
        return 1 - c[0] * (1 - mu)
    if law == "quadratic":
        return 1 - c[0] * (1 - mu) - c[1] * (1 - mu) ** 2
    if law == "squareroot":
        return 1 - c[0] * (1 - mu) - c[1] * (1 - np.sqrt(mu))
    if law == "logarithmic":
        # The log term diverges at mu = 0, so the limb is clipped. This matches
        # what any numerical evaluation of the law has to do.
        safe = np.clip(mu, 1e-12, None)
        return 1 - c[0] * (1 - safe) - c[1] * safe * np.log(safe)
    if law == "exponential":
        # This law has a pole at mu = 0. The intensity runs away at the very
        # limb instead of settling at a finite value, which makes it a poor
        # choice for occultations that spend real time crossing the limb.
        safe = np.clip(mu, 1e-12, None)
        return 1 - c[0] * (1 - safe) - c[1] / (1 - np.exp(safe))
    if law == "power2":
        safe = np.clip(mu, 1e-12, None)
        return 1 - c[0] * (1 - safe ** c[1])
    if law == "nonlinear":
        safe = np.clip(mu, 0.0, None)
        return (
            1
            - c[0] * (1 - safe**0.5)
            - c[1] * (1 - safe)
            - c[2] * (1 - safe**1.5)
            - c[3] * (1 - safe**2)
        )
    raise ValueError(f"unhandled limb darkening law {law!r}")


def _initial_guess(law):
    return {
        "linear": [0.5],
        "quadratic": [0.4, 0.2],
        "squareroot": [0.3, 0.3],
        "logarithmic": [0.5, 0.2],
        "exponential": [0.5, 0.05],
        "power2": [0.5, 0.6],
        "nonlinear": [0.5, 0.1, 0.1, -0.1],
    }.get(law, [0.3] * n_coefficients(law))


def fit_intensity_profile(mu, flux, law="quadratic", mu_min=0.05, normalise=True):
    """
    Fit a limb darkening law to a measured intensity profile.

    Parameters
    ----------
    mu : array_like
        Cosine of the angle from the surface normal for each sample.
    flux : array_like
        Measured intensity at each ``mu``. Any overall scale is allowed when
        ``normalise`` is true.
    law : str, optional
        Law to fit.
    mu_min : float, optional
        Samples below this ``mu`` are dropped. The extreme limb is where
        seeing, scattered light and the instrument point spread function do the
        most damage, and where the laws themselves fit worst.
    normalise : bool, optional
        When true the profile is rescaled so that the fitted intensity is one
        at disk centre. Turn this off if the input is already normalised.

    Returns
    -------
    dict
        With keys ``coeffs``, ``scale``, ``residual_rms``, ``law`` and
        ``n_points``.

    Notes
    -----
    The fit is a plain least squares fit with no weighting. Solar profiles have
    very many points and very small formal errors, so the residual scatter is
    dominated by how well the law describes the Sun rather than by photon
    noise. Look at ``residual_rms`` to compare laws against each other.
    """
    if law not in LIMB_DARKENING_LAWS:
        raise ValueError(f"unknown law {law!r}")

    mu = np.asarray(mu, dtype=np.float64).ravel()
    flux = np.asarray(flux, dtype=np.float64).ravel()
    if mu.shape != flux.shape:
        raise ValueError("mu and flux must have the same shape")

    keep = np.isfinite(mu) & np.isfinite(flux) & (mu >= mu_min)
    mu, flux = mu[keep], flux[keep]
    if mu.size < n_coefficients(law) + 1:
        raise ValueError("not enough finite samples above mu_min to fit")

    scale0 = float(np.median(flux[mu > 0.9])) if np.any(mu > 0.9) else float(flux.max())
    if scale0 == 0:
        raise ValueError("profile has zero intensity near disk centre")

    def residual(p):
        if normalise:
            scale, c = p[0], p[1:]
        else:
            scale, c = 1.0, p
        return scale * intensity(mu, law, c) - flux

    p0 = [scale0, *_initial_guess(law)] if normalise else _initial_guess(law)
    result = least_squares(residual, p0, method="lm", max_nfev=20000)

    scale = float(result.x[0]) if normalise else 1.0
    coeffs = np.asarray(result.x[1:] if normalise else result.x, dtype=np.float64)
    rms = float(np.sqrt(np.mean(result.fun**2)) / scale)

    return {
        "law": law,
        "coeffs": coeffs,
        "scale": scale,
        "residual_rms": rms,
        "n_points": int(mu.size),
    }
