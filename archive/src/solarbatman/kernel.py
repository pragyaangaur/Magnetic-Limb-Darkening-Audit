"""Direct access to the batman occultation kernels.

``batman.TransitModel`` takes Keplerian orbital elements and turns them into a
sky projected separation between the two disk centres. That separation is then
handed to a C kernel which returns the occulted flux. The orbit is only used to
produce the separation, so any occultation geometry can use the same kernels as
long as the separation is supplied directly.

Solar occultations need exactly that. The Moon seen from SDO does not follow an
orbit that ``TransitModel`` can express, but its projected separation from the
solar disk centre is easy to get from an ephemeris. This module exposes the
kernels so that separation can be passed straight in.

All separations are in units of the occulted disk radius, which for solar work
means solar radii. All radius ratios are occulter radius over solar radius.
"""

from __future__ import annotations

import numpy as np
from batman import (
    _exponential_ld,
    _logarithmic_ld,
    _nonlinear_ld,
    _power2_ld,
    _quadratic_ld,
    _uniform_ld,
)

__all__ = [
    "LIMB_DARKENING_LAWS",
    "n_coefficients",
    "occultation_flux",
    "choose_step_size",
]

# Number of limb darkening coefficients each law takes. The names match the
# strings accepted by batman.TransitModel so that code can move between the two.
LIMB_DARKENING_LAWS = {
    "uniform": 0,
    "linear": 1,
    "quadratic": 2,
    "squareroot": 2,
    "logarithmic": 2,
    "exponential": 2,
    "power2": 2,
    "nonlinear": 4,
}

# Laws that batman evaluates by numerical integration rather than analytically.
# These take a step size argument, the others do not.
_NUMERICAL_LAWS = frozenset(
    {"squareroot", "logarithmic", "exponential", "power2", "nonlinear"}
)

# Smallest step size batman itself will consider. Going below this buys no
# accuracy and costs a lot of time.
_STEP_FLOOR = 5.0e-4


def n_coefficients(law):
    """
    Return the number of limb darkening coefficients a law takes.

    Parameters
    ----------
    law : str
        One of the keys of `LIMB_DARKENING_LAWS`.

    Returns
    -------
    int

    Examples
    --------
    >>> n_coefficients("quadratic")
    2
    >>> n_coefficients("nonlinear")
    4
    """
    try:
        return LIMB_DARKENING_LAWS[law]
    except KeyError:
        raise ValueError(
            f"unknown limb darkening law {law!r}, expected one of "
            f"{sorted(LIMB_DARKENING_LAWS)}"
        ) from None


def _check_inputs(z, ratio, law, coeffs):
    z = np.ascontiguousarray(z, dtype=np.float64)
    if np.any(z < 0):
        raise ValueError("projected separation must be non-negative")
    if not np.all(np.isfinite(z)):
        raise ValueError("projected separation contains non-finite values")

    ratio = float(ratio)
    if ratio < 0:
        raise ValueError("radius ratio must be non-negative")

    expected = n_coefficients(law)
    coeffs = np.atleast_1d(np.asarray(coeffs, dtype=np.float64)).ravel()
    if coeffs.size != expected:
        raise ValueError(
            f"limb darkening law {law!r} takes {expected} coefficients, "
            f"got {coeffs.size}"
        )
    return z, ratio, coeffs


def choose_step_size(ratio, law, coeffs, max_error_ppm=1.0, nthreads=1):
    """
    Pick an integration step size for the numerically integrated laws.

    batman evaluates the square root, logarithmic, exponential, power-2 and
    nonlinear laws by integrating over the occulted region. The step size sets
    the accuracy. This function bisects for the largest step that stays within
    ``max_error_ppm`` of the most accurate step batman supports, which is the
    same approach ``batman.TransitModel`` uses internally.

    Parameters
    ----------
    ratio : float
        Occulter radius over solar radius.
    law : str
        Limb darkening law.
    coeffs : array_like
        Limb darkening coefficients.
    max_error_ppm : float, optional
        Tolerated error in parts per million of the unocculted flux.
    nthreads : int, optional
        Threads to use during the search.

    Returns
    -------
    float
        Step size to pass as ``step`` to `occultation_flux`.

    Notes
    -----
    A lunar occulter covers most of the solar disk, so the integrated region
    is far larger than in any exoplanet case and the step size matters more.
    Call this once per geometry and reuse the result across a fit.
    """
    if law not in _NUMERICAL_LAWS:
        return _STEP_FLOOR

    probe = np.linspace(0.0, 1.0 + ratio, 1000)
    reference = occultation_flux(
        probe, ratio, law, coeffs, step=_STEP_FLOOR, nthreads=nthreads
    )

    lo, hi = _STEP_FLOOR, 1.0
    for _ in range(64):
        mid = 0.5 * (lo + hi)
        trial = occultation_flux(
            probe, ratio, law, coeffs, step=mid, nthreads=nthreads
        )
        error_ppm = np.abs(trial - reference).max() * 1e6
        if error_ppm > max_error_ppm:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-6:
            break
    return lo


def occultation_flux(z, ratio, law="quadratic", coeffs=(), step=None, nthreads=1):
    """
    Flux of a limb darkened disk partly covered by an opaque circle.

    Parameters
    ----------
    z : array_like
        Sky projected separation between the centre of the occulted disk and
        the centre of the occulter, in units of the occulted disk radius. For
        solar work this is in solar radii. Values above ``1 + ratio`` give an
        unocculted flux of one.
    ratio : float
        Ratio of angular radii, occulter over occulted disk. Seen from an
        Earth orbit this is about 0.97 to 1.03 for the Moon, 0.031 for Venus
        at transit and 0.006 for Mercury at transit. The physical radius ratio
        is a different and much smaller number, and it is not what goes here.
    law : str, optional
        Limb darkening law. See `LIMB_DARKENING_LAWS`.
    coeffs : array_like, optional
        Limb darkening coefficients for the chosen law.
    step : float, optional
        Integration step size for the numerically integrated laws. When left
        as `None` a value is chosen with `choose_step_size`, which costs an
        extra search. Pass a value directly inside a fitting loop.
    nthreads : int, optional
        Number of OpenMP threads. Only useful for long arrays.

    Returns
    -------
    numpy.ndarray
        Flux relative to the unocculted disk, same shape as ``z``.

    Notes
    -----
    The returned flux is normalised so that an unocculted disk gives exactly
    one. The central depth is deeper than the geometric area ratio, because the
    occulter covers the bright centre of the disk rather than an average patch
    of it. That excess is what makes an occultation light curve a measurement
    of limb darkening.

    Examples
    --------
    Venus at the centre of the disk, with rough photospheric limb darkening
    coefficients:

    >>> import numpy as np
    >>> flux = occultation_flux(np.array([0.0]), 0.0311, "quadratic", (0.42, 0.24))
    >>> float(np.round((1 - flux[0]) * 1e6))
    1179.0

    The geometric area ratio alone would give 967 parts per million, so limb
    darkening deepens the transit by about a fifth. That excess is the part of
    the signal that measures the limb darkening rather than the size.
    """
    z, ratio, coeffs = _check_inputs(z, ratio, law, coeffs)

    if ratio == 0.0:
        return np.ones_like(z)

    if law in _NUMERICAL_LAWS and step is None:
        step = choose_step_size(ratio, law, coeffs, nthreads=nthreads)
    if step is None:
        step = _STEP_FLOOR

    nthreads = int(nthreads)

    if law == "uniform":
        return _uniform_ld._uniform_ld(z, ratio, nthreads)
    if law == "linear":
        return _quadratic_ld._quadratic_ld(z, ratio, coeffs[0], 0.0, nthreads)
    if law == "quadratic":
        return _quadratic_ld._quadratic_ld(z, ratio, coeffs[0], coeffs[1], nthreads)
    if law == "nonlinear":
        return _nonlinear_ld._nonlinear_ld(
            z, ratio, coeffs[0], coeffs[1], coeffs[2], coeffs[3], step, nthreads
        )
    if law == "squareroot":
        # batman maps the square root law onto the nonlinear kernel with the
        # two coefficients swapped and the other two set to zero.
        return _nonlinear_ld._nonlinear_ld(
            z, ratio, coeffs[1], coeffs[0], 0.0, 0.0, step, nthreads
        )
    if law == "logarithmic":
        return _logarithmic_ld._logarithmic_ld(
            z, ratio, coeffs[0], coeffs[1], step, nthreads
        )
    if law == "exponential":
        return _exponential_ld._exponential_ld(
            z, ratio, coeffs[0], coeffs[1], step, nthreads
        )
    if law == "power2":
        return _power2_ld._power2_ld(z, ratio, coeffs[0], coeffs[1], step, nthreads)

    raise ValueError(f"unhandled limb darkening law {law!r}")
