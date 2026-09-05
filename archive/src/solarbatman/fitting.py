"""Fitting an occultation light curve.

The model is a projected separation, a radius ratio, a limb darkening law and a
flat baseline. The separation usually comes from an ephemeris and is held
fixed, which leaves the radius ratio, the limb darkening coefficients and the
baseline to be fitted. That is a very different situation from an exoplanet
fit, where the geometry is the unknown and the limb darkening is often fixed to
a theoretical value.

Nothing here needs sunpy.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from .kernel import choose_step_size, n_coefficients, occultation_flux

__all__ = ["fit_occultation", "occultation_depth"]


def occultation_depth(ratio, law="quadratic", coeffs=(0.42, 0.24)):
    """
    Depth of an occultation at closest approach to the disk centre.

    Parameters
    ----------
    ratio : float
        Ratio of angular radii.
    law : str, optional
        Limb darkening law.
    coeffs : array_like, optional
        Coefficients for the law.

    Returns
    -------
    float
        Fraction of the light removed with the occulter at disk centre.

    Examples
    --------
    >>> round(float(occultation_depth(0.0311) * 1e6))
    1179
    """
    return float(1.0 - occultation_flux(np.array([0.0]), ratio, law, coeffs)[0])


def fit_occultation(
    z,
    flux,
    ratio_guess,
    law="quadratic",
    coeffs_guess=None,
    fit_ratio=True,
    fit_coeffs=True,
    fit_baseline=True,
    uncertainty=None,
):
    """
    Fit a radius ratio and limb darkening to a measured occultation.

    Parameters
    ----------
    z : array_like
        Projected separation at each sample, in solar radii. This is normally
        fixed by the ephemeris and is not fitted.
    flux : array_like
        Measured relative flux, close to one outside the occultation.
    ratio_guess : float
        Starting value for the ratio of angular radii.
    law : str, optional
        Limb darkening law to fit.
    coeffs_guess : array_like, optional
        Starting limb darkening coefficients. Defaults to rough photospheric
        values for the two coefficient laws and to a flat guess otherwise.
    fit_ratio : bool, optional
        Whether the radius ratio is free. Turn it off to fit limb darkening
        alone with the size fixed by the ephemeris, which is the strongest
        test the data can give.
    fit_coeffs : bool, optional
        Whether the limb darkening coefficients are free.
    fit_baseline : bool, optional
        Whether a multiplicative baseline is fitted alongside the model. Leave
        this on unless the light curve is already normalised to a known level.
    uncertainty : array_like, optional
        Per sample flux uncertainty. When given, the residuals are divided by
        it and the returned parameter errors mean something.

    Returns
    -------
    dict
        With keys ``ratio``, ``coeffs``, ``baseline``, ``residual_rms``,
        ``model``, ``success`` and ``errors``. ``errors`` is `None` when no
        uncertainty was supplied or when the covariance could not be formed.

    Notes
    -----
    Radius ratio and limb darkening are correlated, because both control how
    much light is lost. In exoplanet work that correlation is a serious problem
    and it is usually broken by fixing the limb darkening to a model. A solar
    occultation can break it the other way round. The ephemeris gives the
    radius ratio to far better than the photometry ever will, so setting
    ``fit_ratio=False`` leaves the limb darkening as the only free shape and
    turns the light curve into a clean measurement of it.
    """
    z = np.asarray(z, dtype=np.float64).ravel()
    flux = np.asarray(flux, dtype=np.float64).ravel()
    if z.shape != flux.shape:
        raise ValueError("z and flux must have the same shape")
    if z.size == 0:
        raise ValueError("no samples were supplied")

    n_coeff = n_coefficients(law)
    if coeffs_guess is None:
        coeffs_guess = {1: [0.5], 2: [0.42, 0.24], 4: [0.5, 0.1, 0.1, -0.1]}.get(
            n_coeff, [0.3] * n_coeff
        )
    coeffs_guess = np.atleast_1d(np.asarray(coeffs_guess, dtype=np.float64)).ravel()
    if coeffs_guess.size != n_coeff:
        raise ValueError(f"law {law!r} takes {n_coeff} coefficients")

    if uncertainty is None:
        weights = np.ones_like(flux)
    else:
        weights = np.asarray(uncertainty, dtype=np.float64).ravel()
        if weights.shape != flux.shape:
            raise ValueError("uncertainty must have the same shape as flux")
        if np.any(weights <= 0):
            raise ValueError("uncertainty must be positive everywhere")

    # The step size for the numerically integrated laws is held fixed through
    # the fit. Re-searching it at every step would make the residual surface
    # jump around and stop the optimiser from converging.
    step = choose_step_size(float(ratio_guess), law, coeffs_guess)

    # Build the free parameter vector and remember how to take it apart again.
    free = []
    if fit_ratio:
        free.append(("ratio", float(ratio_guess)))
    if fit_coeffs:
        free.extend((f"c{i}", float(v)) for i, v in enumerate(coeffs_guess))
    if fit_baseline:
        free.append(("baseline", float(np.median(flux))))
    if not free:
        raise ValueError("nothing is free to fit")

    names = [name for name, _ in free]
    start = np.array([value for _, value in free], dtype=np.float64)

    def unpack(p):
        values = dict(zip(names, p, strict=True))
        ratio = values.get("ratio", float(ratio_guess))
        if fit_coeffs:
            coeffs = np.array([values[f"c{i}"] for i in range(n_coeff)])
        else:
            coeffs = coeffs_guess
        baseline = values.get("baseline", 1.0)
        return abs(ratio), coeffs, baseline

    def model_for(p):
        ratio, coeffs, baseline = unpack(p)
        return baseline * occultation_flux(z, ratio, law, coeffs, step=step)

    def residual(p):
        return (model_for(p) - flux) / weights

    result = least_squares(residual, start, method="trf", max_nfev=20000)
    ratio, coeffs, baseline = unpack(result.x)
    model = model_for(result.x)

    errors = None
    if uncertainty is not None and result.jac.size:
        try:
            covariance = np.linalg.inv(result.jac.T @ result.jac)
            errors = dict(
                zip(names, np.sqrt(np.diag(covariance)), strict=True)
            )
        except np.linalg.LinAlgError:
            errors = None

    return {
        "ratio": float(ratio),
        "coeffs": coeffs,
        "baseline": float(baseline),
        "law": law,
        "model": model,
        "residual_rms": float(np.sqrt(np.mean((model - flux) ** 2))),
        "success": bool(result.success),
        "errors": errors,
        "n_free": len(names),
        "free_parameters": names,
    }
