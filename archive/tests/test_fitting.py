"""Tests for the occultation fitter.

Every test here builds a light curve from the model itself and then tries to
get the input parameters back out. That checks the fitter and the plumbing
around it. It says nothing about whether the model describes the real Sun,
which is what a fit to real data is for.
"""

import numpy as np
import pytest

from solarbatman import occultation_flux
from solarbatman.fitting import fit_occultation, occultation_depth
from solarbatman.geometry import linear_track

MERCURY = 0.0058
VENUS = 0.0311
MOON = 0.97
SOLAR_U = (0.42, 0.24)


def _light_curve(ratio, coeffs=SOLAR_U, law="quadratic", impact=0.3, n=600, noise=0.0, seed=1):
    t = np.linspace(-1.6, 1.6, n)
    z = linear_track(t, impact_parameter=impact)
    flux = occultation_flux(z, ratio, law, coeffs)
    if noise:
        flux = flux + np.random.default_rng(seed).normal(0, noise, flux.size)
    return z, flux


def test_depth_matches_the_kernel():
    assert occultation_depth(VENUS) == pytest.approx(
        1 - occultation_flux(np.array([0.0]), VENUS, "quadratic", SOLAR_U)[0]
    )


def test_depth_grows_with_the_ratio():
    assert occultation_depth(MERCURY) < occultation_depth(VENUS) < occultation_depth(MOON)


@pytest.mark.parametrize("ratio", [VENUS, 0.5, MOON])
def test_recovers_the_ratio_from_a_clean_curve(ratio):
    z, flux = _light_curve(ratio)
    result = fit_occultation(z, flux, ratio_guess=ratio * 1.3, fit_coeffs=False)
    assert result["success"]
    assert result["ratio"] == pytest.approx(ratio, rel=1e-4)


def test_recovers_ratio_and_limb_darkening_together():
    z, flux = _light_curve(MOON)
    result = fit_occultation(
        z, flux, ratio_guess=0.9, coeffs_guess=(0.3, 0.3)
    )
    assert result["ratio"] == pytest.approx(MOON, rel=1e-3)
    np.testing.assert_allclose(result["coeffs"], SOLAR_U, atol=5e-3)
    assert result["baseline"] == pytest.approx(1.0, rel=1e-4)


def test_fixing_the_ratio_still_recovers_limb_darkening():
    # This is the solar case. The ephemeris knows the size far better than the
    # photometry does, so the size is fixed and the limb darkening is the only
    # free shape left.
    z, flux = _light_curve(MOON)
    result = fit_occultation(
        z, flux, ratio_guess=MOON, fit_ratio=False, coeffs_guess=(0.2, 0.4)
    )
    assert result["ratio"] == pytest.approx(MOON)
    np.testing.assert_allclose(result["coeffs"], SOLAR_U, atol=1e-4)


def test_a_free_baseline_absorbs_a_scale_error():
    z, flux = _light_curve(VENUS)
    result = fit_occultation(z, 1.05 * flux, ratio_guess=VENUS)
    assert result["baseline"] == pytest.approx(1.05, rel=1e-4)
    assert result["ratio"] == pytest.approx(VENUS, rel=1e-2)


def test_noise_widens_the_answer_without_biasing_it():
    z, flux = _light_curve(MOON, noise=1e-3, seed=7)
    result = fit_occultation(z, flux, ratio_guess=0.9)
    assert result["ratio"] == pytest.approx(MOON, abs=5e-3)
    assert result["residual_rms"] == pytest.approx(1e-3, rel=0.3)


def test_uncertainties_come_back_when_they_are_supplied():
    sigma = 1e-3
    z, flux = _light_curve(MOON, noise=sigma, seed=11)
    result = fit_occultation(
        z, flux, ratio_guess=0.9, uncertainty=np.full(flux.size, sigma)
    )
    assert result["errors"] is not None
    assert result["errors"]["ratio"] > 0
    # The recovered value should sit within a few error bars of the truth.
    assert abs(result["ratio"] - MOON) < 5 * result["errors"]["ratio"]


@pytest.mark.parametrize("law", ["linear", "quadratic", "power2", "nonlinear"])
def test_works_for_several_laws(law):
    truth = {
        "linear": (0.5,),
        "quadratic": (0.42, 0.24),
        "power2": (0.5, 0.6),
        "nonlinear": (0.5, 0.1, 0.1, -0.1),
    }[law]
    z, flux = _light_curve(MOON, coeffs=truth, law=law)
    result = fit_occultation(
        z, flux, ratio_guess=MOON, law=law, fit_ratio=False, fit_coeffs=True
    )
    np.testing.assert_allclose(result["coeffs"], truth, atol=2e-3)


def test_mismatched_shapes_are_rejected():
    with pytest.raises(ValueError, match="same shape"):
        fit_occultation(np.zeros(10), np.zeros(11), ratio_guess=0.1)


def test_wrong_coefficient_count_is_rejected():
    z, flux = _light_curve(VENUS)
    with pytest.raises(ValueError, match="takes 2 coefficients"):
        fit_occultation(z, flux, ratio_guess=VENUS, coeffs_guess=(0.3,))


def test_nothing_free_is_rejected():
    z, flux = _light_curve(VENUS)
    with pytest.raises(ValueError, match="nothing is free"):
        fit_occultation(
            z,
            flux,
            ratio_guess=VENUS,
            fit_ratio=False,
            fit_coeffs=False,
            fit_baseline=False,
        )


def test_negative_uncertainty_is_rejected():
    z, flux = _light_curve(VENUS)
    with pytest.raises(ValueError, match="must be positive"):
        fit_occultation(z, flux, ratio_guess=VENUS, uncertainty=np.zeros(flux.size))
