"""Tests for the occultation kernels.

The important test here is `test_matches_brute_force_at_solar_ratios`. Exoplanet
work lives below a radius ratio of about 0.2 and that is where batman is tested.
The Moon seen from an Earth orbit has an angular radius ratio near one, which is
far outside that range and includes the total eclipse case, so the kernel is
checked against a direct integration across the whole span.
"""

import numpy as np
import pytest

from solarbatman import brute_force_flux, occultation_flux
from solarbatman.kernel import LIMB_DARKENING_LAWS

# Ratios of angular radii for the three bodies that really do cross the solar
# disk as seen from an Earth orbit. These are not physical radius ratios. The
# lunar value moves between about 0.97 and 1.03 through the year, which is the
# difference between an annular and a total eclipse.
MERCURY = 0.0058
VENUS = 0.0311
MOON_ANNULAR = 0.97
MOON_TOTAL = 1.03
MOON = MOON_ANNULAR

# Illustrative photospheric coefficients near the HMI continuum wavelength.
SOLAR_U = (0.42, 0.24)

COEFFS = {
    "uniform": (),
    "linear": (0.5,),
    "quadratic": (0.42, 0.24),
    "squareroot": (0.3, 0.3),
    "logarithmic": (0.5, 0.2),
    "exponential": (0.5, 0.05),
    "power2": (0.5, 0.6),
    "nonlinear": (0.5, 0.1, 0.1, -0.1),
}


@pytest.mark.parametrize("law", sorted(LIMB_DARKENING_LAWS))
def test_every_law_runs_and_is_bounded(law):
    z = np.linspace(0.0, 1.5, 200)
    flux = occultation_flux(z, MOON, law, COEFFS[law])
    assert flux.shape == z.shape
    assert np.all(np.isfinite(flux))
    assert np.all(flux <= 1.0 + 1e-12)
    assert np.all(flux > 0.0)


@pytest.mark.parametrize("law", sorted(LIMB_DARKENING_LAWS))
def test_flux_is_one_outside_contact(law):
    z = np.array([1.0 + MOON + 1e-6, 2.0, 10.0])
    flux = occultation_flux(z, MOON, law, COEFFS[law])
    np.testing.assert_allclose(flux, 1.0, atol=1e-12)


# The exponential law has a term in 1 / (1 - exp(mu)), which is infinite at
# mu = 0. The intensity therefore blows up at the very limb instead of falling
# to a finite value, so the light curve is not monotonic while the occulter is
# crossing the limb. This is a property of the law, not a defect in batman, and
# it is the reason the exponential law should not be used for solar work.
@pytest.mark.parametrize(
    "law", sorted(set(LIMB_DARKENING_LAWS) - {"exponential"})
)
def test_flux_rises_monotonically_with_separation(law):
    z = np.linspace(0.0, 1.5, 400)
    flux = occultation_flux(z, MOON, law, COEFFS[law])
    assert np.all(np.diff(flux) >= -1e-9)


def test_exponential_law_is_not_monotonic_at_the_limb():
    # Recorded so that a future change in batman which fixes or worsens this is
    # noticed. The dip sits right where the occulter crosses the limb.
    z = np.linspace(0.0, 1.5, 400)
    flux = occultation_flux(z, VENUS, "exponential", COEFFS["exponential"])
    assert np.diff(flux).min() < -1e-6


def test_zero_ratio_gives_no_occultation():
    z = np.linspace(0.0, 1.0, 10)
    np.testing.assert_array_equal(occultation_flux(z, 0.0, "quadratic", SOLAR_U), 1.0)


def test_uniform_depth_is_the_geometric_area_ratio():
    # With no limb darkening the central depth has to be the area ratio.
    for ratio in (MERCURY, VENUS, MOON):
        depth = 1.0 - occultation_flux(np.array([0.0]), ratio, "uniform", ())[0]
        assert depth == pytest.approx(ratio**2, rel=1e-9)


def test_limb_darkening_deepens_the_eclipse():
    # The occulter covers the bright centre of the disk, so a limb darkened
    # disk always loses more than the geometric area ratio.
    plain = occultation_flux(np.array([0.0]), MOON, "uniform", ())[0]
    darkened = occultation_flux(np.array([0.0]), MOON, "quadratic", SOLAR_U)[0]
    assert darkened < plain


@pytest.mark.parametrize("ratio", [MERCURY, VENUS, 0.5, MOON_ANNULAR, 1.0, MOON_TOTAL])
def test_matches_brute_force_at_solar_ratios(ratio):
    z = np.array([0.0, 0.15, 0.35, 0.55, 0.75, 0.9, 1.0, 1.1])
    analytic = occultation_flux(z, ratio, "quadratic", SOLAR_U)
    numeric = brute_force_flux(z, ratio, "quadratic", SOLAR_U, npix=3000)
    # The direct integration converges as one over the grid side, and its error
    # grows with the length of the occulter edge, so a large ratio needs a
    # looser bound than a small one. Both are far below any realistic solar
    # photometric noise floor.
    tolerance = 5e-6 if ratio < 0.5 else 1e-5
    np.testing.assert_allclose(analytic, numeric, atol=tolerance)


def test_matches_brute_force_for_a_numerically_integrated_law():
    z = np.array([0.0, 0.4, 0.8, 1.05])
    analytic = occultation_flux(z, MOON, "nonlinear", COEFFS["nonlinear"])
    numeric = brute_force_flux(z, MOON, "nonlinear", COEFFS["nonlinear"], npix=3000)
    np.testing.assert_allclose(analytic, numeric, atol=1e-5)


def test_linear_is_quadratic_with_a_zero_second_coefficient():
    z = np.linspace(0.0, 1.2, 50)
    a = occultation_flux(z, MOON, "linear", (0.42,))
    b = occultation_flux(z, MOON, "quadratic", (0.42, 0.0))
    np.testing.assert_allclose(a, b, atol=1e-14)


def test_wrong_coefficient_count_is_rejected():
    with pytest.raises(ValueError, match="takes 2 coefficients"):
        occultation_flux(np.array([0.0]), MOON, "quadratic", (0.4,))


def test_unknown_law_is_rejected():
    with pytest.raises(ValueError, match="unknown limb darkening law"):
        occultation_flux(np.array([0.0]), MOON, "parabolic", (0.1, 0.2))


def test_negative_separation_is_rejected():
    with pytest.raises(ValueError, match="non-negative"):
        occultation_flux(np.array([-0.1]), MOON, "quadratic", SOLAR_U)


def test_step_size_search_returns_something_usable():
    from solarbatman.kernel import choose_step_size

    step = choose_step_size(MOON, "nonlinear", COEFFS["nonlinear"], max_error_ppm=1.0)
    assert 5.0e-4 <= step <= 1.0
    coarse = occultation_flux(
        np.linspace(0, 1.2, 100), MOON, "nonlinear", COEFFS["nonlinear"], step=step
    )
    fine = occultation_flux(
        np.linspace(0, 1.2, 100), MOON, "nonlinear", COEFFS["nonlinear"], step=5.0e-4
    )
    assert np.abs(coarse - fine).max() * 1e6 < 2.0
