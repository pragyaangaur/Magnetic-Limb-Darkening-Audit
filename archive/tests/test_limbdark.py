"""Tests for the limb darkening laws and the profile fitter."""

import numpy as np
import pytest

from solarbatman.kernel import LIMB_DARKENING_LAWS
from solarbatman.limbdark import SOLAR_REFERENCE, fit_intensity_profile, intensity

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

# The exponential law does not give one at disk centre, because its second term
# stays finite there. Every other law does.
_NORMALISED_AT_CENTRE = [law for law in COEFFS if law != "exponential"]


@pytest.mark.parametrize("law", sorted(_NORMALISED_AT_CENTRE))
def test_intensity_is_one_at_disk_centre(law):
    assert intensity(1.0, law, COEFFS[law]) == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("law", sorted(LIMB_DARKENING_LAWS))
def test_intensity_is_finite_across_the_disk(law):
    mu = np.linspace(0.0, 1.0, 200)
    values = intensity(mu, law, COEFFS[law])
    assert np.all(np.isfinite(values))


@pytest.mark.parametrize("law", sorted(_NORMALISED_AT_CENTRE))
def test_intensity_falls_towards_the_limb(law):
    if law == "uniform":
        pytest.skip("a uniform disk has no limb darkening by definition")
    mu = np.linspace(0.05, 1.0, 200)
    values = intensity(mu, law, COEFFS[law])
    assert values[0] < values[-1]


def test_solar_reference_darkens_less_at_longer_wavelength():
    # Limb darkening weakens towards the infrared, because the source function
    # varies less over the height range the continuum is formed in.
    limb = {w: intensity(0.1, "quadratic", u) for w, u in SOLAR_REFERENCE.items()}
    wavelengths = sorted(limb)
    values = [limb[w] for w in wavelengths]
    assert values == sorted(values)


@pytest.mark.parametrize("law", ["linear", "quadratic", "squareroot", "power2", "nonlinear"])
def test_fit_recovers_its_own_coefficients(law):
    truth = np.asarray(COEFFS[law], dtype=float)
    mu = np.linspace(0.05, 1.0, 500)
    profile = 1234.5 * intensity(mu, law, truth)

    result = fit_intensity_profile(mu, profile, law=law)

    np.testing.assert_allclose(result["coeffs"], truth, atol=1e-6)
    assert result["scale"] == pytest.approx(1234.5, rel=1e-6)
    assert result["residual_rms"] < 1e-8


def test_fit_survives_noise():
    rng = np.random.default_rng(20260903)
    truth = (0.42, 0.24)
    mu = np.linspace(0.05, 1.0, 4000)
    profile = intensity(mu, "quadratic", truth) + rng.normal(0, 2e-3, mu.size)

    result = fit_intensity_profile(mu, profile, law="quadratic")

    np.testing.assert_allclose(result["coeffs"], truth, atol=0.01)


def test_fit_drops_the_extreme_limb():
    mu = np.linspace(0.0, 1.0, 100)
    profile = intensity(mu, "quadratic", (0.42, 0.24))
    result = fit_intensity_profile(mu, profile, mu_min=0.3)
    assert result["n_points"] == int(np.sum(mu >= 0.3))


def test_fit_rejects_mismatched_shapes():
    with pytest.raises(ValueError, match="same shape"):
        fit_intensity_profile(np.linspace(0, 1, 10), np.linspace(0, 1, 11))


def test_fit_rejects_an_unknown_law():
    with pytest.raises(ValueError, match="unknown law"):
        fit_intensity_profile(np.linspace(0.1, 1, 10), np.ones(10), law="parabolic")
