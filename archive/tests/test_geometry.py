"""Tests for the occultation geometry.

The ephemeris functions need sunpy and a network call to JPL Horizons, so they
are not covered here. Everything that is pure arithmetic is, which is most of
what can go wrong.
"""

import numpy as np
import pytest

from solarbatman.geometry import (
    BODY_RADII_KM,
    SOLAR_RADIUS_KM,
    angular_radius_ratio,
    contact_separations,
    impact_parameter_from_track,
    linear_track,
)

# Rough distances in kilometres, good enough to check the arithmetic.
EARTH_MOON = 3.844e5
EARTH_SUN = 1.496e8


def test_moon_ratio_from_earth_is_close_to_one():
    # The Moon and the Sun look almost the same size from an Earth orbit, which
    # is why solar eclipses are sometimes total and sometimes annular. The
    # physical radius ratio is 0.0025 and is not the number an occultation uses.
    ratio = angular_radius_ratio("moon", EARTH_MOON, EARTH_SUN)
    assert ratio == pytest.approx(0.972, abs=0.01)


def test_moon_ratio_crosses_one_over_the_lunar_orbit():
    # Perigee near 363000 km gives a total eclipse, apogee near 405000 km an
    # annular one.
    perigee = angular_radius_ratio("moon", 3.63e5, EARTH_SUN)
    apogee = angular_radius_ratio("moon", 4.05e5, EARTH_SUN)
    assert apogee < 1.0 < perigee


def test_ratio_accepts_a_radius_in_kilometres():
    by_name = angular_radius_ratio("moon", EARTH_MOON, EARTH_SUN)
    by_value = angular_radius_ratio(BODY_RADII_KM["moon"], EARTH_MOON, EARTH_SUN)
    assert by_name == pytest.approx(by_value)


def test_ratio_scales_the_right_way_with_distance():
    near = angular_radius_ratio("moon", EARTH_MOON * 0.9, EARTH_SUN)
    far = angular_radius_ratio("moon", EARTH_MOON * 1.1, EARTH_SUN)
    assert near > far


def test_mercury_and_venus_ratios_are_in_the_expected_range():
    # Mercury at inferior conjunction is about 0.6 au from the Sun side of
    # Earth, and Venus about 0.28 au.
    mercury = angular_radius_ratio("mercury", 0.60 * EARTH_SUN, EARTH_SUN)
    venus = angular_radius_ratio("venus", 0.28 * EARTH_SUN, EARTH_SUN)
    assert mercury == pytest.approx(0.0058, abs=1e-3)
    assert venus == pytest.approx(0.0311, abs=2e-3)


def test_unknown_body_is_rejected():
    with pytest.raises(ValueError, match="unknown body"):
        angular_radius_ratio("pluto", EARTH_MOON, EARTH_SUN)


def test_solar_radius_matches_the_iau_nominal_value():
    assert SOLAR_RADIUS_KM == pytest.approx(695700.0)


def test_linear_track_minimum_is_the_impact_parameter():
    t = np.linspace(-5, 5, 1001)
    z = linear_track(t, t0=0.0, impact_parameter=0.42, crossing_time=1.0)
    assert z.min() == pytest.approx(0.42)


def test_linear_track_is_symmetric_about_mid_transit():
    t = np.linspace(-3, 3, 601)
    z = linear_track(t, t0=0.0, impact_parameter=0.2)
    np.testing.assert_allclose(z, z[::-1], atol=1e-12)


def test_linear_track_offset_shifts_the_minimum():
    t = np.linspace(0, 10, 1001)
    z = linear_track(t, t0=3.5, impact_parameter=0.1)
    assert t[np.argmin(z)] == pytest.approx(3.5, abs=0.01)


def test_linear_track_crossing_time_stretches_the_event():
    t = np.linspace(-2, 2, 401)
    fast = linear_track(t, crossing_time=0.5)
    slow = linear_track(t, crossing_time=2.0)
    assert fast.max() > slow.max()


def test_linear_track_rejects_a_non_positive_crossing_time():
    with pytest.raises(ValueError, match="must be positive"):
        linear_track(np.array([0.0]), crossing_time=0.0)


def test_contacts_bracket_the_occultation():
    contacts = contact_separations(0.0311)
    assert contacts["first"] == pytest.approx(1.0311)
    assert contacts["second"] == pytest.approx(0.9689)
    assert contacts["third"] == contacts["second"]
    assert contacts["fourth"] == contacts["first"]


def test_inner_contacts_vanish_for_an_oversized_occulter():
    contacts = contact_separations(1.5)
    assert contacts["second"] is None
    assert contacts["first"] == pytest.approx(2.5)


def test_contacts_reject_a_negative_ratio():
    with pytest.raises(ValueError, match="non-negative"):
        contact_separations(-0.1)


def test_impact_parameter_reads_the_sampled_minimum():
    z = np.array([1.4, 0.9, 0.31, 0.55, 1.2])
    assert impact_parameter_from_track(z) == pytest.approx(0.31)


def test_impact_parameter_rejects_an_empty_track():
    with pytest.raises(ValueError, match="empty"):
        impact_parameter_from_track(np.array([]))


def test_ephemeris_helpers_explain_the_missing_dependency():
    from solarbatman import geometry

    try:
        import sunpy.coordinates  # noqa: F401
    except ImportError:
        with pytest.raises(ImportError, match="need sunpy"):
            geometry.separation_from_ephemeris(["2019-11-11"], "Mercury")
    else:
        pytest.skip("sunpy is importable, so the guard cannot be exercised")
