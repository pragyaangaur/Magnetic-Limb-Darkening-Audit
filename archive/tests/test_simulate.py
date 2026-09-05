"""Tests for the synthetic disk renderer and the direct integration reference."""

import numpy as np
import pytest

from solarbatman import add_spot, brute_force_flux, disk_image, occulted_disk_image

VENUS = 0.0311
MOON = 0.97  # angular radius ratio from an Earth orbit, near an annular eclipse
SOLAR_U = (0.42, 0.24)


def test_disk_image_shape_and_radius():
    image, radius = disk_image(npix=200, radius_fraction=0.9)
    assert image.shape == (200, 200)
    assert radius == pytest.approx(90.0)


def test_disk_image_is_dark_outside_the_disk():
    image, radius = disk_image(npix=200, radius_fraction=0.8)
    corner = image[:10, :10]
    assert np.all(corner == 0.0)


def test_disk_image_is_brightest_at_the_centre():
    image, _ = disk_image(npix=201, radius_fraction=0.9)
    centre = image[100, 100]
    assert centre == pytest.approx(image.max())


def test_uniform_disk_has_no_radial_gradient():
    image, radius = disk_image(npix=200, law="uniform", coeffs=())
    on = image > 0
    assert np.allclose(image[on], 1.0)


def test_occulter_removes_flux():
    image, radius = disk_image(npix=400, radius_fraction=0.9)
    blocked = occulted_disk_image(image, radius, z=0.0, ratio=MOON)
    assert blocked.sum() < image.sum()
    lost = 1 - blocked.sum() / image.sum()
    # The Moon covers almost the whole disk, so almost all the light goes.
    assert lost > 0.95


def test_occulter_outside_the_disk_removes_nothing():
    image, radius = disk_image(npix=300, radius_fraction=0.5)
    blocked = occulted_disk_image(image, radius, z=1.0 + MOON + 0.05, ratio=MOON)
    np.testing.assert_array_equal(blocked, image)


def test_occulter_position_angle_does_not_change_the_flux():
    image, radius = disk_image(npix=600, radius_fraction=0.6)
    fluxes = [
        occulted_disk_image(image, radius, z=0.5, ratio=MOON, angle=a).sum()
        for a in (0.0, 0.7, 2.1, 4.9)
    ]
    assert np.ptp(fluxes) / fluxes[0] < 2e-3


def test_brute_force_matches_the_uniform_area_ratio():
    # Venus fits well inside the disk, so the central loss for a uniform disk
    # has to be exactly the area ratio.
    flux = brute_force_flux([0.0], VENUS, "uniform", (), npix=3000)[0]
    # Venus covers a disk 47 pixels across on this grid, so the staircase edge
    # costs about a part in a thousand of the depth.
    assert 1 - flux == pytest.approx(VENUS**2, rel=3e-3)


def test_brute_force_converges_with_grid_size():
    coarse = brute_force_flux([0.3], MOON, "quadratic", SOLAR_U, npix=500)[0]
    mid = brute_force_flux([0.3], MOON, "quadratic", SOLAR_U, npix=1500)[0]
    fine = brute_force_flux([0.3], MOON, "quadratic", SOLAR_U, npix=4000)[0]
    # Each refinement has to move the answer less than the one before it.
    assert abs(mid - fine) < abs(coarse - fine)
    assert abs(mid - fine) < 1e-4


def test_brute_force_is_one_outside_contact():
    flux = brute_force_flux([1.0 + MOON + 0.01], MOON, "quadratic", SOLAR_U, npix=800)
    assert flux[0] == pytest.approx(1.0, abs=1e-12)


def test_spot_darkens_only_inside_its_radius():
    image, radius = disk_image(npix=400, radius_fraction=0.9)
    spotted = add_spot(image, radius, x_frac=0.3, y_frac=0.0, spot_radius_frac=0.1)
    assert spotted.sum() < image.sum()
    # Far from the spot nothing changes.
    assert spotted[10, 10] == image[10, 10]


def test_spot_contrast_controls_how_much_is_lost():
    image, radius = disk_image(npix=400, radius_fraction=0.9)
    dark = add_spot(image, radius, 0.0, 0.0, 0.1, contrast=0.1)
    light = add_spot(image, radius, 0.0, 0.0, 0.1, contrast=0.8)
    assert dark.sum() < light.sum() < image.sum()


def test_disk_image_rejects_a_tiny_grid():
    with pytest.raises(ValueError, match="at least 8"):
        disk_image(npix=4)
