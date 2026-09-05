"""Synthetic limb darkened disks, and a slow reference occultation model.

Two jobs are done here. The first is to make images of a limb darkened disk
with an opaque circle in front of it, which is useful for testing the whole
pipeline without downloading anything. The second is to integrate those images
directly, which gives an occultation flux that depends on no analytic result at
all. That direct integration is the reference the batman kernels are checked
against.

Nothing in this module needs sunpy.
"""

from __future__ import annotations

import numpy as np

from .limbdark import intensity

__all__ = [
    "disk_image",
    "occulted_disk_image",
    "brute_force_flux",
    "add_spot",
]


def disk_image(npix=512, law="quadratic", coeffs=(0.42, 0.24), radius_fraction=0.9):
    """
    Render a limb darkened disk on a square grid.

    Parameters
    ----------
    npix : int, optional
        Side length of the output array in pixels.
    law : str, optional
        Limb darkening law.
    coeffs : array_like, optional
        Coefficients for the law.
    radius_fraction : float, optional
        Disk radius as a fraction of half the image side. A value below one
        leaves off disk margin, the way a real full disk image does.

    Returns
    -------
    image : numpy.ndarray
        Intensity, zero outside the disk.
    radius_pixels : float
        Disk radius in pixels, so callers can convert to solar radii.

    Notes
    -----
    Pixel centres are used, so the disk edge is a hard cut at the pixel centre
    rather than an area weighted edge. That makes the total flux converge as
    ``1 / npix`` rather than faster, which matters when this is used as an
    accuracy reference. Use a large ``npix`` for that purpose.
    """
    if npix < 8:
        raise ValueError("npix must be at least 8")
    radius = 0.5 * npix * float(radius_fraction)

    axis = np.arange(npix, dtype=np.float64) - 0.5 * (npix - 1)
    x, y = np.meshgrid(axis, axis, indexing="xy")
    r = np.hypot(x, y) / radius

    on_disk = r < 1.0
    mu = np.zeros_like(r)
    mu[on_disk] = np.sqrt(1.0 - r[on_disk] ** 2)

    image = np.zeros_like(r)
    image[on_disk] = intensity(mu[on_disk], law, coeffs)
    return image, radius


def occulted_disk_image(image, radius_pixels, z, ratio, angle=0.0):
    """
    Put an opaque circle in front of a rendered disk.

    Parameters
    ----------
    image : numpy.ndarray
        Disk image from `disk_image`.
    radius_pixels : float
        Disk radius in pixels.
    z : float
        Projected separation of the two centres, in disk radii.
    ratio : float
        Occulter radius over disk radius.
    angle : float, optional
        Position angle of the occulter in radians, measured from the positive
        x axis. The flux does not depend on this, so it only matters when the
        image itself is wanted.

    Returns
    -------
    numpy.ndarray
        The image with the occulted pixels set to zero.
    """
    npix = image.shape[0]
    axis = np.arange(npix, dtype=np.float64) - 0.5 * (npix - 1)
    x, y = np.meshgrid(axis, axis, indexing="xy")

    cx = z * radius_pixels * np.cos(angle)
    cy = z * radius_pixels * np.sin(angle)
    blocked = np.hypot(x - cx, y - cy) < ratio * radius_pixels

    out = image.copy()
    out[blocked] = 0.0
    return out


def brute_force_flux(z, ratio, law="quadratic", coeffs=(0.42, 0.24), npix=4000):
    """
    Occultation flux from direct summation over a rendered image.

    This is the reference the batman kernels are tested against. It is slow and
    its accuracy is set by ``npix``, but it makes no analytic assumption beyond
    the limb darkening law itself.

    Parameters
    ----------
    z : array_like
        Projected separations in disk radii.
    ratio : float
        Occulter radius over disk radius.
    law : str, optional
        Limb darkening law.
    coeffs : array_like, optional
        Coefficients for the law.
    npix : int, optional
        Grid side length. The error scales roughly as ``1 / npix``, so 4000
        gives a few parts in ten million and 1000 gives a few parts in a
        million.

    Returns
    -------
    numpy.ndarray
        Flux relative to the unocculted disk.
    """
    z = np.atleast_1d(np.asarray(z, dtype=np.float64))

    # A radius of exactly half the image keeps the whole disk in frame while
    # wasting as few pixels as possible.
    axis = (np.arange(npix, dtype=np.float64) + 0.5) / npix * 2.0 - 1.0
    x, y = np.meshgrid(axis, axis, indexing="xy")
    r2 = x**2 + y**2

    on_disk = r2 < 1.0
    mu = np.zeros_like(r2)
    mu[on_disk] = np.sqrt(1.0 - r2[on_disk])
    image = np.zeros_like(r2)
    image[on_disk] = intensity(mu[on_disk], law, coeffs)

    total = image.sum()
    if total <= 0:
        raise ValueError("rendered disk has non-positive total flux")

    out = np.empty(z.shape, dtype=np.float64)
    for i, zi in enumerate(z):
        blocked = (x - zi) ** 2 + y**2 < ratio * ratio
        out[i] = (image * ~blocked).sum() / total
    return out


def add_spot(image, radius_pixels, x_frac, y_frac, spot_radius_frac, contrast=0.3):
    """
    Darken a circular patch of a rendered disk, standing in for a sunspot.

    Parameters
    ----------
    image : numpy.ndarray
        Disk image from `disk_image`.
    radius_pixels : float
        Disk radius in pixels.
    x_frac, y_frac : float
        Spot centre in disk radii, measured from disk centre.
    spot_radius_frac : float
        Spot radius in disk radii.
    contrast : float, optional
        Intensity inside the spot as a fraction of the unspotted intensity at
        the same place. A large sunspot umbra sits near 0.2 in the continuum.

    Returns
    -------
    numpy.ndarray
        A copy of the image with the patch darkened.

    Notes
    -----
    The spot is placed in projection, so it is a circle on the image rather
    than a circle on the sphere. That is wrong near the limb and fine near disk
    centre. It is good enough for testing how a spot crossing biases an
    occultation fit, which is what this is for.
    """
    npix = image.shape[0]
    axis = np.arange(npix, dtype=np.float64) - 0.5 * (npix - 1)
    x, y = np.meshgrid(axis, axis, indexing="xy")

    inside = (
        np.hypot(x - x_frac * radius_pixels, y - y_frac * radius_pixels)
        < spot_radius_frac * radius_pixels
    )
    out = image.copy()
    out[inside] *= float(contrast)
    return out
