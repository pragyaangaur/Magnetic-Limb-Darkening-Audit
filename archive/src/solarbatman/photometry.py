"""Sun as a star photometry from full disk images.

An occultation light curve is a plot of total solar brightness against time. A
full disk image already contains that number, as the sum over every pixel on
the disk. The work is in the normalisation rather than in the summing, because
the raw sum changes with exposure time, with the number of pixels the disk
covers and with the response of the instrument.

Everything here needs sunpy, which is imported when the functions are called.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "disk_integrated_intensity",
    "light_curve",
    "radial_profile",
]


def _require_sunpy():
    try:
        import astropy.units as u
        from sunpy.map.maputils import all_coordinates_from_map, coordinate_is_on_solar_disk
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError(
            "photometry needs sunpy. Install it with: pip install 'solarbatman[solar]'"
        ) from exc
    return u, all_coordinates_from_map, coordinate_is_on_solar_disk


def _on_disk_mask(smap):
    """
    Boolean mask of the pixels that fall on the solar disk.

    sunpy gained an ``on_disk_mask`` helper in ``sunpy.map.maputils``. It is
    used when it is there, and rebuilt from its two parts when it is not, so
    this package works against older sunpy releases as well.
    """
    _, all_coordinates_from_map, coordinate_is_on_solar_disk = _require_sunpy()
    try:
        from sunpy.map.maputils import on_disk_mask
    except ImportError:
        return coordinate_is_on_solar_disk(all_coordinates_from_map(smap))
    return on_disk_mask(smap)


def disk_integrated_intensity(smap, normalise=True, mask=None):
    """
    Total brightness of the solar disk in one image.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A full disk image. The map has to contain the whole disk, because a
        partial image measures a different quantity every time the pointing
        moves.
    normalise : bool, optional
        When true the sum is divided by the exposure time in seconds and by
        the number of on disk pixels, which gives a mean intensity per pixel
        per second. When false the raw sum is returned.
    mask : numpy.ndarray, optional
        A precomputed on disk mask. Building the mask means computing a
        coordinate for every pixel, which is by far the slowest part of this
        function, so pass the same mask for every image in a sequence when the
        pointing does not move.

    Returns
    -------
    float

    Notes
    -----
    Three things make the difference between a light curve that shows an
    occultation and one that shows only instrument drift.

    The exposure time has to divide out, because it is not always constant.
    Automatic exposure control changes it, and during a deep lunar eclipse it
    can change a great deal.

    The pixel count has to divide out as well. The Sun to observer distance
    changes by about three percent over the year, so the disk covers about six
    percent more pixels in January than in July. That is a far larger effect
    than any transit.

    Instrument degradation is not handled here at all. For AIA, run the images
    through the degradation correction in ``aiapy`` before calling this.
    """
    u, _, _ = _require_sunpy()

    if mask is None:
        mask = _on_disk_mask(smap)
    if not np.any(mask):
        raise ValueError("no pixels fall on the solar disk in this map")

    total = float(np.nansum(smap.data[mask]))
    if not normalise:
        return total

    n_pixels = int(np.count_nonzero(mask))
    exposure = smap.exposure_time
    if exposure is None:
        raise ValueError(
            "the map has no exposure time, so it cannot be normalised. Pass "
            "normalise=False and handle the scaling yourself."
        )
    seconds = float(exposure.to_value(u.s))
    if seconds <= 0:
        raise ValueError(f"exposure time is {seconds} s, which cannot be used")

    return total / (seconds * n_pixels)


def light_curve(maps, normalise=True, reuse_mask=False):
    """
    Build a Sun as a star light curve from a sequence of full disk images.

    Parameters
    ----------
    maps : iterable of `~sunpy.map.GenericMap`
        Images in time order.
    normalise : bool, optional
        Passed to `disk_integrated_intensity`.
    reuse_mask : bool, optional
        When true the on disk mask from the first image is used for all of
        them. This is much faster and is safe only when the pointing and the
        apparent solar radius do not change over the sequence, so use it for a
        few hours of data and not for a few months.

    Returns
    -------
    `~astropy.table.QTable`
        With columns ``time`` and ``flux``. The flux is normalised to its own
        median, so a quiet sequence sits at one.

    Notes
    -----
    Normalising to the median makes the output a relative light curve, which
    is what `solarbatman.occultation_flux` returns. It also throws away the
    absolute calibration, so do not use this to measure solar irradiance.
    """
    from astropy.table import QTable

    times, values, mask = [], [], None
    for smap in maps:
        if reuse_mask and mask is None:
            mask = _on_disk_mask(smap)
        values.append(
            disk_integrated_intensity(smap, normalise=normalise, mask=mask)
        )
        times.append(smap.date)

    if not values:
        raise ValueError("no maps were supplied")

    values = np.asarray(values, dtype=np.float64)
    median = np.median(values)
    if median <= 0:
        raise ValueError("the light curve has a non-positive median")

    return QTable({"time": times, "flux": values / median})


def radial_profile(smap, n_bins=200, mu_min=0.0):
    """
    Intensity against ``mu`` for a full disk image.

    This is the measurement that feeds
    `solarbatman.limbdark.fit_intensity_profile`, which turns it into limb
    darkening coefficients.

    Parameters
    ----------
    smap : `~sunpy.map.GenericMap`
        A full disk image. Use a continuum or a photospheric band. The extreme
        ultraviolet channels see optically thin coronal emission on top of the
        disk and do not follow a photospheric limb darkening law.
    n_bins : int, optional
        Number of equally spaced bins in ``mu``.
    mu_min : float, optional
        Bins below this ``mu`` are dropped.

    Returns
    -------
    mu : numpy.ndarray
        Bin centres.
    intensity : numpy.ndarray
        Median intensity in each bin. The median is used rather than the mean
        so that sunspots and bright points do not drag the profile around.
    counts : numpy.ndarray
        Number of pixels in each bin, useful for weighting a fit.

    Notes
    -----
    ``mu`` is computed from the projected distance from disk centre, as
    ``sqrt(1 - r**2)`` with ``r`` in units of the apparent solar radius. That
    is the correct relation for a sphere seen from far away, and the error
    from the finite observer distance is well below anything else in the
    measurement.
    """
    u, all_coordinates_from_map, _ = _require_sunpy()
    from sunpy.map.maputils import solar_angular_radius

    coords = all_coordinates_from_map(smap)
    radius = solar_angular_radius(coords).to_value(u.arcsec)

    offset = np.hypot(
        coords.Tx.to_value(u.arcsec), coords.Ty.to_value(u.arcsec)
    )
    r = offset / radius

    on_disk = (r < 1.0) & np.isfinite(smap.data)
    if not np.any(on_disk):
        raise ValueError("no usable on disk pixels in this map")

    mu_pixels = np.sqrt(1.0 - r[on_disk] ** 2)
    values = smap.data[on_disk]

    edges = np.linspace(0.0, 1.0, int(n_bins) + 1)
    index = np.clip(np.digitize(mu_pixels, edges) - 1, 0, n_bins - 1)

    centres = 0.5 * (edges[:-1] + edges[1:])
    profile = np.full(n_bins, np.nan)
    counts = np.zeros(n_bins, dtype=int)
    for b in range(n_bins):
        selected = values[index == b]
        counts[b] = selected.size
        if selected.size:
            profile[b] = np.median(selected)

    keep = (centres >= mu_min) & np.isfinite(profile)
    return centres[keep], profile[keep], counts[keep]
