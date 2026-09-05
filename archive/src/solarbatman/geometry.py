"""Turning an occultation geometry into a projected separation.

The batman kernels need one input that carries all the geometry: the sky
projected separation between the centre of the solar disk and the centre of the
occulter, in solar radii. Everything in this module produces that number.

There are two ways to get it. The first uses an ephemeris, which is what a real
transit needs, and it goes through sunpy. The second builds a straight track
across the disk, which is enough for testing and for understanding the shape of
a light curve, and it needs nothing but numpy.

The functions that need sunpy import it when they are called, so the rest of the
package works without it.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "BODY_RADII_KM",
    "SOLAR_RADIUS_KM",
    "angular_radius_ratio",
    "linear_track",
    "contact_separations",
    "impact_parameter_from_track",
    "separation_from_ephemeris",
    "radius_ratio_from_ephemeris",
]

# Mean radii in kilometres. The solar value is the IAU 2015 nominal solar
# radius, which is the one sunpy and astropy also use.
SOLAR_RADIUS_KM = 695700.0

BODY_RADII_KM = {
    "moon": 1737.4,
    "mercury": 2439.7,
    "venus": 6051.8,  # solid body, the cloud tops sit about 70 km higher
    "earth": 6371.0,
    "mars": 3389.5,
    "phobos": 11.1,
    "deimos": 6.2,
}


def angular_radius_ratio(body, body_distance_km, solar_distance_km):
    """
    Radius ratio of an occulter against the Sun, as seen by one observer.

    Parameters
    ----------
    body : str or float
        Name of a body in `BODY_RADII_KM`, or a radius in kilometres.
    body_distance_km : float or array_like
        Observer to occulter distance in kilometres.
    solar_distance_km : float or array_like
        Observer to Sun distance in kilometres.

    Returns
    -------
    float or numpy.ndarray
        Occulter angular radius over solar angular radius. This is what the
        kernels call ``ratio``.

    Notes
    -----
    An occultation depends on the ratio of angular sizes. The physical radius
    ratio is a different number and it is not the one the kernels want. The
    Moon is about 400 times smaller than the Sun and about 400 times closer, so
    its angular ratio comes out near one and it very nearly covers the whole
    solar disk. The physical ratio is 0.0025, which would give a depth of six
    parts per million instead of a total eclipse.

    SDO sits in a geosynchronous orbit, so its distances are close to the
    Earth centred values and the ratio runs from about 0.97 to about 1.03
    across the year. Below one the Moon fits inside the solar disk and the
    eclipse is annular. Above one it covers the disk completely.

    Examples
    --------
    The Moon from roughly the Earth to Moon distance:

    >>> round(float(angular_radius_ratio("moon", 3.844e5, 1.496e8)), 3)
    0.972
    """
    if isinstance(body, str):
        try:
            radius_km = BODY_RADII_KM[body.lower()]
        except KeyError:
            raise ValueError(
                f"unknown body {body!r}, expected one of {sorted(BODY_RADII_KM)} "
                "or a radius in kilometres"
            ) from None
    else:
        radius_km = float(body)

    body_angular = np.arcsin(radius_km / np.asarray(body_distance_km, dtype=float))
    solar_angular = np.arcsin(
        SOLAR_RADIUS_KM / np.asarray(solar_distance_km, dtype=float)
    )
    return body_angular / solar_angular


def linear_track(t, t0=0.0, impact_parameter=0.0, crossing_time=1.0):
    """
    Projected separation for an occulter moving in a straight line.

    Over the few hours an occultation lasts, the occulter moves across the disk
    on a track that is very close to straight. This is the toy geometry that
    gives the familiar transit shape, and it is enough for tests and for first
    guesses in a fit.

    Parameters
    ----------
    t : array_like
        Times, in any unit, as long as ``t0`` and ``crossing_time`` match.
    t0 : float, optional
        Time of closest approach to the disk centre.
    impact_parameter : float, optional
        Smallest separation reached, in solar radii. Zero is a central
        crossing and one grazes the limb.
    crossing_time : float, optional
        Time taken to move one solar radius along the track.

    Returns
    -------
    numpy.ndarray
        Projected separation in solar radii.

    Examples
    --------
    >>> import numpy as np
    >>> float(linear_track(np.array([0.0]), impact_parameter=0.3)[0])
    0.3
    """
    t = np.asarray(t, dtype=np.float64)
    if crossing_time <= 0:
        raise ValueError("crossing_time must be positive")
    along = (t - t0) / crossing_time
    return np.hypot(along, float(impact_parameter))


def contact_separations(ratio):
    """
    The four contact separations of an occultation.

    Returns
    -------
    dict
        ``first`` and ``fourth`` are the separations at which the occulter
        touches the solar limb from outside, both equal to ``1 + ratio``.
        ``second`` and ``third`` are where the occulter is fully inside the
        disk, both equal to ``1 - ratio``. For a ratio above one the occulter
        is larger than the Sun and the inner contacts do not exist, so they
        come back as `None`.

    Notes
    -----
    Between first and second contact the occulter is crossing the limb, and
    that is where the light curve is steepest. Between second and third the
    occulter is fully on the disk and the light curve only changes because
    limb darkening changes under it. The depth of that slow part is the
    measurement of limb darkening.
    """
    ratio = float(ratio)
    if ratio < 0:
        raise ValueError("ratio must be non-negative")
    inner = 1.0 - ratio if ratio <= 1.0 else None
    return {
        "first": 1.0 + ratio,
        "second": inner,
        "third": inner,
        "fourth": 1.0 + ratio,
    }


def impact_parameter_from_track(z):
    """
    Smallest separation reached by a track, taken straight from the samples.

    Parameters
    ----------
    z : array_like
        Projected separation samples.

    Returns
    -------
    float

    Notes
    -----
    This is the sampled minimum, so it is only as good as the cadence. For a
    proper impact parameter, fit a parabola to the few points around the
    minimum instead.
    """
    z = np.asarray(z, dtype=np.float64)
    if z.size == 0:
        raise ValueError("z is empty")
    return float(np.min(z))


def _require_sunpy():
    try:
        import astropy.units as u
        import sunpy.coordinates
        from sunpy.coordinates import Helioprojective, get_horizons_coord
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError(
            "the ephemeris functions need sunpy and astroquery. Install them "
            "with: pip install 'solarbatman[solar]' astroquery"
        ) from exc
    return u, sunpy.coordinates, Helioprojective, get_horizons_coord


def separation_from_ephemeris(times, body, observer="SDO"):
    """
    Projected separation of a body from the solar disk centre, in solar radii.

    Parameters
    ----------
    times : `~astropy.time.Time` or array_like
        Times to evaluate. Anything sunpy can parse as a time works.
    body : str
        Body name understood by JPL Horizons, such as ``"Moon"``,
        ``"Mercury"`` or ``"Venus"``.
    observer : str, optional
        Observer name understood by JPL Horizons. ``"SDO"`` is the Solar
        Dynamics Observatory, which is where the useful transit data comes
        from.

    Returns
    -------
    numpy.ndarray
        Projected separation in solar radii, ready to pass to
        `solarbatman.occultation_flux`.

    Notes
    -----
    The separation is computed in the helioprojective frame of the observer.
    In that frame the Sun sits at the origin by construction, so the
    helioprojective longitude and latitude of the body are already the sky
    offsets from the solar disk centre. Dividing the total offset by the
    apparent solar angular radius gives the separation in solar radii.

    This is the step that makes solar occultations work with a model written
    for exoplanets. No Keplerian orbit is fitted or assumed anywhere. The
    ephemeris supplies the geometry and the kernel supplies the photometry.
    """
    u, coords, Helioprojective, get_horizons_coord = _require_sunpy()

    observer_coord = get_horizons_coord(observer, times)
    body_coord = get_horizons_coord(body, times)

    frame = Helioprojective(observer=observer_coord, obstime=observer_coord.obstime)
    projected = body_coord.transform_to(frame)

    offset = np.hypot(
        projected.Tx.to_value(u.arcsec), projected.Ty.to_value(u.arcsec)
    )
    solar_radius = coords.sun.angular_radius(observer_coord.obstime).to_value(u.arcsec)
    return offset / solar_radius


def radius_ratio_from_ephemeris(times, body, observer="SDO"):
    """
    Angular radius ratio of a body against the Sun, from an ephemeris.

    Parameters
    ----------
    times : `~astropy.time.Time` or array_like
        Times to evaluate.
    body : str
        Body name, which must also appear in `BODY_RADII_KM` so that its
        physical radius is known.
    observer : str, optional
        Observer name understood by JPL Horizons.

    Returns
    -------
    numpy.ndarray
        Radius ratio at each time.

    Notes
    -----
    The ratio changes through an occultation, because the observer to body
    distance changes. For the Moon seen from SDO the change across a single
    crossing is small but not negligible, so use the value at mid transit
    rather than at the start.
    """
    u, coords, _, get_horizons_coord = _require_sunpy()

    observer_coord = get_horizons_coord(observer, times)
    body_coord = get_horizons_coord(body, times)

    separation = observer_coord.separation_3d(body_coord).to_value(u.km)
    solar_distance = observer_coord.radius.to_value(u.km)
    return angular_radius_ratio(body, separation, solar_distance)
