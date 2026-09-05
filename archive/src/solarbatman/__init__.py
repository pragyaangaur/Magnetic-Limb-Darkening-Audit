"""Occultation light curve modelling for the Sun, built on the batman transit model.

The batman package models a limb darkened stellar disk occulted by an opaque
circle. The Sun is a limb darkened disk, and the Moon, Mercury and Venus are
opaque circles that really do cross it. This package supplies the solar side of
that problem: the geometry, the disk integrated photometry and the limb
darkening measurement, and it calls the batman kernels directly so that no
Keplerian orbit has to be invented.
"""

from .fitting import fit_occultation, occultation_depth
from .geometry import (
    BODY_RADII_KM,
    SOLAR_RADIUS_KM,
    angular_radius_ratio,
    contact_separations,
    linear_track,
)
from .kernel import LIMB_DARKENING_LAWS, choose_step_size, n_coefficients, occultation_flux
from .limbdark import SOLAR_REFERENCE, fit_intensity_profile, intensity
from .simulate import add_spot, brute_force_flux, disk_image, occulted_disk_image

__version__ = "0.1.0"

__all__ = [
    "BODY_RADII_KM",
    "LIMB_DARKENING_LAWS",
    "SOLAR_RADIUS_KM",
    "SOLAR_REFERENCE",
    "angular_radius_ratio",
    "contact_separations",
    "fit_occultation",
    "linear_track",
    "occultation_depth",
    "add_spot",
    "brute_force_flux",
    "choose_step_size",
    "disk_image",
    "fit_intensity_profile",
    "intensity",
    "n_coefficients",
    "occultation_flux",
    "occulted_disk_image",
    "__version__",
]
